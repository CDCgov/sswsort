use foldhash::{HashMap, HashMapExt};
use serde_derive::Deserialize;
use std::{
    cmp::max,
    collections::{BTreeMap, BTreeSet},
    fmt,
    fs::read_to_string,
    io::Error,
    path::{Path, PathBuf},
};
use zoe::{
    alignment::{LocalProfiles, ProfileError, ProfileSets},
    data::{
        WeightMatrix,
        err::ResultWithErrorContext,
        fasta::{FastaNT, FastaNTAnnot},
    },
    iter_utils::ProcessResultsExt,
    prelude::*,
};

// Weights are fixed presently
const GAP_OPEN: i8 = -10;
const GAP_EXTEND: i8 = -1;
const MISMATCH: i8 = -5;
const MATCH: i8 = 2;
const MATRIX: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(MATCH, MISMATCH, Some(b'N'));

/// The strand of the reference to which the query best aligned.
#[derive(Hash, Eq, PartialEq, Clone, Copy, Debug)]
pub enum Strand {
    /// The query aligned to the forward reference sequence.
    Plus,
    /// The query aligned to the reverse complement of the reference sequence.
    Minus,
}

/// Removes repeating `N` residues at the beginning and end of the sequence.
///
/// Since `N` residues are treated as 0-cost mismatches, and local alignment is
/// being done, this will speed up alignment time without affecting the score.
fn trim_n(sequence: &Nucleotides) -> NucleotidesView<'_> {
    let mut iter = sequence.iter();
    let left = iter.by_ref().take_while(|b| **b == b'N').count();
    let right = iter.rev().take_while(|b| **b == b'N').count();
    sequence.slice(left..sequence.len() - right)
}

impl fmt::Display for Strand {
    /// Displays the [`Strand`], using "+" and "-" for positive and negative
    /// strand.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}

/// The various classifications that could be given by SSWSORT.
///
/// ## Parameters
///
/// The lifetime parameter is for the stored taxa, if applicable.
#[derive(PartialEq, Debug)]
pub enum ClassificationResult<'a> {
    /// A classification that is not recognizable or not applicable.
    ///
    /// This is created by [`ClassificationResult::none`], in which case the
    /// internal score field is not populated. It can also be constructed with a
    /// `best_score`, which is used when a score was generated but did not pass
    /// the `score_minimum` or `norm_score_minimum` requirements.
    Unrecognizable { best_score: Option<u32> },
    /// A successful classification with a single taxon.
    Classification {
        taxon:      &'a str,
        best_score: u32,
        strand:     Strand,
    },
    /// A chimeric classification that could belong to any of the stored taxa.
    ///
    /// This is only used when `detect_chimera` is true. A classification is
    /// chimeric if multiple alignments against distinct taxa score at least
    /// 800.
    Chimeric { taxa: Vec<&'a str> },
    /// A classification issued when the query is at least twice as long as the
    /// longest reference for the given taxon.
    UnusuallyLong {
        taxon:      &'a str,
        best_score: u32,
        strand:     Strand,
    },
    /// A classification issued when the classified scores are equal to each
    /// other for distinct taxa.
    Unresolvable { taxa: Vec<&'a str>, best_score: u32 },
}

impl<'a> ClassificationResult<'a> {
    /// Takes the top score and taxon (or taxa) for a given query and returns a
    /// [`ClassificationResult`].
    fn classify_top(
        best_score: u32, taxa: Vec<&'a str>, strand: Strand, module: &'a SSWSortModule, query_len: usize,
    ) -> ClassificationResult<'a> {
        if best_score >= module.score_minimum || best_score as f32 / query_len as f32 >= module.norm_score_minimum {
            if let &[taxon] = taxa.as_slice() {
                if let Some(&max_length) = module.length_by_annot.get(taxon)
                    && query_len >= max_length
                {
                    ClassificationResult::UnusuallyLong {
                        taxon,
                        best_score,
                        strand,
                    }
                } else {
                    ClassificationResult::Classification {
                        taxon,
                        best_score,
                        strand,
                    }
                }
            } else {
                ClassificationResult::Unresolvable { taxa, best_score }
            }
        } else {
            ClassificationResult::Unrecognizable {
                best_score: Some(best_score),
            }
        }
    }

    /// Constructs a new [`Unrecognizable`] classification without a score.
    ///
    /// This is used when an error occurs during alignment, the sequence length
    /// is shorter than `length_minimum`, or there is no applicable
    /// classification data (e.g., for secondary classifications).
    ///
    /// [`Unrecognizable`]: ClassificationResult::Unrecognizable
    pub fn none() -> Self {
        ClassificationResult::Unrecognizable { best_score: None }
    }
}

/// Takes a reference and query profile to calculate a Smith Waterman alignment
/// score. `None` is returned if no alignment is found (or if the alignment
/// overflows `i32`).
///
/// The `reference` argument should contain a tuple with the original record and
/// the pre-computed reverse complement.
#[inline]
fn calculate_alignment_score<'a>(
    reference: &'a (FastaNTAnnot, Nucleotides), query_profile: &LocalProfiles<64, 32, 16, 5>,
) -> Option<(u32, &'a str, Strand)> {
    let (reference, reference_seq_revcomp) = reference;

    let score_pos = query_profile.sw_score_from_i8(&reference.sequence).get();
    let score_neg = query_profile.sw_score_from_i8(reference_seq_revcomp).get();

    match (score_pos, score_neg) {
        (Some(pos), Some(neg)) => {
            if pos >= neg {
                Some((pos, &reference.taxon, Strand::Plus))
            } else {
                Some((neg, &reference.taxon, Strand::Minus))
            }
        }
        (Some(pos), None) => Some((pos, &reference.taxon, Strand::Plus)),
        (None, Some(neg)) => Some((neg, &reference.taxon, Strand::Minus)),
        (None, None) => None,
    }
}

/// Takes an iterator of `(score, taxon, Strand)` and finds the top result,
/// allowing for ties, and ensuring no repeated taxa.
///
/// The taxa are returned in alphabetical order.
///
/// `None` is returned if the input iterator is empty. If there is a tie
/// (multiple returned taxa), then the returned [`Strand`] corresponds to the
/// first taxon in alphabetical order.
fn top_with_ties<'a, I>(iter: I) -> Option<(u32, Vec<&'a str>, Strand)>
where
    I: IntoIterator<Item = (u32, &'a str, Strand)>, {
    let mut top_score = None;
    let mut top_taxa_strand = BTreeMap::new();

    for (score, taxon, strand) in iter {
        match top_score {
            None => {
                top_score = Some(score);
                top_taxa_strand.insert(taxon, strand);
            }
            Some(tscore) if score > tscore => {
                top_score = Some(score);
                top_taxa_strand.clear();
                top_taxa_strand.insert(taxon, strand);
            }
            Some(tscore) if score == tscore => {
                top_taxa_strand.insert(taxon, strand);
            }
            _ => {}
        }
    }

    let strand = *top_taxa_strand.values().next()?;
    let score = top_score?;
    let taxa = top_taxa_strand.keys().copied().collect();

    Some((score, taxa, strand))
}

/// Takes an iterator of `(score, taxon, Strand)` and finds the top two results,
/// allowing for ties, and ensuring no repeated taxa in the top two.
fn top_two_with_ties<'a, I>(iter: I) -> [Option<(u32, Vec<&'a str>, Strand)>; 2]
where
    I: IntoIterator<Item = (u32, &'a str, Strand)>, {
    let mut top_score = None;
    let mut second_score = None;

    let mut top_taxa_strand = BTreeMap::new();
    let mut second_taxa_strand = BTreeMap::new();

    for (score, taxon, strand) in iter {
        match top_score {
            None => {
                top_score = Some(score);
                top_taxa_strand.insert(taxon, strand);
            }
            Some(tscore) if score > tscore => {
                second_score = top_score;
                second_taxa_strand = std::mem::take(&mut top_taxa_strand);

                top_score = Some(score);
                top_taxa_strand.insert(taxon, strand);
            }
            Some(tscore) if score == tscore => {
                top_taxa_strand.insert(taxon, strand);
                top_score = Some(tscore);
            }
            // Ensure that new secondary classifications do not match an
            // existing primary classification
            Some(tscore) if score < tscore && !top_taxa_strand.contains_key(taxon) => match second_score {
                None => {
                    second_score = Some(score);
                    second_taxa_strand.insert(taxon, strand);
                }
                Some(sscore) if score > sscore => {
                    second_score = Some(score);
                    second_taxa_strand.clear();
                    second_taxa_strand.insert(taxon, strand);
                }
                Some(sscore) if score == sscore => {
                    second_taxa_strand.insert(taxon, strand);
                    second_score = Some(sscore);
                }
                _ => {}
            },
            _ => {}
        }
    }

    // Ensure that if the primary classification is resolvable, then the
    // secondary classification does not include the primary taxon
    if top_taxa_strand.len() == 1
        && let Some((taxon, _)) = top_taxa_strand.iter().next()
    {
        second_taxa_strand.remove(taxon);
    }

    let top_classification = (top_score, top_taxa_strand);
    let second_classification = (second_score, second_taxa_strand);

    [top_classification, second_classification].map(|(score_opt, taxa_strand)| {
        // Get first strand, or return None
        let strand = *taxa_strand.values().next()?;

        // Get the score, or return None
        let score = score_opt?;

        // Collect all the taxa in sorted order
        let taxa = taxa_strand.keys().copied().collect();

        Some((score, taxa, strand))
    })
}

/// Takes the top two `(score, taxa, Strand)` and transforms them into
/// [`ClassificationResult`].
///
/// If either is `None`, then [`ClassificationResult::none`] is used.
fn classify_top_two_helper<'a>(
    top_2: [Option<(u32, Vec<&'a str>, Strand)>; 2], module: &'a SSWSortModule, query_len: usize,
) -> [ClassificationResult<'a>; 2] {
    let mut classified = top_2.map(|info| {
        if let Some((best_score, taxa, strand)) = info {
            ClassificationResult::classify_top(best_score, taxa, strand, module, query_len)
        } else {
            ClassificationResult::none()
        }
    });

    // ensures that if first classification is unresolvable, we don't provide a
    // second classification, similar to chimera
    if matches!(classified[0], ClassificationResult::Unresolvable { .. }) {
        classified[1] = ClassificationResult::none();
    }

    classified
}

impl SSWSortModule {
    /// Returns the top classification for a single sequence.
    ///
    /// ## Errors
    ///
    /// Any errors when building the alignment profile are propagated.
    pub fn classify(&self, query: &FastaNT) -> Result<ClassificationResult<'_>, ProfileError> {
        let FastaNT { sequence, .. } = query;
        let sequence = trim_n(sequence);

        if sequence.len() < self.length_minimum {
            return Ok(ClassificationResult::none());
        }

        let query_profile = LocalProfiles::new_with_w512(&sequence, &MATRIX, GAP_OPEN, GAP_EXTEND)?;

        let taxa = self
            .references
            .iter()
            .filter_map(|reference| calculate_alignment_score(reference, &query_profile));

        let mut valid_chimera_taxa = BTreeSet::new();
        let hit_threshold = 800;

        let best = if self.detect_chimera {
            let taxa_iter = taxa.inspect(|&(score, taxon, _)| {
                if score >= hit_threshold {
                    valid_chimera_taxa.insert(taxon);
                }
            });

            let opt = top_with_ties(taxa_iter);

            if valid_chimera_taxa.len() > 1 {
                let taxa = valid_chimera_taxa.into_iter().collect();
                return Ok(ClassificationResult::Chimeric { taxa });
            }
            opt
        } else {
            top_with_ties(taxa)
        };

        if let Some((best_score, taxa, strand)) = best {
            Ok(ClassificationResult::classify_top(
                best_score,
                taxa,
                strand,
                self,
                sequence.len(),
            ))
        } else {
            Ok(ClassificationResult::none())
        }
    }

    /// Returns the top two classifications for a single sequence.
    ///
    /// ## Errors
    ///
    /// Any errors when building the alignment profile are propagated.
    pub fn classify_top_two(&self, query: &FastaNT) -> Result<[ClassificationResult<'_>; 2], ProfileError> {
        let FastaNT { sequence, .. } = query;
        let sequence = trim_n(sequence);

        if sequence.len() < self.length_minimum {
            return Ok([ClassificationResult::none(), ClassificationResult::none()]);
        }

        let query_profile = LocalProfiles::new_with_w512(&sequence, &MATRIX, GAP_OPEN, GAP_EXTEND)?;

        let taxa = self
            .references
            .iter()
            .filter_map(|reference| calculate_alignment_score(reference, &query_profile));

        // Does not allocate unless data is pushed
        let mut valid_chimera_taxa = BTreeSet::new();
        let hit_threshold = 800;

        let top_2 = if self.detect_chimera {
            let taxa_iter = taxa.inspect(|&(score, taxon, _)| {
                if score >= hit_threshold {
                    valid_chimera_taxa.insert(taxon);
                }
            });

            let top_2 = top_two_with_ties(taxa_iter);
            if valid_chimera_taxa.len() > 1 {
                let taxa = valid_chimera_taxa.into_iter().collect();
                return Ok([ClassificationResult::Chimeric { taxa }, ClassificationResult::none()]);
            }
            top_2
        } else {
            top_two_with_ties(taxa)
        };

        Ok(classify_top_two_helper(top_2, self, sequence.len()))
    }
}

/// A module containing all information needed for performing classification
/// with `sswsort`.
pub struct SSWSortModule {
    /// The name of the module (e.g., `flu`, `cov`, `spike`, `rsv`).
    pub name:           String,
    /// The references along with their reverse complements.
    references:         Vec<(FastaNTAnnot, Nucleotides)>,
    /// The minimum normalized score for deciding unrecognizability.
    norm_score_minimum: f32,
    /// The minimum absolute score for deciding unrecognizability.
    score_minimum:      u32,
    /// The minimum required length for deciding unrecognizability.
    length_minimum:     usize,
    /// Whether or not to detect chimeric sequences.
    detect_chimera:     bool,
    /// A map from taxa to its maximum reference length.
    length_by_annot:    HashMap<String, usize>,
}

/// For reading config options from `sswsort_res/config.toml`
#[derive(Deserialize, Debug)]
pub struct TomlConfig {
    classification_module: Vec<ModuleParameters>,
}

/// Configuration parameters for a module.
///
/// This can be parsed from the `config.toml` file (with
/// [`ModuleParameters::load`]), or it can be constructed manually.
///
/// To perform classification using [`ModuleParameters`], first turn it into a
/// [`SSWSortModule`] (which has all reference sequence loaded into memory) via
/// [`SSWSortModule::new`].
///
/// ## Description of Thresholds
///
/// If the query is less than `minimum_length` in length after trimming leading
/// and trailing `N`s, then [`ClassificationResult::Unrecognizable`] is returned
/// by classification functions. Or, if the optimal alignment yields a
/// normalized score less than `norm_score_minimum` _and_ an absolute score less
/// than `score_minimum`, then [`ClassificationResult::Unrecognizable`] is also
/// returned.
#[derive(Deserialize, Debug, Default)]
pub struct ModuleParameters {
    /// The name of the module (e.g., `flu`, `cov`, `spike`, `rsv`).
    pub name:                String,
    /// An optional version for the module to aid in version control.
    pub version:             Option<String>,
    /// Any alternative names that can be used to refer to the module.
    pub alternative_names:   Vec<String>,
    /// The minimum normalized score for deciding unrecognizability.
    pub norm_score_minimum:  f32,
    /// The minimum absolute score for deciding unrecognizability.
    pub score_minimum:       u32,
    /// The minimum required length for deciding unrecognizability.
    pub length_minimum:      usize,
    /// A path to the FASTA file of reference sequences.
    pub reference_sequences: PathBuf,
    /// Whether or not to detect chimeric sequences.
    pub detect_chimera:      bool,
}

impl SSWSortModule {
    /// Creates a new [`SSWSortModule`] (with references loaded and ready for
    /// classification).
    ///
    /// ## Errors
    ///
    /// Any errors loading the references are propagated.
    pub fn new(params: ModuleParameters) -> std::io::Result<Self> {
        let ModuleParameters {
            mut name,
            version,
            alternative_names: _,
            norm_score_minimum,
            score_minimum,
            length_minimum,
            reference_sequences,
            detect_chimera,
        } = params;

        let length_factor = 2;
        let mut length_by_annot = HashMap::new();

        // Create references vector and populate length_by_annot in single pass
        let references = FastaReader::from_path(reference_sequences)?
            .map(|res| res.and_then(FastaNTAnnot::try_from))
            .process_results(|iter| {
                iter.map(|fa| {
                    let this_len = fa.sequence.len() * length_factor;
                    length_by_annot
                        .entry(fa.taxon.clone())
                        .and_modify(|len| *len = max(*len, this_len))
                        .or_insert(this_len);

                    let reference_sequence_revcomp = fa.sequence.to_reverse_complement();

                    (fa, reference_sequence_revcomp)
                })
                .collect::<Vec<_>>()
            })?;

        if let Some(v) = version {
            name.push_str(" v");
            name.push_str(&v);
        }

        Ok(SSWSortModule {
            name,
            references,
            norm_score_minimum,
            score_minimum,
            length_minimum,
            detect_chimera,
            length_by_annot,
        })
    }
}

/// A convenience function for loading a module from the TOML file, which is
/// automatically-located.
///
/// The TOML file is searched for in the same directory as the executable, as
/// well as the grandparent directory.
///
/// ## Errors
///
/// - The path to the executable must be able to be determined
/// - The TOML file must successfully parse and contain the specified module
/// - Loading the references must succeed
pub fn get_sswsort_module(preset_module: &str) -> std::io::Result<SSWSortModule> {
    let suffix = "sswsort_res/config.toml";
    let mut exe_path = std::env::current_exe().with_context("Failed to find the path to the executable")?;
    exe_path.pop();

    // Check same directory
    let mut toml_path = exe_path.join(suffix);
    if !toml_path.exists() {
        // otherwise check grandparent
        exe_path.pop();
        exe_path.pop();
        toml_path = exe_path.join(suffix);
    }

    let params =
        ModuleParameters::load(&toml_path, preset_module).with_path_context("Failed to parse TOML file", toml_path)?;

    SSWSortModule::new(params)
}

impl ModuleParameters {
    /// Loads the `preset_module` from the TOML file.
    ///
    /// ## Errors
    ///
    /// Parsing errors are propogated with type context. If the module is not
    /// found, an error is also returned with the module name included.
    pub fn load(toml: impl AsRef<Path>, preset_module: &str) -> Result<Self, Error> {
        let raw_toml = read_to_string(&toml)?;

        let config = toml::from_str::<TomlConfig>(&raw_toml).with_type_context::<TomlConfig>()?;

        let sswsort_res = toml.as_ref().parent().ok_or(std::io::Error::other(
            "The TOML path must be contained in a folder containing the reference files",
        ))?;

        config
            .classification_module
            .into_iter()
            .find_map(|mut config| {
                if preset_module.eq_ignore_ascii_case(&config.name)
                    || config
                        .alternative_names
                        .iter()
                        .any(|alt| alt.eq_ignore_ascii_case(preset_module))
                {
                    config.reference_sequences = sswsort_res.join(config.reference_sequences);

                    Some(config)
                } else {
                    None
                }
            })
            .ok_or_else(|| {
                std::io::Error::other(format!(
                    "Could not find a configuration for: {preset_module}",
                    preset_module = preset_module.to_ascii_lowercase()
                ))
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_top_two_with_ties() {
        let data: Vec<(u32, &str, Strand)> = vec![
            (105, "A", Strand::Plus),
            (108, "B", Strand::Plus),
            (108, "C", Strand::Minus),
            (108, "D", Strand::Plus),
            (110, "D", Strand::Plus),
            (105, "E", Strand::Plus),
        ];

        let params = ModuleParameters {
            name:                "flu".to_owned(),
            version:             Some("2".to_owned()),
            alternative_names:   Vec::new(),
            norm_score_minimum:  1.0,
            score_minimum:       100,
            length_minimum:      25,
            reference_sequences: PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../sswsort_res/flu.fasta"),
            detect_chimera:      true,
        };
        let module = SSWSortModule::new(params).unwrap();

        let top_2 = top_two_with_ties(data);
        let classifications = classify_top_two_helper(top_2, &module, 50);

        assert_eq!(
            classifications[0],
            ClassificationResult::Classification {
                taxon:      "D",
                best_score: 110,
                strand:     Strand::Plus,
            }
        );
        assert_eq!(
            classifications[1],
            ClassificationResult::Unresolvable {
                taxa:       vec!["B", "C"],
                best_score: 108,
            }
        );
    }

    #[test]
    fn test_redundant_top_two() {
        // from William's counterexample
        let data: Vec<(u32, &str, Strand)> = vec![
            (110, "C", Strand::Plus),
            (120, "A", Strand::Plus),
            (150, "B", Strand::Minus),
            (150, "A", Strand::Plus),
            (90, "E", Strand::Plus),
        ];

        let params = ModuleParameters {
            name:                "flu".to_owned(),
            version:             Some("2".to_owned()),
            alternative_names:   Vec::new(),
            norm_score_minimum:  1.0,
            score_minimum:       100,
            length_minimum:      25,
            reference_sequences: PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../sswsort_res/flu.fasta"),
            detect_chimera:      true,
        };
        let module = SSWSortModule::new(params).unwrap();

        let top_2 = top_two_with_ties(data);
        let classifications = classify_top_two_helper(top_2, &module, 50);

        assert_eq!(
            classifications[0],
            ClassificationResult::Unresolvable {
                taxa:       vec!["A", "B"],
                best_score: 150,
            }
        );
        assert_eq!(classifications[1], ClassificationResult::none());
    }
}
