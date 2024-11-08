use foldhash::{HashMap, HashMapExt};
use serde_derive::Deserialize;
use std::{cmp::max, collections::BTreeSet, fmt, fs::read_to_string, io::Error, io::ErrorKind, path::PathBuf};
use toml::from_str;
use zoe::{
    alignment::LocalProfiles,
    data::{WeightMatrix, fasta::FastaNT, fasta::FastaNTAnnot},
    prelude::*,
};

// Weights are fixed presently
const GAP_OPEN: i8 = -10;
const GAP_EXTEND: i8 = -1;
const MISMATCH: i8 = -5;
const MATCH: i8 = 2;
const MATRIX: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(MATCH, MISMATCH, Some(b'N'));

#[derive(Hash, Eq, PartialEq, Clone, Copy)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
            Strand::Unknown => write!(f, "*"),
        }
    }
}

#[non_exhaustive]
pub enum ClassificationResult {
    Unrecognizable,
    Classification {
        taxon:      String,
        best_score: u64,
        strand:     Strand,
    },
    Chimeric {
        taxa:   Vec<String>,
        strand: Strand,
    },
    UnusuallyLong {
        taxon:      String,
        best_score: u64,
        strand:     Strand,
    },
}

impl fmt::Display for ClassificationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ClassificationResult::Unrecognizable => write!(f, "was UNRECOGNIZABLE!"),
            ClassificationResult::Classification {
                taxon,
                best_score,
                strand,
            } => {
                write!(f, "{taxon}\t{best_score}\t{strand}")
            }
            ClassificationResult::Chimeric { taxa, strand } => {
                write!(f, "Chimeric: {}\t{}", taxa.join("+"), strand)
            }
            ClassificationResult::UnusuallyLong {
                taxon,
                best_score,
                strand,
            } => {
                write!(f, "Unusually Long: {taxon}\t{best_score}\t{strand}")
            }
        }
    }
}

#[inline]
fn calculate_alignment_score<'a>(
    reference: &'a (FastaNTAnnot, Nucleotides), query_profile: &LocalProfiles<64, 32, 16, 8, 5>,
) -> (u64, &'a str, Strand) {
    let (reference, reference_seq_revcomp) = reference;

    let score_pos = query_profile.smith_waterman_score_from_i8(&reference.sequence);
    let score_neg = query_profile.smith_waterman_score_from_i8(reference_seq_revcomp);

    match (score_pos, score_neg) {
        (Some(pos), Some(neg)) if neg > pos => (neg, reference.taxon.as_str(), Strand::Minus),
        (Some(pos), Some(_)) => (pos, reference.taxon.as_str(), Strand::Plus),
        _ => panic!("Alignment failed!"),
    }
}

pub fn classify(query: &FastaNT, params: &SSWSortArgs) -> ClassificationResult {
    if query.sequence.len() < params.length_minimum {
        return ClassificationResult::Unrecognizable;
    }

    let Ok(query_profile) = LocalProfiles::new_with_w512(&query.sequence, &MATRIX, GAP_OPEN, GAP_EXTEND) else {
        eprintln!("WARNING: '{id}' was empty or had bad scoring weights.", id = query.name);
        return ClassificationResult::Unrecognizable;
    };

    let taxa = params
        .references
        .iter()
        .map(|reference| calculate_alignment_score(reference, &query_profile));

    // Does not allocate unless data is pushed
    let mut valid_chimera_taxa = BTreeSet::new();
    let hit_threshold = 800;

    let best = if params.detect_chimera {
        let opt = taxa
            .inspect(|&(score, taxon, _)| {
                if score >= hit_threshold {
                    valid_chimera_taxa.insert(taxon);
                }
            })
            .max_by_key(|&(score, _, _)| score);

        if valid_chimera_taxa.len() > 1 {
            let taxa = valid_chimera_taxa.into_iter().rev().map(str::to_string).collect();
            return ClassificationResult::Chimeric {
                taxa,
                strand: Strand::Unknown,
            };
        }
        opt
    } else {
        taxa.max_by_key(|&(score, _, _)| score)
    };

    if let Some((best_score, taxon, strand)) = best
        && (best_score >= params.score_minimum
            || best_score as f32 / query.sequence.len() as f32 >= params.norm_score_minimum)
    {
        let taxon: String = taxon.to_string();

        if let Some(&max_length) = params.length_by_annot.get(&taxon)
            && query.sequence.len() >= max_length
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
        ClassificationResult::Unrecognizable
    }
}

pub struct SSWSortArgs {
    pub name:               String,
    pub references:         Vec<(FastaNTAnnot, Nucleotides)>,
    pub norm_score_minimum: f32,
    pub score_minimum:      u64,
    pub length_minimum:     usize,
    pub detect_chimera:     bool,
    pub length_by_annot:    HashMap<String, usize>,
}

#[derive(Deserialize, Debug)]
pub struct TomlConfig {
    // Can also place other top-level items in here
    pub classification_module: Vec<ModuleParameters>,
}

#[derive(Deserialize, Debug, Default)]
pub struct ModuleParameters {
    pub name:                String,
    pub alternative_names:   Vec<String>,
    pub norm_score_minimum:  f32,
    pub score_minimum:       u64,
    pub length_minimum:      usize,
    pub reference_sequences: PathBuf,
    pub detect_chimera:      bool,
}

impl TryFrom<ModuleParameters> for SSWSortArgs {
    type Error = std::io::Error;
    fn try_from(params: ModuleParameters) -> Result<Self, Self::Error> {
        let ModuleParameters {
            name,
            alternative_names: _,
            norm_score_minimum,
            score_minimum,
            length_minimum,
            reference_sequences,
            detect_chimera,
        } = params;

        let reference_annots = FastaReader::from_filename(reference_sequences)?
            .filter_map(|f| f.ok())
            .map(FastaNTAnnot::try_from)
            .collect::<Result<Vec<FastaNTAnnot>, Error>>()?;

        let length_factor = 2;
        let mut length_by_annot = HashMap::new();

        // Create references vector and populate length_by_annot in single pass
        let references = reference_annots
            .into_iter()
            .map(|fa| {
                let this_len = fa.sequence.len() * length_factor;
                length_by_annot
                    .entry(fa.taxon.clone())
                    .and_modify(|len| *len = max(*len, this_len))
                    .or_insert(this_len);

                let reference_sequence_revcomp = fa.sequence.to_reverse_complement();

                (fa, reference_sequence_revcomp)
            })
            .collect();

        Ok(SSWSortArgs {
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

pub fn get_sswsort_module_args(preset_module: &str) -> Result<SSWSortArgs, Error> {
    let suffix = "sswsort_res/config.toml";
    let preset_module = preset_module.to_ascii_lowercase();
    let mut exe_path: PathBuf = std::env::current_exe()?;
    exe_path.pop();

    // Check same directory
    let mut toml_path = exe_path.join(suffix);
    if !toml_path.exists() {
        // otherwise check grandparent
        exe_path.pop();
        exe_path.pop();
        toml_path = exe_path.join(suffix);
    }

    let raw_toml = read_to_string(toml_path)?;
    let configs: TomlConfig = from_str(&raw_toml).map_err(|e| std::io::Error::new(ErrorKind::InvalidData, e))?;

    let selected_config = configs
        .classification_module
        .into_iter()
        .find_map(|mut config| {
            if preset_module == config.name.to_ascii_lowercase()
                || config
                    .alternative_names
                    .iter()
                    .any(|alt| alt.to_ascii_lowercase() == preset_module)
            {
                config.reference_sequences = exe_path.join("sswsort_res/").join(config.reference_sequences);

                #[cfg(debug_assertions)]
                println!("{config:?}");

                Some(config)
            } else {
                None
            }
        })
        .unwrap_or_else(|| panic!("Could not find a configuration for: {preset_module}"));

    // Update reference sequence with full path from above
    let params: SSWSortArgs = selected_config.try_into()?;

    Ok(params)
}
