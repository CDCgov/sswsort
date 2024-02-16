#![allow(unused_variables)]
#![feature(iterator_try_collect)]
use serde_derive::Deserialize;
use std::path::PathBuf;
use std::{
    collections::HashSet,
    fmt,
    fs::read_to_string,
    io::{Error, ErrorKind},
};
use toml::from_str;
use zoe::{alignment::sw::sw_score, data::fasta::FastaNT, prelude::*};

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
    Classification { taxon: String, score: i32, strand: Strand },
    Chimeric { taxa: Vec<String>, strand: Strand },
    UnusuallyLong { taxon: String, score: i32, strand: Strand },
}

impl fmt::Display for ClassificationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ClassificationResult::Unrecognizable => write!(f, "was UNRECOGNIZABLE!"),
            ClassificationResult::Classification { taxon, score, strand } => {
                write!(f, "{}\t{}\t{}", taxon, score, strand)
            }
            ClassificationResult::Chimeric { taxa, strand } => {
                write!(f, "Chimeric: {}\t{}", taxa.join("+"), strand)
            }
            ClassificationResult::UnusuallyLong { taxon, score, strand } => {
                write!(f, "Unusually Long: {}\t{}\t{}", taxon, score, strand)
            }
        }
    }
}

#[inline]
fn calculate_alignment_score<'a>(
    reference: &'a FastaNTAnnot, query_pos: &'a FastaNT, query_neg: &'a Nucleotides,
) -> (i32, &'a str, Strand) {
    let score_pos: i32 = get_alignment_score(reference.sequence.as_bytes(), query_pos.sequence.as_bytes());
    let score_neg: i32 = get_alignment_score(reference.sequence.as_bytes(), query_neg.as_bytes());
    if score_neg > score_pos {
        (score_neg, reference.taxon.as_str(), Strand::Minus)
    } else {
        (score_pos, reference.taxon.as_str(), Strand::Plus)
    }
}

pub fn classify(query_pos: &FastaNT, params: &SSWSortArgs) -> ClassificationResult {
    let hit_threshold = 800;
    if query_pos.sequence.len() < params.length_minimum {
        return ClassificationResult::Unrecognizable;
    }

    let query_neg: Nucleotides = query_pos.sequence.reverse_complement();
    let taxa = params
        .references
        .iter()
        .map(|reference| calculate_alignment_score(reference, query_pos, &query_neg))
        .collect::<Vec<_>>();

    if params.detect_chimera {
        classify_chimeric(taxa)
    } else {
        classify_non_chimeric(taxa, params, query_pos)
    }
}

#[inline]
fn classify_chimeric(taxa: Vec<(i32, &str, Strand)>) -> ClassificationResult {
    let mut seen_strings = HashSet::new();
    let mut result = Vec::new();

    for (num, string, strand) in taxa {
        if seen_strings.insert(string) && num >= 800 {
            result.push((num, string, strand));
        }
    }

    result.sort_by_key(|&(score, _, _)| score);

    match result.len() {
        0 => ClassificationResult::Unrecognizable,
        1 => ClassificationResult::Classification {
            score:  result[0].0,
            taxon:  result[0].1.to_owned(),
            strand: result[0].2,
        },
        _ => {
            let taxa: Vec<String> = result.iter().map(|(_, taxon, _)| taxon.to_string()).collect();
            ClassificationResult::Chimeric {
                taxa,
                strand: Strand::Unknown,
            }
        }
    }
}

#[inline]
fn classify_non_chimeric(taxa: Vec<(i32, &str, Strand)>, params: &SSWSortArgs, query_pos: &FastaNT) -> ClassificationResult {
    if let Some((best_score, best_taxon, best_strand)) = taxa.iter().max_by_key(|&&(score, _, _)| score) {
        let taxon = best_taxon.to_string();
        let strand = *best_strand;

        let norm_score: f32 = *best_score as f32 / query_pos.sequence.len() as f32;

        if norm_score < params.norm_score_minimum && *best_score < params.score_minimum {
            ClassificationResult::Unrecognizable
        } else {
            ClassificationResult::Classification {
                taxon,
                score: *best_score,
                strand,
            }
        }
    } else {
        ClassificationResult::Unrecognizable
    }
}

// SSS will fix this
pub fn get_alignment_score(reference: &[u8], query: &[u8]) -> i32 {
    let (_, _, score) = sw_score(reference, query, GAP_OPEN, GAP_EXTEND, &SM);
    score
}

// Weights are fixed presently
const GAP_OPEN: i32 = -10;
const GAP_EXTEND: i32 = -1;
const MISMATCH: i32 = -5;
const MATCH: i32 = 2;

// Push up to Zoe and make 32 x 32
const SM: [[i32; 256]; 256] = {
    let mut bytes = [[0i32; 256]; 256];
    let bases = b"AGCT";

    let mut i = 0;
    while i < bases.len() {
        let mut j = 0;
        while j < bases.len() {
            bytes[bases[i] as usize][bases[j] as usize] = if i == j { MATCH } else { MISMATCH };
            j += 1;
        }
        i += 1;
    }
    bytes
};

pub struct FastaNTAnnot {
    pub name:     String,
    pub sequence: Nucleotides,
    pub taxon:    String,
}

pub struct SSWSortArgs {
    pub name:               String,
    pub references:         Vec<FastaNTAnnot>,
    pub norm_score_minimum: f32,
    pub score_minimum:      i32,
    pub length_minimum:     usize,
    pub detect_chimera:     bool,
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
    pub score_minimum:       i32,
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

        let references = FastaReader::from_filename(reference_sequences)?
            .filter_map(|f| f.ok())
            .map(|fa| {
                if let Some((id, taxon)) = fa.get_id_taxon() {
                    Ok(FastaNTAnnot {
                        name:     id.to_string(),
                        sequence: fa.sequence.filter_to_dna(),
                        taxon:    taxon.to_string(),
                    })
                } else {
                    Err(Error::new(
                        ErrorKind::InvalidData,
                        format!("No Reference Taxon: {id}", id = fa.name),
                    ))
                }
            })
            .try_collect()?;

        Ok(SSWSortArgs {
            name,
            references,
            norm_score_minimum,
            score_minimum,
            length_minimum,
            detect_chimera,
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
    let configs: TomlConfig = from_str(&raw_toml)?;

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
