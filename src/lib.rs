#![allow(unused_variables)]
#![feature(iterator_try_collect)]
use serde_derive::Deserialize;
use std::path::PathBuf;
use std::{
    fs::read_to_string,
    io::{Error, ErrorKind},
};
use toml::from_str;
use zoe::{alignment::sw::sw_score, data::fasta::FastaNT, prelude::*};

const UNRECOGNIZABLE: &str = "UNRECOGNIZABLE";

// Can we return invalid states differently?
pub fn classify_query(query_pos: &FastaNT, params: &SSWSortArgs) -> (String, i32, u8) {
    let query_neg = query_pos.sequence.reverse_complement();

    let (best_score, best_taxon, best_strand) = params
        .references
        .iter()
        .map(|reference| {
            let score_pos = get_alignment_score(reference.sequence.as_bytes(), query_pos.sequence.as_bytes());
            let score_neg = get_alignment_score(reference.sequence.as_bytes(), query_pos.sequence.as_bytes());
            if score_neg > score_pos {
                (score_neg, reference.taxon.as_str(), b'-')
            } else {
                (score_pos, reference.taxon.as_str(), b'+')
            }
        })
        .max_by_key(|&(score, _, _)| score)
        .unwrap_or((0, UNRECOGNIZABLE, b'?'));

    let best_taxon = best_taxon.to_string();
    if query_pos.sequence.len() < params.length_minimum || best_score < params.score_minimum {
        (best_taxon, best_score, best_strand)
    } else {
        let norm_score: f32 = if !query_pos.sequence.is_empty() {
            best_score as f32 / query_pos.sequence.len() as f32
        } else {
            0.0
        };

        if norm_score < params.norm_score_minimum {
            (UNRECOGNIZABLE.to_string(), best_score, best_strand)
        } else {
            (best_taxon, best_score, best_strand)
        }
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
