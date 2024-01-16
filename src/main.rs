#![allow(dead_code, clippy::single_element_loop, unused_variables)]

use clap::Parser;
use std::path::PathBuf;
use zoe::{data::fasta::FastaNT, prelude::*};

// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
/// Uses Smith-Waterman to classify sequences into a simple compound type.
pub struct ClassifierArgs {
    // Positional arguments
    /// Name of the classification module
    module:     String,
    /// Name of the nucleotide sequences to classify in FASTA format.
    fasta_file: PathBuf,

    // TO-DO: create flag for STDOUT. If it exists, this positional shall be ignored
    // Mutually exclusive arguments
    /// Name of the tab-separated-value file for classifier results.
    output_file: PathBuf,

    // Optional arguments
    /// Number of threads for shared-multiprocessing execution.
    #[arg(short = 'T', long, default_value_t = 1)]
    threads: u8,
    // TO-DO: Implement a grid mode that takes the total processes and the process number
}

fn main() {
    let args = ClassifierArgs::parse();

    // SSS will fix this
    let query_reader = FastaReader::from_filename(args.fasta_file)
        .unwrap_or_die("Cannot open FASTA file!")
        .filter_map(|f| f.ok())
        .map(|seq| seq.filter_to_dna());

    // Wiley:
    // 1. need module configurations
    // 2. need to ingest them
    // 3. need reference data

    for query in query_reader {
        // Will print to file, for now, print to STDOUT
        let (taxon, score, strand) = classify_query(&query);
        // TO-DO: either create a struct for output and write a Display OR just print it out
        // TO-DO: can we use arrow2 or the like to write to parquet OPTIONALLY
        // COuld interact poorly with standard out, so make that exclusive as well
        println!(
            "{PROGRAM_VERSION}\t{module_name}\t{id}\t{taxon}\t{score}\t{seq_length}\t{strand}",
            module_name = args.module,
            id = query.name,
            seq_length = query.sequence.len()
        );
    }
}

// TO-DO: need to pass reference array in some manner
// You are welcome to fix this, otherwise I will
// To-DO: classify query needs other parameters similar to the classifySAM script
// What is global or in the main program versus what makes sense to keep in local scope?
fn classify_query(query: &FastaNT) -> (String, i32, char) {
    for reference in [b"AGCT"] {
        // Classify data
        let score = get_alignment_score(reference, query.sequence.as_bytes());
    }

    ("UNRECOGNIZABLE".to_owned(), 0, '+')
}

// SSS will fix this
fn get_alignment_score(reference: &[u8], query: &[u8]) -> i32 {
    1
}
