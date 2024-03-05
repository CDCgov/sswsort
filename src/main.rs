#![allow(dead_code, clippy::single_element_loop, unused_variables)]

//use ahash::RandomState;
use clap::Parser;
use rayon::iter::ParallelBridge;
use rayon::prelude::ParallelIterator;
use sswsort::*;
use std::path::PathBuf;
use zoe::prelude::*;

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
    /// Run in simultaneous multi-threaded mode.
    #[arg(short = 'T', long)]
    threads: bool,
    // TO-DO: Implement a grid mode that takes the total processes and the process number
}

fn main() {
    let args = ClassifierArgs::parse();

    let query_reader = FastaReader::from_filename(args.fasta_file)
        .unwrap_or_die("Cannot open FASTA file!")
        .filter_map(|f| f.ok())
        .map(|seq| seq.filter_to_dna());

    let params = get_sswsort_module_args(&args.module).unwrap_or_die("Failed to load module data!");

    let results: Vec<(ClassificationResult, String, usize)> = if args.threads {
        query_reader
            .par_bridge()
            .map(|query| (classify(&query, &params), query.name, query.sequence.len()))
            .collect()
    } else {
        query_reader
            .map(|query| (classify(&query, &params), query.name, query.sequence.len()))
            .collect()
    };

    for (classification, query_name, seq_length) in results.into_iter() {
        println!(
            "{PROGRAM_VERSION}\t{module_name}\t{query_name}\t{classification}",
            module_name = args.module,
        );
    }
}
