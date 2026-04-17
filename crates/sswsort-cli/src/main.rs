#![feature(int_format_into)]

use clap::Parser;
use rayon::{ThreadPoolBuilder, iter::ParallelBridge, prelude::ParallelIterator};
use sswsort::*;
use std::{io::Write, num::NonZero, path::PathBuf};
use zoe::{
    data::fasta::{FastaNT, FastaSeq},
    define_whichever,
    iter_utils::ProcessResultsExt,
    prelude::*,
};

mod app;
use app::*;

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
/// Uses striped Smith-Waterman to classify sequences into a simple compound type.
pub struct ClassifierArgs {
    // Positional arguments
    /// Name of the classification module
    module:     String,
    /// Name of the nucleotide sequences to classify in FASTA format.
    fasta_file: PathBuf,

    // Optional arguments
    /// Name of the tab-separated-value file for classifier results. If none are
    /// provided, STDOUT is used. If a directory is specified, a default
    /// filename of `sswsort_output.tsv`
    output_file: Option<PathBuf>,

    /// Number of threads to use. Defaults to number of physical cores
    /// otherwise.
    #[arg(short = 'T', long)]
    threads: Option<NonZero<usize>>,

    /// Execute as a partitioned task in a grid job, for use with: --submit-grid-job
    ///
    /// Detects the array size + task id and write out the data
    /// to a file at the prefixed location for downstream collation.
    ///
    /// An `output_file` is required and will be used with a partition suffix.
    #[arg(short = 'G', long, requires = "output_file", conflicts_with_all = ["threads", "submit_grid_job"])]
    is_grid_task: bool,

    /// Submits and blocks on a grid job of the specified array size.
    #[arg(short='S', long, value_name = "SIZE", conflicts_with_all = ["threads", "is_grid_task"])]
    submit_grid_job: Option<usize>,
}

define_whichever! {
    /// An enum representing the allowable output types.
    ///
    /// The types are not wrapped in a [`BufWriter`], so it is advised to do
    /// that before using [`AnyOutput`].
    ///
    /// [`BufWriter`]: std::io::BufWriter
    pub enum AnyOutput {
        File(std::fs::File),
        Stdout(std::io::Stdout),
    }

    impl Write for AnyOutput {}
}

fn main() {
    let args = ClassifierArgs::parse();
    let use_stderr = args.output_file.is_none();

    let query_reader = FastaReader::from_path(&args.fasta_file)
        .unwrap_or_die("Cannot open FASTA file!")
        .map(|res| res.map(FastaSeq::filter_to_dna));

    let module = get_sswsort_module(&args.module).unwrap_or_die("Failed to load module data!");
    let canonical_module = module.name.as_str();

    if let Some(n) = args.submit_grid_job
        && let Some(output) = args.output_file
    {
        submit_job_sync(n, canonical_module, &args.fasta_file, output).unwrap_or_die("Qsub job submission failed!");
        return;
    }

    time_stamp(
        &format!(
            "Doing sort on '{query}' with module '{canonical_module}'.",
            query = args.fasta_file.file_name().unwrap_or(args.fasta_file.as_os_str()).display(),
        ),
        use_stderr,
    );

    let mut grid = None;
    let mut w = std::io::BufWriter::new(if let Some(mut path) = args.output_file {
        if args.is_grid_task {
            let g = Grid::from_env();
            let id = g.task_id;
            grid = Some(g);

            path.set_file_name(get_partition_filename(&path, id));
        } else if path.is_dir() {
            path.set_file_name("sswsort_output.tsv");
        }
        AnyOutput::File(std::fs::File::create(path).unwrap_or_die("Cannot create output file!"))
    } else {
        AnyOutput::Stdout(std::io::stdout())
    });

    query_reader
        .process_results(|queries| {
            if let Some(g) = grid {
                // 0-based skip so we start on our parition
                let offset = (g.task_id - g.task_first) / g.task_stepsize;
                // modulus for interleaved partitioning
                let array_size = (g.task_last - g.task_first + 1).div_ceil(g.task_stepsize);

                let iter_results = queries.skip(offset).step_by(array_size).map(|query| {
                    (
                        classify_top_two_with_warning(&module, &query),
                        query.name,
                        query.sequence.len(),
                    )
                });

                write_results(&mut w, iter_results, canonical_module).unwrap_or_die("Could not write to file!");
            } else if args.threads.is_none_or(|n| n.get() > 1) {
                let t = if let Some(n) = args.threads {
                    n.get()
                } else if let Some(v) = std::env::var("IFX_LOCAL_PROCS").ok().and_then(|v| v.parse::<usize>().ok()) {
                    v
                } else {
                    num_cpus::get_physical()
                };
                ThreadPoolBuilder::new().num_threads(t).build_global().unwrap();

                let results: Vec<([ClassificationResult; 2], String, usize)> = queries
                    .par_bridge()
                    .map(|query| {
                        (
                            classify_top_two_with_warning(&module, &query),
                            query.name,
                            query.sequence.len(),
                        )
                    })
                    .collect();

                write_results(&mut w, results.into_iter(), canonical_module).unwrap_or_die("Could not write to file!");
            } else {
                let iter_results = queries.map(|query| {
                    (
                        classify_top_two_with_warning(&module, &query),
                        query.name,
                        query.sequence.len(),
                    )
                });
                write_results(&mut w, iter_results, canonical_module).unwrap_or_die("Could not write to file!");
            }
        })
        .unwrap_or_die("Failed to read input file");

    w.flush().expect("Flushing failed, there may be missing data");

    time_stamp("finished", use_stderr);
}

/// Performs classification with [`SSWSortModule::classify_top_two`], printing
/// to `stderr` in case of empty input.
fn classify_top_two_with_warning<'a>(module: &'a SSWSortModule, query: &FastaNT) -> [ClassificationResult<'a>; 2] {
    if query.sequence.is_empty() {
        eprintln!("WARNING: '{name}' was empty.", name = query.name);
    }

    module.classify_top_two(&query.sequence)
}
