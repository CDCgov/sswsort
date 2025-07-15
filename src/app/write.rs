use jiff::Timestamp;
use sswsort::*;
// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

// Takes an iterator of (ClassificationResult, String, usize) and print lints to
// the provided writer.
pub(crate) fn write_results<'a, I, W>(writer: &mut W, results: I, module_name: &str) -> std::io::Result<()>
where
    I: Iterator<Item = (ClassificationResult<'a>, String, usize)>,
    W: std::io::Write, {
    for (classification, query_name, seq_len) in results {
        write!(writer, "{PROGRAM_VERSION}\t{module_name}\t{query_name}\t",)?;
        match classification {
            ClassificationResult::Classification {
                taxon,
                best_score,
                strand,
            } => {
                writeln!(writer, "{taxon}\t{best_score}\t{seq_len}\t{strand}")?;
            }
            ClassificationResult::Unrecognizable { best_score, strand } => {
                writeln!(writer, "UNRECOGNIZABLE\t{best_score}\t{seq_len}\t{strand}")?;
            }
            ClassificationResult::Chimeric { taxa, strand } => {
                writeln!(
                    writer,
                    "*Chimeric: {taxa_list}\t\\N\t{seq_len}\t{strand}",
                    taxa_list = taxa.join("+")
                )?;
            }
            ClassificationResult::UnusuallyLong {
                taxon,
                best_score,
                strand,
            } => {
                writeln!(writer, "*Unusually Long: {taxon}\t{best_score}\t{seq_len}\t{strand}")?;
            }
        };
    }
    Ok(())
}

pub(crate) fn time_stamp(message: &str, use_stderr: bool) {
    let shlvl = std::env::var("SHLVL").ok().and_then(|v| v.parse::<usize>().ok()).unwrap_or(0);
    let pad = "  ".repeat(shlvl.saturating_sub(1));

    if use_stderr {
        eprintln!(
            "[{now}] {pad}SSWSORT :: {message}",
            now = Timestamp::now().strftime("%Y-%m-%d %k:%M:%S")
        );
    } else {
        println!(
            "[{now}] {pad}SSWSORT :: {message}",
            now = Timestamp::now().strftime("%Y-%m-%d %k:%M:%S")
        );
    }
}
