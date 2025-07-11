use sswsort::*;

// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

// Takes an iterator of (ClassificationResult, String, usize) and print lints to
// the provided writer.
pub(crate) fn write_results<I, W>(writer: &mut W, results: I, module_name: &str) -> std::io::Result<()>
where
    I: Iterator<Item = (ClassificationResult, String, usize)>,
    W: std::io::Write, {
    for (classification, query_name, seq_len) in results {
        writeln!(
            writer,
            "{PROGRAM_VERSION}\t{module_name}\t{query_name}\t{classification}\t{seq_len}",
        )?;
    }
    Ok(())
}
