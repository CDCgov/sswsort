use jiff::Zoned;
use sswsort::*;

// If we later want to match the shell script, we can use:
// <https://docs.rs/git-version/latest/git_version/macro.git_describe.html>
const PROGRAM_VERSION: &str = concat!(env!("CARGO_PKG_NAME"), " v", env!("CARGO_PKG_VERSION"));

/// Output the top two classification results for the given query.
///
/// This is a helper function for [`write_results`]. Argument `t` is a `&str`,
/// intended to be either empty or contain a `\t`, depending if the
/// `ClassificationResult` is primary (outputs a sequence length between the
/// score and strand) or secondary (does not, in which case `seq_len` should
/// also be an empty string).
fn write_result<W>(
    writer: &mut W, classification_result: &ClassificationResult<'_>, t: &str, seq_len: &str,
) -> std::io::Result<()>
where
    W: std::io::Write, {
    match classification_result {
        ClassificationResult::Classification {
            taxon,
            best_score,
            strand,
        } => write!(writer, "{taxon}\t{best_score}{t}{seq_len}\t{strand}"),
        ClassificationResult::Unrecognizable { best_score } => {
            if let Some(best_score) = best_score {
                write!(writer, "UNRECOGNIZABLE\t{best_score}{t}{seq_len}\t*")
            } else {
                write!(writer, "UNRECOGNIZABLE\t\\N{t}{seq_len}\t*")
            }
        }
        ClassificationResult::Chimeric { taxa } => {
            let taxa_list = taxa.join("+");
            write!(writer, "*Chimeric: {taxa_list}\t\\N{t}{seq_len}\t*")
        }
        ClassificationResult::UnusuallyLong {
            taxon,
            best_score,
            strand,
        } => write!(writer, "*Unusually Long: {taxon}\t{best_score}{t}{seq_len}\t{strand}"),
        ClassificationResult::Unresolvable { taxa, best_score } => {
            let taxa_list = taxa.join(",");
            write!(writer, "*Unresolvable: {taxa_list}\t{best_score}{t}{seq_len}\t*")
        }
    }
}

/// Takes an iterator of classification results `[(usize, &str, Strand); 2]` and
/// prints to the provided `writer`, in tab-separated format with fields:
///
/// `program_version, module_name, query_name`
///
/// Followed by the classification results with the following fields:
///
/// `classification1, score1, sequence_length, strand1, classification2, score2,
/// strand2`
pub(crate) fn write_results<'a, I, W>(writer: &mut W, results: I, module_name: &str) -> std::io::Result<()>
where
    I: Iterator<Item = ([ClassificationResult<'a>; 2], String, usize)>,
    W: std::io::Write, {
    for (classification, query_name, seq_len) in results {
        write!(writer, "{PROGRAM_VERSION}\t{module_name}\t{query_name}\t")?;

        let mut buff = core::fmt::NumBuffer::new();
        let seq_len = seq_len.format_into(&mut buff);

        write_result(writer, &classification[0], "\t", seq_len)?;
        write!(writer, "\t")?;
        write_result(writer, &classification[1], "", "")?;
        writeln!(writer)?;
    }
    Ok(())
}

/// Creates a time stamp to be used on initialization and completion of
/// classification
pub(crate) fn time_stamp(message: &str, use_stderr: bool) {
    // Pad the timestamp based on the shell level using environment variable
    // SHLVL
    let shlvl = std::env::var("SHLVL").ok().and_then(|v| v.parse::<usize>().ok()).unwrap_or(0);
    let pad = "  ".repeat(shlvl.saturating_sub(1));

    if use_stderr {
        eprintln!(
            "[{now}] {pad}SSWSORT :: {message}",
            now = Zoned::now().strftime("%Y-%m-%d %k:%M:%S")
        );
    } else {
        println!(
            "[{now}] {pad}SSWSORT :: {message}",
            now = Zoned::now().strftime("%Y-%m-%d %k:%M:%S")
        );
    }
}
