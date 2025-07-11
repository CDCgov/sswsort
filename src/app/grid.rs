use std::env;
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Grid Engine environment variables for array job processing.
#[derive(Debug, Clone)]
pub(crate) struct Grid {
    pub task_id:       usize,
    pub task_first:    usize,
    pub task_last:     usize,
    pub task_stepsize: usize,
}

impl Grid {
    /// Extract SGE environment variables from the current environment.
    ///
    /// ## Panics
    ///
    /// Panics if required variables `SGE_TASK_ID`` or `SGE_TASK_LAST` are not set.
    /// Uses default values for optional variables if not present.
    pub fn from_env() -> Self {
        Self {
            task_id:       env::var("SGE_TASK_ID")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .expect("Invalid SGE_TASK_ID!"),
            task_first:    env::var("SGE_TASK_FIRST")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(1),
            task_last:     env::var("SGE_TASK_LAST")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .expect("Invalid SGE_TASK_LAST!"),
            task_stepsize: env::var("SGE_TASK_STEPSIZE")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or(1),
        }
    }
}

pub(crate) fn get_partition_filename(path: &Path, id: usize) -> String {
    if let Some(stem) = path.file_stem() {
        let stem = stem.to_string_lossy();
        if let Some(ext) = path.extension() {
            format!("{stem}_{id:03}.{e}", e = ext.to_string_lossy())
        } else {
            format!("{stem}_{id:03}.tsv")
        }
    } else {
        format!("sswsort_output_{id:03}.tsv")
    }
}

/// Submit a job to Grid Engine and block until completion.
///
/// # Arguments
/// * `task_count` - Number of array tasks to submit
/// * `input_path` - Path to the input file
/// * `output_path` - Path to the output file
///
/// # Panics
/// Panics if the qsub command fails or returns a non-zero exit code.
pub fn submit_job_sync(
    task_count: usize, module: &str, input_path: PathBuf, mut output_path: PathBuf,
) -> Result<(), std::io::Error> {
    let current_exe = env::current_exe().expect("Failed to get current executable path");

    let mut log_path = output_path.clone();
    let log_file = log_path
        .file_stem()
        .map(|s| format!("{}_log.txt", s.to_string_lossy()))
        .unwrap_or_else(|| "sswsort_log.txt".to_string());
    log_path.set_file_name(log_file);
    let (exec, input, output, log) = (
        current_exe.to_str().expect("Invalid executable path"),
        input_path.to_str().expect("Invalid input path"),
        output_path.to_str().expect("Invalid output path"),
        log_path.to_str().expect("Log invallid and failed"),
    );

    let mut cmd = Command::new("qsub");
    cmd.args([
        "-t",
        &format!("1-{task_count}"),
        "-sync",
        "yes",
        "-cwd",
        "-j",
        "yes",
        "-o",
        log,
        "-b",
        "yes",
        exec,
        module,
        input,
        output,
        "--is-grid-task",
    ])
    .stdout(Stdio::piped())
    .stderr(Stdio::piped());

    eprintln!("Executing: {cmd:#?}");

    let output_cmd = cmd.output()?;
    if !output_cmd.status.success() {
        let stderr = String::from_utf8_lossy(&output_cmd.stderr);
        panic!("Grid Engine submission failed: {stderr}",);
    }

    let mut partition_file = output_path.clone();
    let mut partition_data = Vec::new();

    if output_path.file_name().is_none() {
        output_path.set_file_name("sswsort_output.tsv");
    }
    let mut output_file = File::create(&output_path)?;

    for id in 1..=task_count {
        partition_file.set_file_name(get_partition_filename(&output_path, id));

        let mut file = File::open(&partition_file)?;
        file.read_to_end(&mut partition_data)?;
        output_file.write_all(&partition_data)?;
        output_file.flush()?;

        std::fs::remove_file(&partition_file)?;
        partition_data.clear();
    }

    Ok(())
}
