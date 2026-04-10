# Migration Guide

This guide describes differences in behavior between SSWSort 1 and 2. Notably, SSWSort is now released as a Rust library, which can be used as a dependency in other projects, as well as a command-line tool, meant to replace the original SSWSort. This guide discusses usage for the command-line executable binary provided.

## Installation

Previously it was sufficient to clone the repository to launch the process. Now there are two options:

1. Download the SSWSort-CLI binary package from the CDCGOV [releases] section and install it.

## Usage

New flags `--threads` (`-T`) for setting threads, `--submit-grid-job` (`-S`)
for blocking a grid engine job of specified array size, and `--is-grid-task`
(`-G`) for detecting array size and submitting have been added.

### Example input

To run SSWSort2, arguments are provided in the format:

```bash
sswsort <MODULE> <FASTA_FILE> [OUTPUT_FILE]
```

For a demo, and to the classifier reference set against itself, run:

```bash
sswsort flu sswsort_res/flu.fasta out.tsv
```

which will classify all references and output primary and secondary classification information to `out.tsv`. If no output path is provided, the results will print to `STDOUT`.

## Output

- Three new columns have been added to the output, representing a secondary
  classification, secondary score, and secondary strand for the second-best
  classification
- A new classification has been added: `*Unresolvable: {taxon},{taxon},...`.
  This will occur in cases where the query aligns to multiple different
  reference taxa with identical scores
- Unusually long data is now denoted: `*Unusually Long: {taxon}` rather than
  `*{taxon}`
- Chimeric data is now denoted `*Chimeric: {taxon1}+{taxon2}+...` rather than
  `*{taxon1}+{taxon2}`
- Chimeric data now uses `\N` for the *best score* (Null in Hive / Impala)
- In cases where the primary classification is unresolvable or chimeric, the
  secondary classification provided will be `UNRECOGNIZABLE` with a `\N` score
- When in STDOUT mode, logging is redirected to STDERR

## Configuration

- For configuration, a module still needs to provided. The main module options are `flu`, `cov`, `spike` (for SARS-CoV-2 spike protein), and `rsv`; each of these has aliases provided in `sswsort_res/config.toml`.

[releases]: https://github.com/CDCgov/sswsort/releases
