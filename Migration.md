# Migration Guide

This guide describes differences in behavior between SSWSORT 1 and 2.

## Installation

Previously it was sufficient to clone the repository to launch the process. Now there are two options:

1. Download the package from the CDCGOV [releases] section and install it.
2. Clone the package as before and run the build script. Please note that you might have a recent nightly version of Rust.

## Usage

New flags `--threads` (`-T`) for multi-threaded mode, `--submit-grid-job` (`-S`)
for blocking a grid engine job of specified array size, and `--is-grid-task`
(`-G`) for detecting array size and submitting have been added.

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

[releases]: https://github.com/CDCgov/sswsort