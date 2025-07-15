# Migration Guide

This guide describes differences in behavior between SSWSORT 1 and 2.

## Installation

Previously it was sufficient to clone the repository to launch the process. Now there are two options:

1. Download the package from the CDCGOV [releases] section and install it.
2. Clone the package as before and run the build script. Please note that you might have a recent nightly version of Rust.

## Usage

## Output

- Unusually long data is now denoted: `*Unusually Long: {taxon}` rather than `*{taxon}`
- Chimeric data is now denoted `*Chimeric: {taxon1}+{taxon2}+...` rather than `*{taxon1}+{taxon2}`
- Chimeric data now uses `\N` for the *best score* (Null in Hive / Impala)
- When in STDOUT mode, logging is redirected to STDERR

## Configuration

[releases]: https://github.com/CDCgov/sswsort