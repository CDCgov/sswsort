# SSWSORT - simple virus gene segment / genome classification via Smith-Waterman

Sorts (classifies) sequences (influenza, SARS-CoV-2, RSV) using presets or DBs and a query sequence (FASTA). The SORT is a best match criterion, and ignores same strand repeats. SSWSORT is vulnerable to anomalous data.

## Usage

Usage:

```bash
SSWSORT <ref> <query> <output_file>
```

Self-tests:

```bash
./SSWSORT flu-ABCD90P presets/flu-ABCD90P.fasta test.txt
./SSWSORT cov-beta presets/cov-beta.fasta test2.txt
./SSWSORT rsv presets/rsv.fasta test3.txt
```

## Component Citations

This package uses a modified version of SSW (`bin/component_licenses_and_source`):

> Zhao M, Lee WP, Garrison EP, Marth GT. SSW library: an SIMD Smith-Waterman C/C++ library for use in genomic applications. PLoS One. 2013;8(12):e82138. Published 2013 Dec 4. doi:10.1371/journal.pone.0082138

This package uses GNU Parallel :

> Tange O. GNU Parallel. Version 20200422. doi:10.5281/zenodo.1146014
