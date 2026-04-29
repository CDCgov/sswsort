# SSWSort 2 - simple virus gene segment / genome classification via striped Smith-Waterman

SSWSort 2 is provided as both a Rust library, which can be used as a dependency
in other projects, and as an executable binary, which can be run via
command-line and replaces SSWSORT. The following documentation discusses usage
and outputs for the binary, `sswsort`. For information on the library, read
the [Rust docs].

Classifies (or sorts) sequences (influenza, SARS-CoV-2, RSV) using presets or
DBs and a query sequence (FASTA). The classification is a best match criterion,
and ignores same-strand repeats. SSWSORT will reject sequences that are
unexpectedly long, chimeric (optional), and beyond the scoring thresholds.

## Usage

```bash
Uses striped Smith-Waterman to classify sequences into a simple compound type

Usage: sswsort [OPTIONS] <MODULE> <FASTA_FILE> [OUTPUT_FILE]

Arguments:
  <MODULE>       Name of the classification module
  <FASTA_FILE>   Name of the nucleotide sequences to classify in FASTA format
  [OUTPUT_FILE]  Name of the tab-separated-value file for classifier results. If none are provided, STDOUT is used. If a directory is specified, a default filename of `sswsort_output.tsv`

Options:
  -T, --threads <THREADS>       Number of threads to use. Defaults to number of physical cores otherwise
  -G, --is-grid-task            Execute as a partitioned task in a grid job, for use with: --submit-grid-job
  -S, --submit-grid-job <SIZE>  Submits and blocks on a grid job of the specified array size
  -h, --help                    Print help (see more with '--help')
  -V, --version                 Print version
```

### Self-tests

```bash
./sswsort flu sswsort_res/flu.fasta
./sswsort cov-beta sswsort_res/cov-beta.fasta
./sswsort rsv sswsort_res/rsv.fasta
```

## Installation

### Compile from source

Clone this repo and, optionally, check out the latest release tag. You must
install [Rust](https://rust-lang.org) *nightly* as a pre-requisite. You can then
simply run the installer script:

```bash
./install.sh
```

### Via Archive

1. Download the latest archive via our [releases page](https://github.com/CDCgov/sswsort/releases)
   for your target platform.
2. Unzip the archive containing `sswsort`.
3. Move the package to your desired location and add the folder to your `PATH`
   - Note: `sswsort_res` and `sswsort` must be in the same folder.

### Via Container

Simply run:

```bash
## From Github Container Repo
docker run --rm -itv $(pwd):/data ghcr.io/cdcgov/sswsort:latest sswsort # more args
```

## Outputs

SSWSORT will provide output in a tab-separated format, with columns representing:

Program version, Reference module, Query name, Primary classified taxon, Primary
classification score, Sequence length, Primary strand, Secondary classified
taxon, Secondary score, Secondary strand

To receive a classification, the query must align with one of the reference
sequences with a score greater than or equal to the `score_minimum`, or a
normalized score greater than or equal to the `norm_score_minimum`, equal to the
score divided by the length of the query. These parameters have default values
that can be changed in the `config.toml`.

Classifications fall into the following categories:

| Classification                                    | Description                                                                                                                                                                                                                                         |
| ------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A single taxon                                    | The taxon with the highest score for the provided query                                                                                                                                                                                             |
| `*Unusually Long` with a single taxon             | The query length is over twice the length of the reference it is matched to                                                                                                                                                                         |
| `*Chimeric` with a `+`-separated list of taxa     | The query matches to multiple different reference taxa with scores above 800 each. Chimera detection can be disabled in `config.toml`                                                                                                               |
| `*Unresolvable` with a `,`-separated list of taxa | The query matches to multiple reference taxa with exactly matching scores                                                                                                                                                                           |
| `UNRECOGNIZABLE`                                  | The query did not match to any reference with a high enough score or normalized score to pass the threshold. Queries with `*Chimeric` or `*Unresolvable` primary classifications are also given `UNRECOGNIZABLE` for their secondary classification |

## Alignment and Scoring

SSWSORT utilizes a Striped Smith-Waterman alignment algorithm to align each
query sequence to a list of references before classifying the query with the
reference taxa that had the highest alignment score. The alignment uses these
default following weights for scoring:

| Parameter  | Weight |
| ---------- | ------ |
| Match      | 2      |
| Mismatch   | -5     |
| Gap open   | -10    |
| Gap extend | -1     |

Defaults can be overridden at a module level within the `config.toml`. 

Example:

```toml
[[classification_module]]
name = "flu"
version = "2.0"
alternative_names = []
norm_score_minimum = 1.0
score_minimum = 100
length_minimum = 25
reference_sequences = "flu.fasta"
detect_chimera = true
weights = { gap_open = -11, gap_extend = -2, mismatch = -3, match_weight = 1 }
```
will override the values for all weights. Note if some but not all weights are
being overridden, the other values still need to be included in the `toml`.

Ambiguous nucleotides (`N`) and unrecognized characters are treated as 0-penalty mismatches.
So, a query with one or more `N` bases will align with an identical score as a
query with those bases missing.

Leading and trailing `N`'s in a query sequence are removed prior to alignment.
This will not affect the alignment score, but will affect the sequence length in
the output, as well as the normalized score, which uses the edited length.

## Notices

### Contact Info

For direct correspondence on the project, feel free to contact: [Samuel S. Shepard](mailto:sshepard@cdc.gov), Centers for Disease Control and Prevention or reach out to other [contributors](CONTRIBUTORS.md).

### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).  All contributions to this repository will be released under the CC0 dedication.  By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

### License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later. This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version. This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details. You should have received a copy of the Apache Software License along with this program. If not, see: <http://www.apache.org/licenses/LICENSE-2.0.html>. The source code forked from other open source projects will inherit its license.

### Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the [Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md). For more information about CDC's privacy policy, please visit <http://www.cdc.gov/other/privacy.html>.

### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git).) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
