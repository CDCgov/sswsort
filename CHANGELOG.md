# SSWSORT Changelog

All notable changes to this project will be documented in this file. The format
is roughly based on [Keep a Changelog], and this project tries to adheres to
[Semantic Versioning].

## [2.1.0] - TBD

### Fixed

- Fixed upstream bug where piped inputs were being read as empty and throwing
  error
- Fixed SSW scoring weights in README
- Made classification methods infallible

### Changed

- Removed redundant `Strand::Unknown` enum variant in library. This does not
  change any outputs in the command-line tool
- Generalized public API in `lib` to take `AsRef<[u8]>` for classification
  functions instead of requiring a `FastaNT`
- Dockerfile now uses hardened images (requires login to dhi.io)

## [2.0.0] - 2026-04-13

### Changed

- Rewrites SSWSort into Rust with a provided binary CLI tool and library. See
  `Migration.md` for details about updating existing pipelines to use SSWSort2,
  and for differences with outputs.

## [1.6.3] 2024-12

### Fixed

- Grid engine execution was broken from v1.4+, now it is fixed.

## [1.6.2] 2024-08

### Changed

- RSV compound types renamed to match the conventions of other tooling, e.g., RSVA becomes RSV_AD.

## [1.6.1] 2024-05

### Added

- adds configurable working directories with `IFX_WORK_DIR`.

## [1.6.0] 2024-04

### Added

- adds an RSV module. Many thanks to C. Paden and M. Mandal!

## [1.5.0] 2022-09

- **Fix:** Flu reference data reduced and refined for public access.
- **Change:** chimera annotation list now sorts by annotation name.

## [1.4.0] 2022-07

- Removed git dependency, re-licensed, cleaned and formatted, better error output.

## [1.3.0] 2020-11

- Added support for SC2 (thanks to K. Lacek and G. Stott). Upgrade GNU Parallel and SSW.

## [1.2.0] 2019-01

- Re-factored to make more self-contained

## [1.1.0] 2018-08

- Cleaned up output and added fields such as program version.

## [1.0.0] 2018-08

- Initial release after alternatives testing. Smith-Waterman gave the most sensitive results for this purose.

<!-- Versions -->
[2.1.0]: https://github.com/CDCgov/sswsort/compare/v2.0.0...v2.1.0
[2.0.0]: https://github.com/CDCgov/sswsort/compare/v1.6.3...v2.0.0

<!-- Links -->
[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html
