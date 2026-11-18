# Changelog

All notable changes to DELT-Hit are documented here. Versions follow the Git tags.

## Unreleased

### Removed
- Internal notes, manuscript files and unpublished data from the repository and its history.
- Unused `python-docx` dependency.

### Changed
- Generated outputs of the supporting-material experiments are no longer tracked.

## 1.1.1 - 2026-09-10

### Changed
- Conda is the default installation method; Pixi is documented as an alternative, including Graphviz plugin registration.
- The example sequencing download shows its progress.
- Updated the Pixi lockfile.

## 1.1.0 - 2026-06-05

### Changed
- Supporting-material workflows select the comparisons that match the input FASTQ files.
- Replaced the old Excel templates with `templates/single-display-two-cycle.xlsx`.

### Removed
- The enrichment benchmark.

## 1.0.0 - 2026-06-01

First public release.

### Added
- End-to-end workflow driven by one Excel workbook: demultiplexing, library enumeration, chemical properties, enrichment analysis (counts, edgeR, z-score) and quality control.
- `visualize enumerate` for reaction graphs, reaction schemes and building-block structures, and `visualize library` for one structure image per library member.
- Dual-display DEL support in library parsing, enumeration and visualization, including strand-aware building-block validation and reverse-complement codons.
- Filtered enumeration from observed counts (`--counts_path`, `--top_n`, `--library_name`) and named libraries for downstream property generation.
- Synthetic benchmarking tools for demultiplexing, including FASTQ generation and converters for DELi inputs.
- Example experiments and published-dataset workflows under `supporting_material/`.

### Changed
- Visualization outputs write one structure per PNG into per-library folders.
- Parallel, chunked count parsing for large sequencing files.
- Environment management moved to Pixi lockfiles.
