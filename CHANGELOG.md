# Changelog

## 2026-04-24 to 2026-04-29

### Added
- New `visualize enumerate` CLI workflow for reviewer-friendly chemistry outputs, including reaction graphs, reaction schemes, building-block grids, and compound structure panels.
- Filtered library enumeration from observed `counts.txt` inputs with `--counts_path`, `--top_n`, and `--library_name`, plus named-library support for downstream property generation.
- Synthetic benchmarking utilities for demultiplexing and enrichment, including FASTQ generation, DELT/DELi input converters, runtime plotting, worker/chunk-size sweeps, and large-library benchmark matrices.
- Additional experiment templates and supporting-material runs for Favalli and related comparison workflows.
- New helper utilities, benchmark tests, synthetic-enrichment tests, FASTQ generator tests, filtered-enumeration tests, and demultiplex parsing/preprocess coverage.

### Improved
- Visualization output layout and export locations were reorganized to produce cleaner library-level output folders and more consistent graph/property save paths.
- Enumeration was migrated toward the new visualization-oriented workflow while preserving the core `library enumerate` command.
- Building-block handling now supports reverse-complement/complement codon configurations needed for B-sheet style inputs.
- Demultiplex and library inputs were made more flexible, including support for simpler input files, direct `counts.txt` enumeration inputs, and cleaner tool output.
- Chemical property generation and CLI/docs now better describe named libraries, property output structure, and enumeration behavior.
- Benchmark orchestration was streamlined with dedicated benchmark folders, reusable run scripts, configurable loops, estimated iteration controls, and result-copyback helpers.
- The project environment was modernized by switching to `pixi`/`uv`-style lockfiles and cleaning out older dependency setup files.

### Documented
- Expanded CLI, benchmarking, z-score, chemical-property, DELi-option, and enumeration-improvement documentation.
- Added reviewer-oriented examples and workflow notes for library visualization and synthetic benchmarking.
