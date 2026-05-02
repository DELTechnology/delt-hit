# Changelog

## 2026-04-30 to 2026-05-02

### Added
- QC codon-hit plotting now exports the manuscript-aligned image set needed for the revision package, with regression coverage to lock in the expected output behavior.
- Enrichment plotting scripts now export PNG outputs directly for the Favalli and pure-DEL supporting-material analyses so revision figures can be reused without extra conversion steps.

### Improved
- Published-dataset supporting-material workflows were aligned and documented more clearly, including streamlined run scripts and revised README guidance for reproducing the revision analyses.
- Favalli CA9 analysis selections were updated to match the current manuscript re-analysis workflow and supporting-material configuration.
- Enrichment rerun documentation and output locations were cleaned up so the generated figures and tables map cleanly onto the revision deliverables.
- Manuscript figure/output handling was aligned with the latest revision pass, including consistent QC image export behavior and updated LaTeX references in `paper/main.tex`.

## 2026-04-24 to 2026-04-29

### Added
- New `visualize enumerate` CLI workflow for reviewer-friendly chemistry outputs, including reaction graphs, reaction schemes, per-building-block structure panels, and configured compound structure panels.
- New `visualize library` CLI workflow for exporting named-library members as individual structure PNGs with filenames derived from their DEL code combinations (for example `B0=1-B1=0.png`).
- Explicit dual-display DEL support across library parsing, enumeration, and visualization workflows, including strand-aware building-block validation and dual-display combination handling.
- Filtered library enumeration from observed `counts.txt` inputs with `--counts_path`, `--top_n`, and `--library_name`, plus named-library support for downstream property generation.
- Synthetic benchmarking utilities for demultiplexing and enrichment, including FASTQ generation, DELT/DELi input converters, runtime plotting, worker/chunk-size sweeps, and large-library benchmark matrices.
- Additional experiment templates and supporting-material runs for Favalli and related comparison workflows.
- New helper utilities, benchmark tests, synthetic-enrichment tests, FASTQ generator tests, filtered-enumeration tests, and demultiplex parsing/preprocess coverage.

### Improved
- Visualization output layout and export locations were reorganized to produce cleaner library-level output folders and more consistent graph/property save paths.
- `visualize enumerate` now writes building-block outputs to `library/visualization/building_blocks/<BUILDING_BLOCK_ID>/<index>.png` instead of one combined grid per building-block family.
- `visualize library` now writes per-entry outputs to `library/visualization/libraries/<library_name>/`, with filenames derived from the `code_*` columns so outputs map directly to DEL building-block combinations.
- Visualization commands now consistently render one structure per PNG for building blocks, named-library members, and configured compounds, avoiding oversized combined structure sheets.
- Enumeration was migrated toward the new visualization-oriented workflow while preserving the core `library enumerate` command.
- Building-block handling now supports reverse-complement/complement codon configurations needed for B-sheet style inputs.
- Demultiplex and library inputs were made more flexible, including support for simpler input files, direct `counts.txt` enumeration inputs, and cleaner tool output.
- Chemical property generation and CLI/docs now better describe named libraries, property output structure, and enumeration behavior.
- Benchmark orchestration was streamlined with dedicated benchmark folders, reusable run scripts, configurable loops, estimated iteration controls, and result-copyback helpers.
- The project environment was modernized by switching to Pixi-managed lockfiles and cleaning out older dependency setup files.

### Documented
- Expanded CLI, benchmarking, z-score, chemical-property, DELi-option, and enumeration-improvement documentation.
- Added reviewer-oriented examples and workflow notes for `visualize enumerate`, `visualize library`, named-library outputs, and synthetic benchmarking.
- Updated supporting-material example run scripts to reflect the new visualization commands and current output paths for filtered named-library workflows.
