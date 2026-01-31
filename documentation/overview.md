# DELT-Hit codebase overview

## Project context
DELT-Hit is an end-to-end, open-source toolkit for DNA-encoded library (DEL) analysis. The protocol outlines a five-module pipeline that converts raw sequencing reads into analysis-ready chemical information: adaptive demultiplexing, reaction-based structure reconstruction, molecular property/descriptor generation, statistical hit ranking, and integrated QC/visualization. The implementation in this repository mirrors those modules through CLI commands and supporting packages, so each step is callable as part of a reproducible workflow.

## How the pipeline maps to the repository

### Core packages (`src/delt_hit`)
- **`cli/`**: CLI entrypoints and subcommands (`delt-hit <group> <method>`). These bind to the modules below and are the primary way users interact with the pipeline.
- **`demultiplex/`**: FASTQ demultiplexing, Cutadapt input generation, and post-processing of counts. The main data products are selection-level count tables in TSV format.
- **`library/`**: Reaction graph construction, building-block enumeration, and library-level property/descriptor generation. Output is a `library.parquet` dataset with computed properties and fingerprints.
- **`quality_control/`**: Plotting and summary reports tied to demultiplexing and hit inspection.
- **`utils.py`**: Shared YAML I/O helpers and small utilities used throughout the CLI.

### Supporting assets
- **`templates/`**: Excel templates used to initialize configuration via `delt-hit init`.
- **`protocols.pdf`**: The reference protocol describing the intended scientific workflow and rationale for each module.
- **`tests/`**: Unit tests that validate core building blocks.

## High-level data flow
1. **Configuration**: Start from an Excel library template, then run `delt-hit init` to produce a `config.yaml` that captures experiment metadata, library structure, reaction catalog, and barcode/primer definitions.
2. **Demultiplexing**: Use `delt-hit demultiplex prepare/run` to generate and execute Cutadapt workflows, then `delt-hit demultiplex process` to convert adapter-tagged reads into per-selection count tables.
3. **Library enumeration**: Use `delt-hit library enumerate` to build the reaction graph and generate SMILES for the complete library. This creates `library.parquet` and reaction graph visualizations.
4. **Property/descriptor computation**: `delt-hit library properties` and `delt-hit library represent` produce cheminformatics features (properties and fingerprints) for downstream modeling.
5. **Analysis & QC**: `delt-hit analyse enrichment` runs statistical enrichment (edgeR or count-based) and exports result tables, while `delt-hit demultiplex report/qc` and `delt-hit dashboard` provide textual/visual inspection of the results.

For CLI specifics (arguments, outputs, and example commands), see [documentation/cli.md](cli.md).
