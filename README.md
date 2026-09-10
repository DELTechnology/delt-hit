# DELT-Hit

DELT-Hit is an open-source, end-to-end computational framework for DNA-encoded chemical library analysis. It connects sequence demultiplexing, chemical structure reconstruction, molecular property calculation, enrichment analysis, and quality control.

## 1. System requirements

The protocol specifies Python 3.12+, R 4.1+, and Cutadapt 4.9+, with Conda for environment management. R is required for all three enrichment methods; edgeR is additionally required for the edgeR method.

### Tested versions

The example workflow and paper analyses were run on **macOS 26.6.2** with:

| Software | Version |
|---|---|
| Python | 3.12.13 |
| Cutadapt | 5.2 |
| R | 4.4.0 |
| edgeR | 4.4.2 |

The complete Python and R dependency list is in [pyproject.toml](pyproject.toml). Python dependencies are installed automatically with DELT-Hit. For Pixi, R and its packages (tidyverse, GGally, edgeR, limma) are installed automatically as well; for Conda, Graphviz/pygraphviz and the R packages are installed separately as described below. The R workflow uses tidyverse and GGally, plus edgeR and limma for the edgeR method.

Minimum hardware: 16 GB RAM, 8 CPU cores, and 50 GB available storage. Recommended: 32 GB RAM, 16 CPU cores, and 100 GB available storage. Large datasets may require additional memory and disk space.

## 2. Installation guide

### Option A: Conda (default)

Install [Miniconda](https://docs.anaconda.com/miniconda/) for your operating system and initialize your shell during installation. Then create and activate an isolated environment, as in Box 1 of the protocol:

```bash
conda create -n delt-hit python=3.12 -y
conda activate delt-hit
conda install -c conda-forge pygraphviz -y
pip install git+https://github.com/DELTechnology/delt-hit.git
```

Activate the environment with `conda activate delt-hit` in each new terminal session.

### Option B: Pixi

Install [Pixi](https://pixi.sh/latest/), then clone the repository and install its configured environment:

```bash
git clone https://github.com/DELTechnology/delt-hit.git
cd delt-hit
pixi install
pixi run dot -c  # registers Graphviz plugins; Pixi skips the post-install step Conda runs automatically
pixi run delt-hit --help
pixi shell
```

`pixi shell` activates the environment so the commands and supporting-material scripts below can be run directly. Use it instead of `conda activate delt-hit` when following this README or the experiment instructions. Alternatively, prefix commands with `pixi run`, including `pixi run bash run.sh` for a complete workflow. The configured Pixi platforms are Linux x86-64 and macOS Apple Silicon.

### R dependencies and verification

For enrichment analysis, R 4.1+ with tidyverse, GGally, edgeR, and limma is required.

The Pixi environment (Option B) installs R and all required R packages automatically; no separate step is needed.

For Conda (Option A), install R 4.1+ if it is not already available (for example, `conda install -c conda-forge r-base -y`), then open the corresponding R installation or RStudio and run:

```r
install.packages(c("tidyverse", "GGally"))
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("edgeR", "limma"))
```

Verify the installation:

```bash
delt-hit --help
Rscript --vanilla -e 'library(tidyverse); library(GGally); library(edgeR); library(limma)'
```

Run these checks in the environment activated using your chosen option (`pixi shell` for Pixi, `conda activate delt-hit` for Conda). The first command prints the CLI help; the R check should complete without missing-package errors.

Allow approximately **15 minutes** for installation on a desktop computer, as estimated in the protocol. Download speed and compilation of dependencies can increase this time. This estimate excludes the demo dataset download.

## 3. Demo

### Run the example data

The [single-display example](supporting_material/experiments/example-single-display/) contains the Excel library definition, analysis YAML, and scripts for the workflow in Box 2 of the protocol. The compressed FASTQ download is approximately **6.6 GB**.

```bash
conda activate delt-hit
git clone https://github.com/DELTechnology/delt-hit.git
cd delt-hit/supporting_material/experiments/example-single-display
bash download.sh
bash run.sh
```

If you already cloned the repository, enter the existing example directory instead. Run `download.sh` and `run.sh` from that directory. The [run script](supporting_material/experiments/example-single-display/run.sh) lists each command, including demultiplexing, selection-focused enumeration, enrichment analysis, and optional full-library processing.

Allow **one to several hours** for the demo on a desktop with 32 GB RAM and 16 CPU cores, consistent with the protocol's overall timing guidance. This is an indicative estimate, not a separately measured demo benchmark, and excludes downloading the sequencing data. Each stage in `run.sh` reports its elapsed time. Full-library enumeration and representation generation can be omitted when only selection analysis is required.

### Expected output

The workflow writes the following folders under the example directory:

```text
campaign/
├── config.yaml
├── demultiplex/
│   ├── cutadapt_input_files/    # Barcode files and demultiplex.sh
│   └── cutadapt_output_files/   # Cutadapt logs, JSON reports and mapped reads
├── qc/                         # report.txt and barcode-recovery plots
├── selections/                 # <selection>/counts.txt and flat count exports
├── library/
│   ├── AG24_4_top_hits.parquet
│   ├── library.parquet         # Full-library enumeration
│   ├── visualization/          # Reaction graphs and molecular structures
│   └── properties/             # Descriptor tables and property plots
├── analysis/
│   ├── counts/condition_vs_control/
│   ├── edgeR/condition_vs_control/
│   └── z_score/                # Separate results for AG24_13, AG24_14, AG24_15
└── representations/            # Molecular representation files
```

Selection count tables contain the recovered building-block combinations and their read counts. After the generated R scripts execute, analysis folders contain `stats.csv` and ranked `hits.csv` files, with normalized tables and replicate-correlation plots where applicable. The protocol's Anticipated Results section describes these outputs and their interpretation.

The current `bert` CLI option in `run.sh` writes Morgan fingerprints to `bert.npz`; it does not yet generate the transformer embeddings described in the protocol. Use `morgan.npz` for Morgan fingerprints until that implementation discrepancy is resolved.

The [dual-display example](supporting_material/experiments/example-dual-display/) demonstrates library configuration and enumeration. Its output uses `smiles_a` and `smiles_b`; property and representation commands currently require the single-display `smiles` format.

## 4. Instructions for use

### Prepare your own library and experiment

Copy [templates/single-display-two-cycle.xlsx](templates/single-display-two-cycle.xlsx) and fill in the sheets for your library and selections:

- `experiment`: experiment name, FASTQ path, output directory, and CPU count.
- `selection`: selection identifiers, multiplexing barcodes, and experimental metadata.
- `structure` and `constant`: DNA-region order, constant sequences, and matching error tolerances.
- `B0`, `B1`, and any additional building-block sheets: DNA codons and, for enumeration, building-block SMILES and reaction assignments.
- `reactions`, `compounds`, and `reaction_graph`: reaction SMIRKS, scaffolds, and any additional reaction steps needed for enumeration.

See the [dual-display example](supporting_material/experiments/example-dual-display/) for strand-specific configuration. Validate reaction definitions on a representative subset before enumerating a complete library.

Generate the project configuration:

```bash
delt-hit init --excel_path=/path/to/library.xlsx
```

This creates `<save_dir>/<name>/config.yaml` using the workbook's experiment fields. Review the generated file, then use its path in subsequent commands. The placeholders below must be replaced with your own paths and selection name.

### Demultiplex and inspect sequencing reads

```bash
CONFIG_PATH=/path/to/output/experiment/config.yaml
PROJECT_DIR=/path/to/output/experiment

delt-hit demultiplex prepare --config_path="$CONFIG_PATH"
bash "$PROJECT_DIR/demultiplex/cutadapt_input_files/demultiplex.sh"
delt-hit demultiplex report --config_path="$CONFIG_PATH"
delt-hit demultiplex qc --config_path="$CONFIG_PATH"
delt-hit demultiplex process --config_path="$CONFIG_PATH"
delt-hit demultiplex process --config_path="$CONFIG_PATH" --as_files=True
```

Inspect read retention and barcode recovery before interpreting selection counts. For interactive exploration:

```bash
SELECTION_NAME=your_selection
delt-hit dashboard --config_path="$CONFIG_PATH" \
  --counts_path="$PROJECT_DIR/selections/$SELECTION_NAME/counts.txt"
```

Open the printed local URL, then stop the dashboard with Ctrl+C when finished.

### Reconstruct selected compounds

```bash
delt-hit visualize enumerate --config_path="$CONFIG_PATH"
delt-hit library enumerate \
  --config_path="$CONFIG_PATH" \
  --counts_path="$PROJECT_DIR/selections/$SELECTION_NAME/counts.txt" \
  --top_n=1000 --library_name=top_hits
delt-hit visualize library --config_path="$CONFIG_PATH" --library_name=top_hits
delt-hit library properties --config_path="$CONFIG_PATH" --library_name=top_hits
```

This enumerates the top combinations by observed count. For full-library processing, omit the count filter:

```bash
delt-hit library enumerate --config_path="$CONFIG_PATH"
delt-hit library properties --config_path="$CONFIG_PATH"
delt-hit library represent --method=morgan --config_path="$CONFIG_PATH"
```

### Compare selections and rank hits

Create a **separate `analysis.yaml`** using the [example analysis configuration](supporting_material/experiments/example-single-display/analysis.yaml). Each entry in its `experiments` list defines a named comparison. Its `selections` list provides each selection's `name`, `counts_path`, and `group` (`condition` or `control`). Replace the example selections and paths with your own, retaining replicate groupings.

```bash
ANALYSIS_CONFIG_PATH=/path/to/analysis.yaml
ANALYSIS_OUTPUT_ROOT="$PROJECT_DIR/analysis"

delt-hit analyse enrichment \
  --analysis_config="$ANALYSIS_CONFIG_PATH" \
  --name=condition_vs_control --method=counts \
  --save_dir="$ANALYSIS_OUTPUT_ROOT"
Rscript --vanilla "$ANALYSIS_OUTPUT_ROOT/counts/condition_vs_control/enrichment_counts.R"

delt-hit analyse enrichment \
  --analysis_config="$ANALYSIS_CONFIG_PATH" \
  --name=condition_vs_control --method=edgeR \
  --save_dir="$ANALYSIS_OUTPUT_ROOT"
Rscript --vanilla "$ANALYSIS_OUTPUT_ROOT/edgeR/condition_vs_control/enrichment_edgeR.R"
```

`--name` must match a comparison in `analysis.yaml`. The CLI generates analysis inputs and an R script; the `Rscript` command performs the analysis. Counts and edgeR compare replicate groups. Normalized z-scores can also be calculated for an individual selection without replicates:

```bash
delt-hit analyse enrichment \
  --config_path="$CONFIG_PATH" \
  --counts="$PROJECT_DIR/selections/$SELECTION_NAME/counts.txt" \
  --method=z_score --name="$SELECTION_NAME" \
  --save_dir="$ANALYSIS_OUTPUT_ROOT"
Rscript --vanilla "$ANALYSIS_OUTPUT_ROOT/z_score/$SELECTION_NAME/enrichment_z_score.R"
```

### Reproduction instructions (optional)

For instructions to reproduce the published-dataset reanalyses reported in the manuscript, including data downloads, workflow commands, and comparisons with the original selection counts, see:

- [Favalli et al. re-analysis](supporting_material/experiments/favalli/README.md)
- [Pure-DEL (Keller et al.) re-analysis](supporting_material/experiments/pure-del/README.md)

### Documentation

- [CLI input/output reference](documentation/input-output.md)
- [Codebase overview](documentation/overview.md)
- [Supporting material and published-data reanalyses](supporting_material/README.md)
- [Demultiplexing benchmarks](benchmarks/demultiplex/README.md)
- [Archived workflow data and outputs](https://doi.org/10.5281/zenodo.20447074)
- [Archived code release](https://doi.org/10.5281/zenodo.20556531)

For command-specific help, use `delt-hit --help` or `delt-hit <group> <command> --help`.

## License

DELT-Hit is distributed under the [MIT License](LICENSE), which permits use, modification, and redistribution subject to its terms. Source code is available in the [DELTechnology/delt-hit repository](https://github.com/DELTechnology/delt-hit).
