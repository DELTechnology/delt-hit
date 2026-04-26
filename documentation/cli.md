# DELT-Hit CLI guide

The CLI is the main entrypoint for running DELT-Hit workflows. It is exposed via the `delt-hit` console script and uses `jsonargparse` to map subcommands to public methods in the CLI classes.

```
delt-hit <group> <method> [--args]
```

The available groups are:
- `init`
- `library`
- `demultiplex`
- `analyse`
- `dashboard`

> Tip: run `delt-hit --help` or `delt-hit <group> --help` to see argument details.

## Configuration inputs
Most commands expect a YAML configuration file. The standard way to create it is:

```
delt-hit init --excel_path <path/to/library.xlsx>
```

The resulting `config.yaml` contains:
- `experiment`: name, save directory, input FASTQ path, CPU cores.
- `selections`: metadata for each selection plus primer identifiers.
- `library`: reaction graph edges and building-block definitions.
- `catalog`: reaction SMARTS and compound definitions.
- `structure`: parsing structure (selection/building block/constant regions).
- `whitelists`: codon lists derived from selections, building blocks, and constants.

The configuration layout is derived directly from the Excel template sheets (see `templates/library.xlsx`) and is parsed by `delt_hit.demultiplex.parser`.

## `init`
Creates a YAML config from an Excel template.

```
delt-hit init --excel_path <path/to/library.xlsx>
```

**Outputs**
- `<save_dir>/<experiment_name>/config.yaml`

## `demultiplex`
Provides FASTQ demultiplexing and post-processing of Cutadapt outputs.

### `prepare`
Generates Cutadapt input files and an executable shell script.

```
delt-hit demultiplex prepare --config_path <path/to/config.yaml>
```

**Outputs**
- `<save_dir>/<experiment_name>/demultiplex/cutadapt_input_files/`
- `demultiplex.sh` shell script that chains Cutadapt steps
- FASTQ barcode files per region (`S*`, `B*`, etc.)

### `run`
Runs the demultiplexing pipeline end-to-end by generating the script and executing it.

```
delt-hit demultiplex run --config_path <path/to/config.yaml>
```

### `process`
Consumes Cutadapt output and computes per-selection barcode counts.

```
delt-hit demultiplex process --config_path <path/to/config.yaml>
```

**Outputs**
- `<save_dir>/<experiment_name>/selections/<SELECTION_NAME>/counts.txt`
  - `code_0`, `code_1`, … columns plus `count`

### `report`
Builds a text summary of Cutadapt statistics.

```
delt-hit demultiplex report --config_path <path/to/config.yaml>
```

**Outputs**
- `<save_dir>/<experiment_name>/qc/report.txt`

### `qc`
Generates QC plots from demultiplexed counts.

```
delt-hit demultiplex qc --config_path <path/to/config.yaml>
```

**Outputs**
- `<save_dir>/<experiment_name>/qc/` (plots)

## `library`
Library and descriptor generation for downstream analysis.

### `enumerate`
Builds the reaction graph, enumerates building block combinations, and generates SMILES.

```
delt-hit library enumerate --config_path <path/to/config.yaml>
```

Useful options:
- `--overwrite` to re-generate an existing library
- `--graph_only` to skip enumeration and only write reaction graph visualizations
- `--building_block_ids` to enumerate a subset of building blocks
- `--counts_path`, `--top_n`, and `--library_name` to enumerate only the top observed combinations from a demultiplex `counts.txt` file

Example filtered-enumeration workflow:

```
delt-hit demultiplex process --config_path <path/to/config.yaml>

delt-hit library enumerate \
  --config_path <path/to/config.yaml> \
  --counts_path <save_dir>/<experiment_name>/selections/<SELECTION_NAME>/counts.txt \
  --top_n 1000 \
  --library_name observed_hits
```

**Outputs**
- `<save_dir>/<experiment_name>/library/library.parquet`
- `<save_dir>/<experiment_name>/library/<library_name>.parquet` for filtered enumeration
- Reaction graph PNGs in the same directory

### `properties`
Computes molecular properties (RDKit/chem-informatics descriptors) and plots their distributions.

```
delt-hit library properties --config_path <path/to/config.yaml>
```

**Outputs**
- `<save_dir>/<experiment_name>/library/properties/properties.parquet`
- Histogram PNGs per property

### `visualize`
Generates reviewer-friendly chemical visualizations: reaction graphs annotated with SMIRKS, 2D reaction scheme panels, and 2D structure grids for products and selected building blocks.

```
delt-hit library visualize --config_path <path/to/config.yaml>
```

Useful options:
- `--counts_path`, `--top_n`, and `--library_name` to visualize the top observed combinations from a demultiplex `counts.txt` file
- `--output_name` to control the PNG filename prefix
- `--library_path` to render structures from an existing parquet library

Example top-hit visualization workflow:

```
delt-hit library visualize \
  --config_path <path/to/config.yaml> \
  --counts_path <save_dir>/<experiment_name>/selections/<SELECTION_NAME>/counts.txt \
  --top_n 24 \
  --library_name observed_hits \
  --output_name reviewer_panel
```

**Outputs**
- Reaction graph PNGs in `<save_dir>/<experiment_name>/library/`
- `<output_name>_reaction_schemes.png`
- `<output_name>_products.png`
- `<output_name>_<BUILDING_BLOCK_ID>.png` for each visualized building-block family

### `represent`
Generates machine-learning representations (fingerprints).

```
delt-hit library represent --config_path <path/to/config.yaml> --method morgan
```

Supported methods:
- `morgan` (stored as a SciPy sparse matrix)
- `bert` (currently routed through the Morgan generator in the CLI wrapper; see implementation)

**Outputs**
- `<save_dir>/<experiment_name>/representations/<method>.npz`

## `analyse`
Statistical analysis over per-selection counts. The analysis config expects an `experiments` list with explicit selection entries and `counts_path` values (see `delt_hit.cli.analyse.api.prepare_data`).

### `enrichment`
Runs count-based or edgeR-based enrichment analysis.

```
delt-hit analyse enrichment --config_path <path/to/config.yaml> --name <experiment-name> --method edgeR
```

**Outputs**
- `<save_dir>/<experiment_name>/edgeR/` (or `counts/`) with statistics tables and normalized counts.

## `dashboard`
Launches a Dash web UI for inspecting a single selection’s counts.

```
delt-hit dashboard --config_path <path/to/config.yaml> --counts_path <path/to/selections/SELECTION_NAME/counts.txt>
```

The dashboard defaults to port `8050` and displays:
- experiment metadata
- selection metadata
- interactive scatter plots for code combinations

## Where to look in the code
- CLI wiring: `src/delt_hit/cli/main.py`
- CLI group implementations: `src/delt_hit/cli/{init,demultiplex,library,analyse,dashboard}`
- Config parsing: `src/delt_hit/demultiplex/parser.py`
- Demultiplexing pipeline: `src/delt_hit/demultiplex/preprocess.py` and `src/delt_hit/demultiplex/postprocess.py`
