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

delt-hit library properties \
  --config_path <path/to/config.yaml> \
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

Useful options:
- `--library_name` to compute properties for a named library parquet such as filtered enumeration output
- `--library_path` to compute properties from an explicit parquet path

**Outputs**
- `<save_dir>/<experiment_name>/library/properties/properties.parquet`
- `<save_dir>/<experiment_name>/library/properties/<library_name>.parquet` for named-library mode
- Histogram PNGs per property

## `visualize`
Visualization commands for chemistry assets and workflow outputs.

### `enumerate`
Generates reviewer-friendly input visualizations for enumeration: reaction graphs, 2D reaction scheme panels from SMIRKS, 2D structure grids for building blocks, and configured compound structures.

```
delt-hit visualize enumerate --config_path <path/to/config.yaml>
```

Useful options:
- `--building_block_ids` to restrict the visualization to selected building-block families
- `--nrow` to control how many molecules are shown per row in structure grids
- `--dpi` to control PNG export resolution
- `--tile_size` to control the RDKit tile size used for each rendered molecule
- `--graph`, `--reactions`, `--building_blocks`, and `--compounds` to select specific visualization outputs

If none of the selector flags are passed, all enumeration visualizations are generated.

Example visualization workflow:

```
delt-hit visualize enumerate \
  --config_path <path/to/config.yaml>
```

Example compounds-only visualization:

```
delt-hit visualize enumerate \
  --config_path <path/to/config.yaml> \
  --compounds true
```

**Outputs**
- Reaction graph PNG files in `<save_dir>/<experiment_name>/library/visualization/` named `reaction_graph.png`, `reaction_graph_additional.png`, and `reaction_graph_building_blocks.png`
- Reaction scheme PNG files in `<save_dir>/<experiment_name>/library/visualization/reactions/`
- `building_blocks_<BUILDING_BLOCK_ID>.png` for each visualized building-block family
- Per-compound PNG files in `<save_dir>/<experiment_name>/library/visualization/compounds/`

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
