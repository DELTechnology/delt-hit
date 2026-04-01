# DELT-Hit

DELT-Hit is an open-source Python CLI for end-to-end analysis of DNA-encoded chemical library (DEL) screening data. Starting from raw FASTQ reads and a library design spreadsheet, it demultiplexes sequencing data, enumerates combinatorial compound libraries, calculates enrichment statistics, and produces ranked hit lists with molecular properties and fingerprint representations. Developed at ETH Zurich and released under the MIT License.

## Quick Start

**Install from GitHub (users):**
```bash
pip install git+https://github.com/DELTechnology/delt-hit.git
```

> RDKit requires Python 3.12+. If `pip install rdkit` fails on your platform, install it via conda first:
> ```bash
> conda install -c conda-forge rdkit
> pip install git+https://github.com/DELTechnology/delt-hit.git --no-deps
> ```

**Development install:**
```bash
git clone git@github.com:DELTechnology/delt-hit.git
cd delt-hit
uv sync                       # installs all core dependencies into .venv/
source .venv/bin/activate
delt-hit --help
```

## Pipeline Overview

The pipeline is organised into five modules that run in sequence:

| Module | Command | Description |
|---|---|---|
| **Init** | `delt-hit init` | Parse an Excel workbook and generate a `config.yaml` |
| **Demultiplex** | `delt-hit demultiplex` | Run cutadapt to assign reads to selections; count barcodes per compound |
| **Library** | `delt-hit library` | Enumerate SMILES via RDKit reactions; compute properties and fingerprints |
| **Analyse** | `delt-hit analyse` | Generate enrichment R scripts (edgeR or simple counts); rank hits |
| **Dashboard** | `delt-hit dashboard` | Launch an interactive Dash web app to explore counts |

For details on each command's flags, run `delt-hit <module> --help`.

## Running the Full Pipeline

### Manual (step by step)

```bash
# 1. Generate config from Excel
delt-hit init --excel_path templates/library.xlsx

# 2. Demultiplex FASTQ and compute counts
delt-hit demultiplex run     --config_path experiments/EXPNAME/config.yaml
delt-hit demultiplex process --config_path experiments/EXPNAME/config.yaml
delt-hit demultiplex report  --config_path experiments/EXPNAME/config.yaml
delt-hit demultiplex qc      --config_path experiments/EXPNAME/config.yaml

# 3. Enumerate library and compute properties (independent of demultiplex)
delt-hit library enumerate   --config_path experiments/EXPNAME/config.yaml
delt-hit library properties  --config_path experiments/EXPNAME/config.yaml
delt-hit library represent   --config_path experiments/EXPNAME/config.yaml --method=morgan

# 4. Generate enrichment R script, then execute it
delt-hit analyse enrichment  --config_path analysis.yaml --name=protein_vs_no_protein --method=counts
Rscript experiments/EXPNAME/counts/enrichment_counts.R

# 5. Explore results interactively
delt-hit dashboard \
  --config_path experiments/EXPNAME/config.yaml \
  --counts_path experiments/EXPNAME/selections/SELECTION_NAME/counts.txt
```

### Snakemake (recommended)

Snakemake resolves step dependencies automatically and runs independent branches in parallel.

```bash
uv sync --extra workflow          # install snakemake into the project venv
snakemake -n                      # dry-run: show what would execute
snakemake --cores 4               # run full pipeline locally
snakemake --dag | dot -Tpng > dag.png   # visualise the DAG

# On a SLURM cluster:
snakemake --profile configs/snakemake/slurm/
```

Edit `configs/snakemake/slurm/config.yaml` to set your cluster partition before submitting.

### Docker

```bash
docker build -t delt-hit:dev .
docker run --rm \
  -v $(pwd)/experiments:/workspace/experiments \
  -v $(pwd)/config.yaml:/workspace/config.yaml \
  delt-hit:dev demultiplex run --config_path=config.yaml
```

## Development

**Set up the environment:**
```bash
uv sync --extra dev     # installs pytest, ruff, mypy, pre-commit
source .venv/bin/activate
pre-commit install      # wire git hooks
```

**Run tests:**
```bash
pytest tests/ -v
```

**Lint:**
```bash
ruff check src/
ruff format src/        # auto-format
```

**Type check:**
```bash
mypy src/delt_hit/
```

**Run a single pipeline step during development:**
```bash
uv run delt-hit demultiplex prepare --config_path=config.yaml --fast_dev_run=true
```

## Documentation

- [Pipeline architecture](docs/references/analysis_pipeline_architecture.md)
- [CLI reference](documentation/cli.md)
- [Codebase overview](documentation/overview.md)

## License

MIT — see [LICENSE](LICENSE).
