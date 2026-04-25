# Benchmarks

This directory contains the synthetic benchmark workflow for DELi and DELT-Hit.

The workflow has three stages:

1. Generate a synthetic FASTQ dataset plus truth tables.
2. Convert that dataset into prepared DELi and DELT-Hit inputs.
3. Run both toolchains and collect split timing outputs.

## Prerequisites

The benchmark commands below assume these environments already exist:

- DELT-Hit: `/Users/adrianomartinelli/projects/delt-hit/.venv`
- DELi: `/Users/adrianomartinelli/projects/delt-hit/other_tools/DELi/.venv`

## 1. Generate A Synthetic Dataset

The canonical generator is [`scripts/generate_synthetic_fastq.py`](/Users/adrianomartinelli/projects/delt-hit/scripts/generate_synthetic_fastq.py).

It writes four files under `benchmarks/data/<dataset>/`:

- `<dataset>.fastq.gz`
- `building_blocks.tsv`
- `expected_counts.tsv`
- `manifest.json`

The total reads are:

```text
building_blocks_per_cycle ** num_cycles * num_reads_per_compound
```

Example: generate the 4-cycle 100M-read dataset we created in this repo:

```bash
./.venv/bin/python scripts/generate_synthetic_fastq.py \
  --num-cycles 4 \
  --building-blocks-per-cycle 10 \
  --num-reads-per-compound 10000 \
  --output-dir benchmarks/data \
  --experiment-name synthetic_4cycle_100m
```

This produces:

- `/Users/adrianomartinelli/projects/delt-hit/benchmarks/data/synthetic_4cycle_100m/synthetic_4cycle_100m.fastq.gz`
- `/Users/adrianomartinelli/projects/delt-hit/benchmarks/data/synthetic_4cycle_100m/building_blocks.tsv`
- `/Users/adrianomartinelli/projects/delt-hit/benchmarks/data/synthetic_4cycle_100m/expected_counts.tsv`
- `/Users/adrianomartinelli/projects/delt-hit/benchmarks/data/synthetic_4cycle_100m/manifest.json`

## 2. Prepare DELi And DELT-Hit Inputs

Use the converter scripts in [`benchmarks/converter`](/Users/adrianomartinelli/projects/delt-hit/benchmarks/converter).

For DELi:

```bash
./.venv/bin/python benchmarks/converter/create_deli_inputs.py \
  --dataset-name synthetic_4cycle_100m
```

This writes prepared inputs under:

- `/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/deli/synthetic_4cycle_100m`

Key files:

- `decode_synthetic.yaml`
- `decode_settings_v02.yaml`
- `deli_config`
- `data/libraries/synthetic_4cycle_100m_library.json`

For DELT-Hit:

```bash
./.venv/bin/python benchmarks/converter/create_delt_inputs.py \
  --dataset-name synthetic_4cycle_100m
```

This writes prepared inputs under:

- `/Users/adrianomartinelli/projects/delt-hit/benchmarks/tools/delt/synthetic_4cycle_100m`

Key files:

- `config.yaml`
- `synthetic_4cycle_100m/demultiplex/cutadapt_input_files/demultiplex.sh`

`create_delt_inputs.py` also runs `delt_hit demultiplex prepare`, so the `demultiplex.sh` script is ready to execute immediately.

## 3. Run The Benchmark

To reproduce the same style of output as:

- `/Users/adrianomartinelli/projects/delt-hit/target/benchmarks/synthetic_3cycle_50m/split_timing`

run [`benchmarks/run_split_timing.py`](/Users/adrianomartinelli/projects/delt-hit/benchmarks/run_split_timing.py):

```bash
./.venv/bin/python benchmarks/run_split_timing.py \
  --dataset-name synthetic_4cycle_100m
```

By default this writes to:

- `/Users/adrianomartinelli/projects/delt-hit/target/benchmarks/synthetic_4cycle_100m/split_timing`

That output directory contains:

- `deli/`
- `delt/`
- `summary.json`

Within those tool-specific directories the runner stores:

- command logs
- temporary cache directories
- DELi decode/count outputs
- DELT-Hit demultiplex and count outputs

## What `run_split_timing.py` Measures

For DELi:

- `decode run`
- `decode collect`
- `decode count`

For DELT-Hit:

- execution of the prepared `demultiplex.sh`
- `delt_hit demultiplex process`

The runner also compares both observed outputs against `expected_counts.tsv` and records whether counts match.

## Manual Command Sequence

If you want to run the prepared inputs manually instead of using `run_split_timing.py`, the equivalent commands are:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/adrianomartinelli/projects/delt-hit"
DATASET="synthetic_4cycle_100m"

DELI_DIR="$ROOT/benchmarks/tools/deli/$DATASET"
DELT_DIR="$ROOT/benchmarks/tools/delt/$DATASET"
DELI_OUT="$DELI_DIR/output_run"
DELT_EXP_DIR="$DELT_DIR/$DATASET"

mkdir -p "$DELI_OUT"
mkdir -p "$ROOT/.tmp/mpl" "$ROOT/.tmp/fontconfig"

export MPLCONFIGDIR="$ROOT/.tmp/mpl"
export XDG_CACHE_HOME="$ROOT/.tmp"
export FC_CACHEDIR="$ROOT/.tmp/fontconfig"
export PATH="$ROOT/.venv/bin:$ROOT/other_tools/DELi/.venv/bin:$PATH"

deli \
  --config-file "$DELI_DIR/deli_config" \
  decode run \
  "$DELI_DIR/decode_synthetic.yaml" \
  --decode-settings-file "$DELI_DIR/decode_settings_v02.yaml" \
  --out-dir "$DELI_OUT" \
  --prefix "$DATASET" \
  --skip-report

deli \
  --config-file "$DELI_DIR/deli_config" \
  decode collect \
  "$DELI_OUT/${DATASET}_decoded.tsv" \
  --out-loc "$DELI_OUT/${DATASET}_collected.ndjson"

deli \
  --config-file "$DELI_DIR/deli_config" \
  decode count \
  "$DELI_OUT/${DATASET}_collected.ndjson" \
  --out-loc "$DELI_OUT/${DATASET}_counts.parquet" \
  --output-format parquet

bash -e \
  "$DELT_EXP_DIR/demultiplex/cutadapt_input_files/demultiplex.sh"

python -m delt_hit.cli.main \
  demultiplex process \
  --config_path "$DELT_DIR/config.yaml"
```

Use the Python runner when you want the standardized `split_timing` directory layout and `summary.json`.
