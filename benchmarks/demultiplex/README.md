# Synthetic Demultiplex Benchmarking

## Prerequisites

The benchmark commands below assume these environments already exist:

- DELT-Hit: [`.venv`](../../.venv)
- DELi: [`other_tools/DELi/.venv`](../../other_tools/DELi/.venv)

## Dataset Matrix

This benchmark runbook has two synthetic FASTQ benchmark families:

- the canonical matrix with `10` building blocks per cycle for `2`, `3`, and `4` cycles
- an extended `2`-cycle matrix with `100` and `1000` building blocks per cycle

The canonical dataset names are:

- `synthetic_2cycle_10bbpc_1m`
- `synthetic_2cycle_10bbpc_10m`
- `synthetic_2cycle_10bbpc_100m`
- `synthetic_2cycle_10bbpc_1000m`
- `synthetic_3cycle_10bbpc_1m`
- `synthetic_3cycle_10bbpc_10m`
- `synthetic_3cycle_10bbpc_100m`
- `synthetic_3cycle_10bbpc_1000m`
- `synthetic_4cycle_10bbpc_1m`
- `synthetic_4cycle_10bbpc_10m`
- `synthetic_4cycle_10bbpc_100m`
- `synthetic_4cycle_10bbpc_1000m`

The extended `2`-cycle dataset names are:

- `synthetic_2cycle_100bbpc_1m`
- `synthetic_2cycle_100bbpc_10m`
- `synthetic_2cycle_100bbpc_100m`
- `synthetic_2cycle_100bbpc_1000m`
- `synthetic_2cycle_1000bbpc_1m`
- `synthetic_2cycle_1000bbpc_10m`
- `synthetic_2cycle_1000bbpc_100m`
- `synthetic_2cycle_1000bbpc_1000m`

Total reads follow:

```text
building_blocks_per_cycle ** num_cycles * num_reads_per_compound
```

The exact `--num-reads-per-compound` values for the target read depths are:

- `2` cycles with `10` BB/cycle: `1m=10000`, `10m=100000`, `100m=1000000`, `1000m=10000000`
- `2` cycles with `100` BB/cycle: `1m=100`, `10m=1000`, `100m=10000`, `1000m=100000`
- `2` cycles with `1000` BB/cycle: `1m=1`, `10m=10`, `100m=100`, `1000m=1000`
- `3` cycles: `1m=1000`, `10m=10000`, `100m=100000`, `1000m=1000000`
- `4` cycles: `1m=100`, `10m=1000`, `100m=10000`, `1000m=100000`

## 1. Single-Dataset Example

The canonical generator is [`generate_synthetic_fastq.py`](./generate_synthetic_fastq.py).

Each dataset is written under `benchmarks/demultiplex/data/<dataset>_err=<num-errors>/` and includes:

- `<dataset>.fastq.gz`
- `building_blocks.tsv`
- `expected_counts.tsv`
- `manifest.json`

Create one synthetic dataset as a concrete example:

```bash
CYCLES=4
BBPC=10
DEPTH=100m
READS=10000
ERR=0
DATASET="synthetic_${CYCLES}cycle_${BBPC}bbpc_${DEPTH}_err=${ERR}"

./.venv/bin/python benchmarks/demultiplex/generate_synthetic_fastq.py \
  --num-cycles "$CYCLES" \
  --building-blocks-per-cycle "$BBPC" \
  --num-reads-per-compound "$READS" \
  --num-errors "$ERR" \
  --output-dir benchmarks/demultiplex/data \
  --dataset-name "$DATASET"
```

Prepare DELi inputs for that dataset:

```bash
./.venv/bin/python benchmarks/demultiplex/converter/create_deli_inputs.py \
  --dataset-name "$DATASET"
```

Prepare DELT-Hit inputs for that dataset:

```bash
./.venv/bin/python benchmarks/demultiplex/converter/create_delt_inputs.py \
  --dataset-name "$DATASET" \
  --num-cores 11
```

Run the benchmark for that dataset:

```bash
./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
  --dataset-name "$DATASET" \
  --tool both
```

## 2. Prepare DELi And DELT-Hit Inputs

Use the converter scripts in [`benchmarks/demultiplex/converter`](./converter).

### DELi

```bash
./.venv/bin/python benchmarks/demultiplex/converter/create_deli_inputs.py \
  --dataset-name synthetic_4cycle_10bbpc_100m_err=0
```

Prepared inputs are written to:

- repo-local `benchmarks/demultiplex/tools/deli/<dataset>/` for prepared DELi inputs and `timing.json`
- scratch-local `$DATA_DIR/tools/deli/<dataset>/` for generated DELi inputs during node execution

### DELT-Hit

```bash
./.venv/bin/python benchmarks/demultiplex/converter/create_delt_inputs.py \
  --dataset-name synthetic_4cycle_10bbpc_100m_err=0 \
  --num-cores 11
```

Prepared inputs are written to:

- repo-local `benchmarks/demultiplex/tools/delt/<dataset>/` for prepared DELT-Hit inputs and `timing.json`
- scratch-local `$DATA_DIR/tools/delt/<dataset>/` for generated DELT-Hit inputs during node execution

Use `--num-cores 11` when preparing DELT-Hit inputs inside a Slurm job with `12` CPUs so one core remains outside the tool's configured worker count.

## 3. Run Split-Timing Benchmarks

Run [`benchmarks/demultiplex/run_split_timing.py`](./run_split_timing.py) against one prepared dataset.

The runner supports:

- `--tool deli`
- `--tool delt`
- `--tool both`

With `--data-dir "$DATA_DIR"`, it writes runtime logs and intermediate artifacts to:

- `$DATA_DIR/runtime/<dataset>/deli`
- `$DATA_DIR/runtime/<dataset>/delt`

It writes the canonical machine-readable report to:

- `benchmarks/demultiplex/tools/deli/<dataset>/timing.json`
- `benchmarks/demultiplex/tools/delt/<dataset>/timing.json`

### DELi Only

```bash
./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
  --dataset-name synthetic_4cycle_10bbpc_100m_err=0 \
  --tool deli
```

### DELT-Hit Only

```bash
./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
  --dataset-name synthetic_4cycle_10bbpc_100m_err=0 \
  --tool delt
```

## 4. Full Experiment Matrix

For the full experiment matrix, define the dataset families once and reuse them in the Slurm submission loops below:

```bash
CANONICAL_DATASETS=(
  synthetic_2cycle_10bbpc_1m
  synthetic_2cycle_10bbpc_10m
  synthetic_2cycle_10bbpc_100m
  synthetic_2cycle_10bbpc_1000m
  synthetic_3cycle_10bbpc_1m
  synthetic_3cycle_10bbpc_10m
  synthetic_3cycle_10bbpc_100m
  synthetic_3cycle_10bbpc_1000m
  synthetic_4cycle_10bbpc_1m
  synthetic_4cycle_10bbpc_10m
  synthetic_4cycle_10bbpc_100m
  synthetic_4cycle_10bbpc_1000m
)

LARGE_2CYCLE_DATASETS=(
  synthetic_2cycle_100bbpc_1m
  synthetic_2cycle_100bbpc_10m
  synthetic_2cycle_100bbpc_100m
  synthetic_2cycle_100bbpc_1000m
  synthetic_2cycle_1000bbpc_1m
  synthetic_2cycle_1000bbpc_10m
  synthetic_2cycle_1000bbpc_100m
  synthetic_2cycle_1000bbpc_1000m
)

ALL_DATASETS=(
  "${CANONICAL_DATASETS[@]}"
  "${LARGE_2CYCLE_DATASETS[@]}"
)

# When you want one job per dataset and per tool, keep the same arrays and submit them from a login node. 
# These examples require `$TMPDIR` to be set, and each job creates its own synthetic dataset there so jobs do not race on shared files.

dataset_time_limit() {
  case "$1" in
    *_1m) echo "04:00:00" ;;
    *_10m) echo "08:00:00" ;;
    *_100m) echo "12:00:00" ;;
    *_1000m) echo "24:00:00" ;;
    *) echo "Unknown dataset size for $1" >&2; return 1 ;;
  esac
}

dataset_spec() {
  case "$1" in
    synthetic_2cycle_10bbpc_1m) echo "2 10 10000" ;;
    synthetic_2cycle_10bbpc_10m) echo "2 10 100000" ;;
    synthetic_2cycle_10bbpc_100m) echo "2 10 1000000" ;;
    synthetic_2cycle_10bbpc_1000m) echo "2 10 10000000" ;;
    synthetic_3cycle_10bbpc_1m) echo "3 10 1000" ;;
    synthetic_3cycle_10bbpc_10m) echo "3 10 10000" ;;
    synthetic_3cycle_10bbpc_100m) echo "3 10 100000" ;;
    synthetic_3cycle_10bbpc_1000m) echo "3 10 1000000" ;;
    synthetic_4cycle_10bbpc_1m) echo "4 10 100" ;;
    synthetic_4cycle_10bbpc_10m) echo "4 10 1000" ;;
    synthetic_4cycle_10bbpc_100m) echo "4 10 10000" ;;
    synthetic_4cycle_10bbpc_1000m) echo "4 10 100000" ;;
    synthetic_2cycle_100bbpc_1m) echo "2 100 100" ;;
    synthetic_2cycle_100bbpc_10m) echo "2 100 1000" ;;
    synthetic_2cycle_100bbpc_100m) echo "2 100 10000" ;;
    synthetic_2cycle_100bbpc_1000m) echo "2 100 100000" ;;
    synthetic_2cycle_1000bbpc_1m) echo "2 1000 1" ;;
    synthetic_2cycle_1000bbpc_10m) echo "2 1000 10" ;;
    synthetic_2cycle_1000bbpc_100m) echo "2 1000 100" ;;
    synthetic_2cycle_1000bbpc_1000m) echo "2 1000 1000" ;;
    *) echo "Unsupported dataset: $1" >&2; return 1 ;;
  esac
}

for dataset in "${ALL_DATASETS[@]}"; do
  time_limit="$(dataset_time_limit "$dataset")"
  read -r cycles bbpc reads < <(dataset_spec "$dataset")
  err=0
  dataset="${dataset}_err=${err}"

  for tool in deli delt; do
    sbatch \
      --job-name="${tool}-${dataset}" \
      --time="$time_limit" \
      --cpus-per-task=12 \
      --mem=128G \
      --output="$HOME/logs/slurm-%j.out" \
      --export=ALL,DATASET="$dataset",ERR="$err",CYCLES="$cycles",BBPC="$bbpc",READS="$reads",TOOL="$tool" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

cd /users/amarti51/projects/delt-hit

TMPDIR="${TMPDIR:?TMPDIR must be set for Slurm benchmark jobs}"
DATA_ROOT="$TMPDIR/delt-hit-benchmarks/${SLURM_JOB_ID}/${DATASET}"
mkdir -p "$DATA_ROOT"

pixi run python benchmarks/demultiplex/generate_synthetic_fastq.py \
  --num-cycles "$CYCLES" \
  --building-blocks-per-cycle "$BBPC" \
  --num-reads-per-compound "$READS" \
  --num-errors "$ERR" \
  --output-dir "$DATA_ROOT/data" \
  --dataset-name "$DATASET"

pixi run python benchmarks/demultiplex/converter/create_deli_inputs.py \
  --dataset-name "$DATASET" \
  --data-dir "$DATA_ROOT"

pixi run python benchmarks/demultiplex/converter/create_delt_inputs.py \
  --dataset-name "$DATASET" \
  --data-dir "$DATA_ROOT" \
  --num-cores 11

pixi run python benchmarks/demultiplex/run_split_timing.py \
  --dataset-name "$DATASET" \
  --data-dir "$DATA_ROOT" \
  --tool "$TOOL"
EOF
  done
done
```

## What `run_split_timing.py` Measures

For DELi:

- `decode run`
- `decode collect`
- `decode count`

For DELT-Hit:

- execution of the prepared `demultiplex.sh`
- `delt_hit demultiplex process`

Each run also compares the observed output against `expected_counts.tsv` and records whether counts match.

## Plot The Current Benchmark Configurations

After the timing runs are available under `benchmarks/demultiplex/tools`, generate the current runtime plots with:

```bash
uv run python benchmarks/demultiplex/plot_benchmark_runtimes.py

uv run python  benchmarks/demultiplex/plot_2cycle_bbpc_comparison.py
```

This writes:

- `benchmarks/demultiplex/plots/synthetic_2cycle_10bbpc_runtime.png`
- `benchmarks/demultiplex/plots/synthetic_3cycle_10bbpc_runtime.png`
- `benchmarks/demultiplex/plots/synthetic_4cycle_10bbpc_runtime.png`
- `benchmarks/demultiplex/plots/synthetic_2cycle_bbpc_comparison_runtime.png`

## Manual Command Sequence

If you want to run the prepared inputs manually instead of using `run_split_timing.py`, the equivalent commands are:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
DATASET="synthetic_4cycle_10bbpc_100m_err=0"

DELI_DIR="$ROOT/benchmarks/demultiplex/tools/deli/$DATASET"
DELT_DIR="$ROOT/benchmarks/demultiplex/tools/delt/$DATASET"
DELI_OUT="$ROOT/target/benchmarks/$DATASET/split_timing/deli"
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
