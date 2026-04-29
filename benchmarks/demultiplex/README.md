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

## 1. Generate Synthetic Datasets

The canonical generator is [`generate_synthetic_fastq.py`](./generate_synthetic_fastq.py).

Each dataset is written under `benchmarks/demultiplex/data/<dataset>/` and includes:

- `<dataset>.fastq.gz`
- `building_blocks.tsv`
- `expected_counts.tsv`
- `manifest.json`

Example: generate `synthetic_4cycle_10bbpc_100m`:

```bash
./.venv/bin/python benchmarks/demultiplex/generate_synthetic_fastq.py \
  --num-cycles 4 \
  --building-blocks-per-cycle 10 \
  --num-reads-per-compound 10 \
  --num-errors 1 \
  --output-dir benchmarks/demultiplex/data \
  --experiment-name synthetic_4cycle_10bbpc_100m
```

Generate the canonical `10`-BB/cycle matrix:

```bash
./.venv/bin/python benchmarks/demultiplex/generate_synthetic_benchmark_matrix.py \
  --profile canonical \
  --num-errors 0 \
  --output-dir benchmarks/demultiplex/data
```

Generate the extended `2`-cycle large-library matrix:

```bash
./.venv/bin/python benchmarks/demultiplex/generate_synthetic_benchmark_matrix.py \
  --profile two_cycle_large_libraries \
  --num-errors 0 \
  --output-dir benchmarks/demultiplex/data
```

Generate every benchmark dataset in the combined matrix:

```bash
./.venv/bin/python benchmarks/demultiplex/generate_synthetic_benchmark_matrix.py \
  --profile all \
  --num-errors 0 \
  --output-dir benchmarks/demultiplex/data
```

Preview the full plan without writing files:

```bash
./.venv/bin/python benchmarks/demultiplex/generate_synthetic_benchmark_matrix.py \
  --profile all \
  --dry-run
```

## 2. Prepare DELi And DELT-Hit Inputs

Use the converter scripts in [`benchmarks/demultiplex/converter`](./converter).

### DELi

```bash
./.venv/bin/python benchmarks/demultiplex/converter/create_deli_inputs.py \
  --dataset-name synthetic_4cycle_10bbpc_100m
```

Prepared inputs are written to:

- repo-local `benchmarks/demultiplex/tools/deli/<dataset>/` for prepared DELi inputs and `timing.json`
- scratch-local `$DATA_DIR/tools/deli/<dataset>/` for generated DELi inputs during node execution

### DELT-Hit

```bash
./.venv/bin/python benchmarks/demultiplex/converter/create_delt_inputs.py \
  --dataset-name synthetic_4cycle_10bbpc_100m \
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
  --dataset-name synthetic_4cycle_10bbpc_100m \
  --tool deli
```

### DELT-Hit Only

```bash
./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
  --dataset-name synthetic_4cycle_10bbpc_100m_err=1 \
  --tool delt
```

## 4. Full-Matrix Bash Loops

For ad hoc local runs, define the full matrix once and reuse it across generation, preparation, and timing steps:

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
```

Generate the full dataset matrix without Slurm:

```bash
for dataset in "${ALL_DATASETS[@]}"; do
  echo "Generating $dataset"
done

./.venv/bin/python benchmarks/demultiplex/generate_synthetic_benchmark_matrix.py \
  --profile all \
  --num-errors 0 \
  --output-dir benchmarks/demultiplex/data
```

Prepare all DELi inputs without Slurm:

```bash
for dataset in "${ALL_DATASETS[@]}"; do
  ./.venv/bin/python benchmarks/demultiplex/converter/create_deli_inputs.py \
    --dataset-name "$dataset"
done
```

Prepare all DELT-Hit inputs without Slurm:

```bash
for dataset in "${ALL_DATASETS[@]}"; do
  ./.venv/bin/python benchmarks/demultiplex/converter/create_delt_inputs.py \
    --dataset-name "$dataset" \
    --num-cores 11
done
```

Run timings for all datasets without Slurm:

```bash
for dataset in "${ALL_DATASETS[@]}"; do
  ./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
    --dataset-name "$dataset" \
    --tool both
done
```

Rerun just the canonical `10bbpc` family:

```bash
for dataset in "${CANONICAL_DATASETS[@]}"; do
  ./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
    --dataset-name "$dataset" \
    --tool both
done
```

## 5. Slurm / `sbatch` Loops

When you want one job per dataset, keep the same arrays and submit them from a login node:

```bash
for dataset in "${ALL_DATASETS[@]}"; do
  case "$dataset" in
    *_1m) time_limit="04:00:00" ;;
    *_10m) time_limit="08:00:00" ;;
    *_100m) time_limit="12:00:00" ;;
    *_1000m) time_limit="24:00:00" ;;
    *) echo "Unknown dataset size for $dataset" >&2; exit 1 ;;
  esac

  sbatch --job-name="bench-${dataset}" --time="$time_limit" --export=ALL,DATASET="$dataset" <<'EOF'
#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --output=slurm-%x-%j.out
set -euo pipefail

cd /path/to/delt-hit

./.venv/bin/python benchmarks/demultiplex/converter/create_deli_inputs.py \
  --dataset-name "$DATASET"

./.venv/bin/python benchmarks/demultiplex/converter/create_delt_inputs.py \
  --dataset-name "$DATASET" \
  --num-cores 11

./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
  --dataset-name "$DATASET" \
  --tool both
EOF
done
```

If the datasets are already generated and prepared, submit timing-only jobs:

```bash
for dataset in "${ALL_DATASETS[@]}"; do
  case "$dataset" in
    *_1m) time_limit="04:00:00" ;;
    *_10m) time_limit="08:00:00" ;;
    *_100m) time_limit="12:00:00" ;;
    *_1000m) time_limit="24:00:00" ;;
    *) echo "Unknown dataset size for $dataset" >&2; exit 1 ;;
  esac

  sbatch --job-name="timing-${dataset}" --time="$time_limit" --export=ALL,DATASET="$dataset" <<'EOF'
#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --output=slurm-%x-%j.out
set -euo pipefail

cd /path/to/delt-hit

./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
  --dataset-name "$DATASET" \
  --tool both
EOF
done
```

If you prefer a single batch job that loops over everything on the compute node:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=synthetic-bench-all
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%x-%j.out
set -euo pipefail

cd /path/to/delt-hit

ALL_DATASETS=(
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
  synthetic_2cycle_100bbpc_1m
  synthetic_2cycle_100bbpc_10m
  synthetic_2cycle_100bbpc_100m
  synthetic_2cycle_100bbpc_1000m
  synthetic_2cycle_1000bbpc_1m
  synthetic_2cycle_1000bbpc_10m
  synthetic_2cycle_1000bbpc_100m
  synthetic_2cycle_1000bbpc_1000m
)

for dataset in "${ALL_DATASETS[@]}"; do
  case "$dataset" in
    *_1m) time_limit="04:00:00" ;;
    *_10m) time_limit="08:00:00" ;;
    *_100m) time_limit="12:00:00" ;;
    *_1000m) time_limit="24:00:00" ;;
    *) echo "Unknown dataset size for $dataset" >&2; exit 1 ;;
  esac

  echo "Running $dataset with recommended wall time $time_limit"
  ./.venv/bin/python benchmarks/demultiplex/run_split_timing.py \
    --dataset-name "$dataset" \
    --tool both
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

## Manual Command Sequence

If you want to run the prepared inputs manually instead of using `run_split_timing.py`, the equivalent commands are:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
DATASET="synthetic_4cycle_10bbpc_100m"

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
