# Benchmarks

This directory contains the synthetic benchmark workflow for DELi and DELT-Hit.

The workflow has three stages:

1. Generate a synthetic FASTQ dataset plus truth tables.
2. Convert that dataset into prepared DELi and DELT-Hit inputs.
3. Run one tool or both tools and collect split timing outputs.

There is also a separate counts-first synthetic enrichment benchmark for validating
the `counts` and `edgeR` analysis methods directly without going through FASTQ
generation or demultiplexing.

## Synthetic Enrichment Benchmark

This benchmark exists to compare ranking recovery of known positive controls under
configurable count noise. It starts from per-selection count tables on purpose:
the goal is to isolate enrichment-analysis behavior rather than test sequencing,
adapter trimming, or barcode decoding.

The generator writes:

- six DELT-Hit-compatible count tables:
  - `no_protein_rep1..3`
  - `protein_rep1..3`
- `config.yaml` ready for `delt-hit analyse enrichment`
- `truth.tsv` with seeded positive tiers and target means
- `manifest.json` with the exact generation parameters

Noise in this benchmark means abundance heterogeneity plus replicate-to-replicate
count variation. It does not model sequencing corruption.

### Generate a synthetic enrichment dataset

```bash
./.venv/bin/python benchmarks/generate_synthetic_enrichment.py \
  --building-blocks-per-cycle 100 \
  --positive-low 50 \
  --positive-medium 25 \
  --positive-high 10 \
  --fc-low 1.5 \
  --fc-medium 3.0 \
  --fc-high 8.0 \
  --dispersion 0.20 \
  --library-size-cv 0.10 \
  --seed 13 \
  --output-dir target/synthetic_enrichment \
  --experiment-name synthetic_enrichment_default
```

The default model uses:

- a heterogeneous per-compound background mean
- explicit low / medium / high positive spike-ins
- negative-binomial replicate sampling
- three `protein` replicates and three `no_protein` replicates

`truth.tsv` includes metadata comment lines describing the distribution family and
the exact parameter values used for that run.

### Run the comparison

```bash
./.venv/bin/python benchmarks/compare_synthetic_enrichment.py \
  --dataset-dir target/synthetic_enrichment/synthetic_enrichment_default \
  --top-k 100
```

This runs both existing enrichment paths:

- `counts`
- `edgeR`

and writes a comparison bundle under `<dataset>/comparison/`:

- `comparison.tsv`: one row per compound with truth labels, scores, and ranks
- `comparison.json`: machine-readable metrics summary
- `top_hits.tsv`: side-by-side top-ranked compounds from each method
- `summary.md`: compact human-readable report

### Interpreting the report

The report is relative rather than threshold-based. The main questions are:

- How many seeded positives appear in the top `K` results?
- Does `edgeR` recover more `positive_low` compounds than raw counts?
- Are negatives staying out of the top-ranked region?
- How do the median and mean positive ranks compare across methods?

## Analysis

Minimal run:

```bash
./.venv/bin/python benchmarks/generate_synthetic_enrichment.py \
  --output-dir target/synthetic_enrichment \
  --experiment-name synthetic_enrichment_default

./.venv/bin/python benchmarks/compare_synthetic_enrichment.py \
  --dataset-dir target/synthetic_enrichment/synthetic_enrichment_default \
  --top-k 100
```

Main output files:

- `target/synthetic_enrichment/synthetic_enrichment_default/comparison/summary.md`
- `target/synthetic_enrichment/synthetic_enrichment_default/comparison/comparison.tsv`

## Prerequisites

The benchmark commands below assume these environments already exist:

- DELT-Hit: [`.venv`](../.venv)
- DELi: [`other_tools/DELi/.venv`](../other_tools/DELi/.venv)

## Dataset Matrix

This benchmark runbook has two synthetic FASTQ benchmark families:

- the canonical matrix with `10` building blocks per cycle for `2`, `3`, and `4` cycles
- an extended `2`-cycle matrix with `100` and `1000` building blocks per cycle

The canonical dataset names are:

- `synthetic_2cycle_1m`
- `synthetic_2cycle_10m`
- `synthetic_2cycle_100m`
- `synthetic_2cycle_1000m`
- `synthetic_3cycle_1m`
- `synthetic_3cycle_10m`
- `synthetic_3cycle_100m`
- `synthetic_3cycle_1000m`
- `synthetic_4cycle_1m`
- `synthetic_4cycle_10m`
- `synthetic_4cycle_100m`
- `synthetic_4cycle_1000m`

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

Each dataset is written under `benchmarks/data/<dataset>/` and includes:

- `<dataset>.fastq.gz`
- `building_blocks.tsv`
- `expected_counts.tsv`
- `manifest.json`

Example: generate `synthetic_4cycle_100m`:

```bash
./.venv/bin/python benchmarks/generate_synthetic_fastq.py \
  --num-cycles 4 \
  --building-blocks-per-cycle 10 \
  --num-reads-per-compound 10 \
  --num-errors 1 \
  --output-dir benchmarks/data \
  --experiment-name synthetic_4cycle_100m
```

Generate the canonical `10`-BB/cycle matrix:

```bash
./.venv/bin/python benchmarks/generate_synthetic_benchmark_matrix.py \
  --profile canonical \
  --num-errors 0 \
  --output-dir benchmarks/data
```

Generate the extended `2`-cycle large-library matrix:

```bash
./.venv/bin/python benchmarks/generate_synthetic_benchmark_matrix.py \
  --profile two_cycle_large_libraries \
  --num-errors 0 \
  --output-dir benchmarks/data
```

Preview the full plan without writing files:

```bash
./.venv/bin/python benchmarks/generate_synthetic_benchmark_matrix.py \
  --profile all \
  --dry-run
```

Generate, prepare, and benchmark the full extended `2`-cycle large-library suite locally:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
cd "$ROOT"

./.venv/bin/python benchmarks/generate_synthetic_benchmark_matrix.py \
  --profile two_cycle_large_libraries \
  --num-errors 0 \
  --output-dir benchmarks/data

for dataset in \
  synthetic_2cycle_100bbpc_1m \
  synthetic_2cycle_100bbpc_10m \
  synthetic_2cycle_100bbpc_100m \
  synthetic_2cycle_100bbpc_1000m \
  synthetic_2cycle_1000bbpc_1m \
  synthetic_2cycle_1000bbpc_10m \
  synthetic_2cycle_1000bbpc_100m \
  synthetic_2cycle_1000bbpc_1000m
do
  ./.venv/bin/python benchmarks/converter/create_deli_inputs.py \
    --dataset-name "$dataset"
  ./.venv/bin/python benchmarks/converter/create_delt_inputs.py \
    --dataset-name "$dataset" \
    --num-cores 11
  ./.venv/bin/python benchmarks/run_split_timing.py \
    --dataset-name "$dataset" \
    --tool both
done
```

If you only want to generate data on a compute node, use `$TMPDIR` as the scratch root:

```bash
ROOT="/users/amarti51/projects/delt-hit"
DATA_DIR="$TMPDIR/benchmarks/synthetic_4cycle_100m"
NUM_ERRORS=0

./.venv/bin/python benchmarks/generate_synthetic_fastq.py \
  --num-cycles 4 \
  --building-blocks-per-cycle 10 \
  --num-reads-per-compound 10000 \
  --num-errors "$NUM_ERRORS" \
  --output-dir "$DATA_DIR/data" \
  --experiment-name synthetic_4cycle_100m
```

## 2. Prepare DELi And DELT-Hit Inputs

Use the converter scripts in [`benchmarks/converter`](./converter).

### DELi

Example:

```bash
./.venv/bin/python benchmarks/converter/create_deli_inputs.py \
  --dataset-name synthetic_4cycle_100m
```

Prepared inputs are written to:

- repo-local `benchmarks/tools/deli/<dataset>/timing.json` for final reports
- scratch-local `$DATA_DIR/tools/deli/<dataset>/` for generated DELi inputs during node execution

Prepare DELi inputs for all 12 datasets in-place:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
cd "$ROOT"

for cycles in 2 3 4; do
  for depth in 1m 10m 100m 1000m; do
    dataset="synthetic_${cycles}cycle_${depth}"
    ./.venv/bin/python benchmarks/converter/create_deli_inputs.py \
      --dataset-name "$dataset"
  done
done
```

Prepare DELi inputs for the extended `2`-cycle large-library suite:

```bash
for dataset in \
  synthetic_2cycle_100bbpc_1m \
  synthetic_2cycle_100bbpc_10m \
  synthetic_2cycle_100bbpc_100m \
  synthetic_2cycle_100bbpc_1000m \
  synthetic_2cycle_1000bbpc_1m \
  synthetic_2cycle_1000bbpc_10m \
  synthetic_2cycle_1000bbpc_100m \
  synthetic_2cycle_1000bbpc_1000m
do
  ./.venv/bin/python benchmarks/converter/create_deli_inputs.py \
    --dataset-name "$dataset"
done
```

### DELT-Hit

Example:

```bash
./.venv/bin/python benchmarks/converter/create_delt_inputs.py \
  --dataset-name synthetic_4cycle_100m \
  --num-cores 11
```

Prepared inputs are written to:

- repo-local `benchmarks/tools/delt/<dataset>/timing.json` for final reports
- scratch-local `$DATA_DIR/tools/delt/<dataset>/` for generated DELT-Hit inputs during node execution

Use `--num-cores 11` when preparing DELT-Hit inputs inside a Slurm job with `12` CPUs so one core remains outside the tool's configured worker count.

Prepare DELT-Hit inputs for all 12 datasets in-place:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
cd "$ROOT"

for cycles in 2 3 4; do
  for depth in 1m 10m 100m 1000m; do
    dataset="synthetic_${cycles}cycle_${depth}"
    ./.venv/bin/python benchmarks/converter/create_delt_inputs.py \
      --dataset-name "$dataset" \
      --num-cores 11
  done
done
```

Prepare DELT-Hit inputs for the extended `2`-cycle large-library suite:

```bash
for dataset in \
  synthetic_2cycle_100bbpc_1m \
  synthetic_2cycle_100bbpc_10m \
  synthetic_2cycle_100bbpc_100m \
  synthetic_2cycle_100bbpc_1000m \
  synthetic_2cycle_1000bbpc_1m \
  synthetic_2cycle_1000bbpc_10m \
  synthetic_2cycle_1000bbpc_100m \
  synthetic_2cycle_1000bbpc_1000m
do
  ./.venv/bin/python benchmarks/converter/create_delt_inputs.py \
    --dataset-name "$dataset" \
    --num-cores 11
done
```

## 3. Run Split-Timing Benchmarks

Run [`benchmarks/run_split_timing.py`](./run_split_timing.py) against one prepared dataset.

The runner supports:

- `--tool deli`
- `--tool delt`
- `--tool both`

With `--data-dir "$DATA_DIR"`, it writes runtime logs and intermediate artifacts to:

- `$DATA_DIR/runtime/<dataset>/deli`
- `$DATA_DIR/runtime/<dataset>/delt`

It writes the canonical machine-readable report to:

- `benchmarks/tools/deli/<dataset>/timing.json`
- `benchmarks/tools/delt/<dataset>/timing.json`

### DELi Only

```bash
./.venv/bin/python benchmarks/run_split_timing.py \
  --dataset-name synthetic_4cycle_100m \
  --tool deli
```

Run DELi across all 12 datasets in-place:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
cd "$ROOT"

for cycles in 2 3 4; do
  for depth in 1m 10m 100m 1000m; do
    dataset="synthetic_${cycles}cycle_${depth}"
    ./.venv/bin/python benchmarks/run_split_timing.py \
      --dataset-name "$dataset" \
      --tool deli
  done
done
```

Run DELi across the extended `2`-cycle large-library suite:

```bash
for dataset in \
  synthetic_2cycle_100bbpc_1m \
  synthetic_2cycle_100bbpc_10m \
  synthetic_2cycle_100bbpc_100m \
  synthetic_2cycle_100bbpc_1000m \
  synthetic_2cycle_1000bbpc_1m \
  synthetic_2cycle_1000bbpc_10m \
  synthetic_2cycle_1000bbpc_100m \
  synthetic_2cycle_1000bbpc_1000m
do
  ./.venv/bin/python benchmarks/run_split_timing.py \
    --dataset-name "$dataset" \
    --tool deli
done
```

Submit one `18h` Slurm job per dataset for end-to-end DELi benchmarking on `$TMPDIR`:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
NUM_ERRORS=0
cd "$ROOT"

for cycles in 2 3 4; do
  for depth in 1m 10m 100m 1000m; do
    case "${cycles}:${depth}" in
      2:1m) reads_per_compound=10000 ;;
      2:10m) reads_per_compound=100000 ;;
      2:100m) reads_per_compound=1000000 ;;
      2:1000m) reads_per_compound=10000000 ;;
      3:1m) reads_per_compound=1000 ;;
      3:10m) reads_per_compound=10000 ;;
      3:100m) reads_per_compound=100000 ;;
      3:1000m) reads_per_compound=1000000 ;;
      4:1m) reads_per_compound=100 ;;
      4:10m) reads_per_compound=1000 ;;
      4:100m) reads_per_compound=10000 ;;
      4:1000m) reads_per_compound=100000 ;;
    esac
    dataset="synthetic_${cycles}cycle_${depth}"
    sbatch --time=18:00:00 --mem=64G --cpus-per-task=12 --job-name="bench_deli_${dataset}" --output="$HOME/logs/%j.out" --wrap "
cd $ROOT &&
DATA_DIR=\$TMPDIR/benchmarks/$dataset &&
./.venv/bin/python benchmarks/generate_synthetic_fastq.py \
  --num-cycles $cycles \
  --building-blocks-per-cycle 10 \
  --num-reads-per-compound $reads_per_compound \
  --num-errors $NUM_ERRORS \
  --output-dir \$DATA_DIR/data \
  --experiment-name $dataset &&
./.venv/bin/python benchmarks/converter/create_deli_inputs.py \
  --dataset-name $dataset \
  --data-dir \$DATA_DIR &&
./.venv/bin/python benchmarks/run_split_timing.py \
  --dataset-name $dataset \
  --data-dir \$DATA_DIR \
  --tool deli"
  done
done
```

Submit one Slurm job per dataset for the extended `2`-cycle large-library DELi suite:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
cd "$ROOT"

for bbpc in 100 1000; do
  for depth in 1m 10m 100m 1000m; do
    dataset="synthetic_2cycle_${bbpc}bbpc_${depth}"
    sbatch --time=18:00:00 --mem=64G --cpus-per-task=12 --job-name="bench_deli_${dataset}" --output="$HOME/logs/%j.out" --wrap "
cd $ROOT &&
DATA_DIR=\$TMPDIR/benchmarks/$dataset &&
./.venv/bin/python benchmarks/generate_synthetic_benchmark_matrix.py \
  --profile two_cycle_large_libraries \
  --num-errors 0 \
  --output-dir \$DATA_DIR/data &&
./.venv/bin/python benchmarks/converter/create_deli_inputs.py \
  --dataset-name $dataset \
  --data-dir \$DATA_DIR &&
./.venv/bin/python benchmarks/run_split_timing.py \
  --dataset-name $dataset \
  --data-dir \$DATA_DIR \
  --tool deli"
  done
done
```

### DELT-Hit Only

```bash
./.venv/bin/python benchmarks/run_split_timing.py \
  --dataset-name synthetic_4cycle_100m_err=1 \
  --tool delt
```

Run DELT-Hit across all 12 datasets in-place:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
NUM_ERRORS=0
cd "$ROOT"
for depth in 1m 10m 100m 1000m; do
  for cycles in 2 3 4; do
    dataset="synthetic_${cycles}cycle_${depth}"
    ./.venv/bin/python benchmarks/run_split_timing.py \
      --dataset-name "$dataset" \
      --tool delt
  done
done
```

Run DELT-Hit across the extended `2`-cycle large-library suite:

```bash
for dataset in \
  synthetic_2cycle_100bbpc_1m \
  synthetic_2cycle_100bbpc_10m \
  synthetic_2cycle_100bbpc_100m \
  synthetic_2cycle_100bbpc_1000m \
  synthetic_2cycle_1000bbpc_1m \
  synthetic_2cycle_1000bbpc_10m \
  synthetic_2cycle_1000bbpc_100m \
  synthetic_2cycle_1000bbpc_1000m
do
  ./.venv/bin/python benchmarks/run_split_timing.py \
    --dataset-name "$dataset" \
    --tool delt
done
```

Submit one `18h` Slurm job per dataset for end-to-end DELT-Hit benchmarking on `$TMPDIR`:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
cd "$ROOT"

for cycles in 2 3 4; do
  for depth in 1m 10m 100m 1000m; do
    case "${cycles}:${depth}" in
      2:1m) reads_per_compound=10000 ;;
      2:10m) reads_per_compound=100000 ;;
      2:100m) reads_per_compound=1000000 ;;
      2:1000m) reads_per_compound=10000000 ;;
      3:1m) reads_per_compound=1000 ;;
      3:10m) reads_per_compound=10000 ;;
      3:100m) reads_per_compound=100000 ;;
      3:1000m) reads_per_compound=1000000 ;;
      4:1m) reads_per_compound=100 ;;
      4:10m) reads_per_compound=1000 ;;
      4:100m) reads_per_compound=10000 ;;
      4:1000m) reads_per_compound=100000 ;;
    esac
    dataset="synthetic_${cycles}cycle_${depth}"
    sbatch --time=18:00:00 --mem=64G --cpus-per-task=12 --job-name="bench_delt_${dataset}" --output="$HOME/logs/%j.out" --wrap "
cd $ROOT &&
DATA_DIR=\$TMPDIR/benchmarks/$dataset &&
./.venv/bin/python benchmarks/generate_synthetic_fastq.py \
  --num-cycles $cycles \
  --building-blocks-per-cycle 10 \
  --num-reads-per-compound $reads_per_compound \
  --num-errors $NUM_ERRORS \
  --output-dir \$DATA_DIR/data \
  --experiment-name $dataset &&
./.venv/bin/python benchmarks/converter/create_delt_inputs.py \
  --dataset-name $dataset \
  --data-dir \$DATA_DIR \
  --num-cores 11 &&
./.venv/bin/python benchmarks/run_split_timing.py \
  --dataset-name $dataset \
  --data-dir \$DATA_DIR \
  --tool delt"
  done
done
```

Submit one Slurm job per dataset for the extended `2`-cycle large-library DELT-Hit suite:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
cd "$ROOT"

for bbpc in 100 1000; do
  for depth in 1m 10m 100m 1000m; do
    dataset="synthetic_2cycle_${bbpc}bbpc_${depth}"
    sbatch --time=18:00:00 --mem=64G --cpus-per-task=12 --job-name="bench_delt_${dataset}" --output="$HOME/logs/%j.out" --wrap "
cd $ROOT &&
DATA_DIR=\$TMPDIR/benchmarks/$dataset &&
./.venv/bin/python benchmarks/generate_synthetic_benchmark_matrix.py \
  --profile two_cycle_large_libraries \
  --num-errors 0 \
  --output-dir \$DATA_DIR/data &&
./.venv/bin/python benchmarks/converter/create_delt_inputs.py \
  --dataset-name $dataset \
  --data-dir \$DATA_DIR \
  --num-cores 11 &&
./.venv/bin/python benchmarks/run_split_timing.py \
  --dataset-name $dataset \
  --data-dir \$DATA_DIR \
  --tool delt"
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

## Manual Command Sequence

If you want to run the prepared inputs manually instead of using `run_split_timing.py`, the equivalent commands are:

```bash
#!/usr/bin/env bash
set -euo pipefail

ROOT="/users/amarti51/projects/delt-hit"
DATASET="synthetic_4cycle_100m"

DELI_DIR="$ROOT/benchmarks/tools/deli/$DATASET"
DELT_DIR="$ROOT/benchmarks/tools/delt/$DATASET"
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
