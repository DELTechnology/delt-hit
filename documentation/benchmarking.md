# Synthetic Benchmarking

This branch now includes a reproducible synthetic benchmark workflow for DELi and DELT-Hit.

The benchmark runner does three things:

1. Regenerates synthetic DELi fixtures for 2-, 3-, and 4-cycle decoding.
2. Regenerates the synthetic 2-cycle DELT-Hit workbook and FASTQ fixture in an isolated output directory.
3. Runs both toolchains and writes a machine-readable summary with timings and count totals.

## Prerequisites

The benchmark assumes both local virtual environments already exist:

- DELT-Hit: [/.venv](/Users/adrianomartinelli/projects/delt-hit/.venv)
- DELi: [/other_tools/DELi/.venv](/Users/adrianomartinelli/projects/delt-hit/other_tools/DELi/.venv)

## Run Everything

```bash
./.venv/bin/python scripts/run_synthetic_benchmarks.py
```

Outputs are written under [`target/benchmarks/synthetic`](/Users/adrianomartinelli/projects/delt-hit/target/benchmarks/synthetic), including:

- regenerated DELi fixtures
- regenerated DELT-Hit synthetic inputs
- per-step command logs
- `benchmark_summary.json`

## Useful Variants

Only DELi:

```bash
./.venv/bin/python scripts/run_synthetic_benchmarks.py --tools deli
```

Only DELT-Hit:

```bash
./.venv/bin/python scripts/run_synthetic_benchmarks.py --tools delt-hit
```

Increase depth per compound:

```bash
./.venv/bin/python scripts/run_synthetic_benchmarks.py --num-reads-per-compound 10
```

Reuse the existing benchmark directory:

```bash
./.venv/bin/python scripts/run_synthetic_benchmarks.py --keep-existing
```

## What Gets Compared

For DELi, the runner measures:

- `decode run`
- `decode collect`
- `decode count`

For DELT-Hit, the runner measures:

- `init`
- `demultiplex prepare`
- execution of the generated `demultiplex.sh`
- `demultiplex process`

For the shared 2-cycle dataset, the runner also normalizes DELi `bb_ids` into integer codes and checks that the final count table matches DELT-Hit’s `code_1`, `code_2`, `count` output exactly.

## Current Caveats

- DELT-Hit is benchmarked only on the 2-cycle synthetic fixture in this branch.
- The DELT-Hit CLI `demultiplex run` path does not finish the count table on its own here, so the benchmark explicitly runs `prepare`, the generated shell script, and `process`.
- DELi already ships a `decode summarize` command, but this benchmark computes its own summary from the counted parquet and decode statistics so the output stays consistent across runs.
