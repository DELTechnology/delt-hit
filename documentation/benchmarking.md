# Synthetic Benchmarking

The synthetic benchmark workflow is described in [`benchmarks/demultiplex/README.md`](../benchmarks/demultiplex/README.md).

Use the demultiplex runbook for:

- dataset generation with [`benchmarks/demultiplex/generate_synthetic_fastq.py`](../benchmarks/demultiplex/generate_synthetic_fastq.py)
- DELi and DELT-Hit input conversion
- scratch-space Slurm execution with `--data-dir "$TMPDIR/..."`
- timing report collection under `benchmarks/demultiplex/tools/<tool>/<dataset>/timing.json`
