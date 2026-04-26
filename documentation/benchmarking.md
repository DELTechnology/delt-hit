# Synthetic Benchmarking

The synthetic benchmark workflow is now maintained in [`benchmarks/README.md`](/Users/adrianomartinelli/projects/delt-hit/benchmarks/README.md).

Use that runbook for:

- dataset generation with [`benchmarks/generate_synthetic_fastq.py`](/Users/adrianomartinelli/projects/delt-hit/benchmarks/generate_synthetic_fastq.py)
- DELi and DELT-Hit input conversion
- scratch-space Slurm execution with `--data-dir "$TMPDIR/..."`
- timing report collection under `benchmarks/tools/<tool>/<dataset>/timing.json`
