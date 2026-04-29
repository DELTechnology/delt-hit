# Synthetic Benchmarking

The synthetic benchmark workflows are now split across:

- [`benchmarks/demultiplex/README.md`](/Users/adrianomartinelli/projects/delt-hit/benchmarks/demultiplex/README.md)
- [`benchmarks/enrichment/README.md`](/Users/adrianomartinelli/projects/delt-hit/benchmarks/enrichment/README.md)

Use the demultiplex runbook for:

- dataset generation with [`benchmarks/demultiplex/generate_synthetic_fastq.py`](/Users/adrianomartinelli/projects/delt-hit/benchmarks/demultiplex/generate_synthetic_fastq.py)
- DELi and DELT-Hit input conversion
- scratch-space Slurm execution with `--data-dir "$TMPDIR/..."`
- timing report collection under `benchmarks/demultiplex/tools/<tool>/<dataset>/timing.json`
