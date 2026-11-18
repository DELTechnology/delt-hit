# Benchmarks

## Prerequisites

- **DELT-Hit:** the project Pixi environment (`pixi install` from the repo root).
- **DELi:** clone into `other_tools/DELi`, check out the benchmarked commit, install from source into its own venv, and append to PATH:

  ```bash
  git clone https://github.com/Popov-Lab-UNC/DELi.git other_tools/DELi
  git -C other_tools/DELi checkout 8b8f8b570917a60ba76a050f2d1428669eec95c5  # tag 0.1.3
  python -m venv other_tools/DELi/.venv
  other_tools/DELi/.venv/bin/pip install -e other_tools/DELi
  export PATH="$PATH:$PWD/other_tools/DELi/.venv/bin"
  ```

  The DELi venv is kept separate from the Pixi environment. Appending to PATH ensures Pixi retains priority.

  The stored DELi results under `demultiplex/tools/deli/` were produced on 2026-05-01. The DELi version was not recorded at run time. Commit `8b8f8b5` (tag `0.1.3`, 2026-04-22) was the head of DELi `main` from before the benchmark setup (2026-04-23) until 2026-07-07, so it is the version those runs used.

---

Benchmark documentation now lives in benchmark-specific subfolders:

- [`benchmarks/demultiplex/README.md`](./demultiplex/README.md) for the synthetic FASTQ demultiplexing benchmarks
