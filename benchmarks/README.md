# Benchmarks

## Prerequisites

- **DELT-Hit:** the project Pixi environment (`pixi install` from the repo root).
- **DELi:** clone into `other_tools/DELi`, check out the benchmarked commit, install from source into its own venv, and append to PATH:

  ```bash
  git clone https://github.com/Popov-Lab-UNC/DELi.git other_tools/DELi
  git -C other_tools/DELi checkout 8b8f8b570917a60ba76a050f2d1428669eec95c5  # assumed commit of the stored results (tag 0.1.3)
  python -m venv other_tools/DELi/.venv
  other_tools/DELi/.venv/bin/pip install -e other_tools/DELi
  export PATH="$PATH:$PWD/other_tools/DELi/.venv/bin"
  ```

  The DELi venv is kept separate from the Pixi environment. Appending to PATH ensures Pixi retains priority.

---

Benchmark documentation now lives in benchmark-specific subfolders:

- [`benchmarks/demultiplex/README.md`](./demultiplex/README.md) for the synthetic FASTQ demultiplexing benchmarks
