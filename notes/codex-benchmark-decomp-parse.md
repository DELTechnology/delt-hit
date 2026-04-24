## Benchmark follow-up notes

- Use `main` as the reference for the current demultiplex implementation rather than `benchmarking/synthetic-deli-setup`.
- Confirm the current `main` implementation uses stdlib `gzip` for the serial `get_counts(..., num_workers=1)` path and `isal.igzip` for the parallel chunk reader in `_iter_byte_chunks`.
- Extend the existing `tests/benchmark_get_counts.py` benchmark so it reports decompression-only and parsing-only timings in addition to the end-to-end serial and parallel `get_counts` timings.
