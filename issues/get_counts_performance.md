# get_counts Performance: Parallel Scaling and C-level Parsing

## Problem

`get_counts` (postprocess.py:99) is the wall-clock bottleneck of the demultiplexing pipeline. Production experiments produce files with **billions of reads** that need to be parsed and counted. The original serial implementation took ~20s on 5M reads; at 1B reads that is ~67 minutes single-threaded.

Two bottlenecks compound:

1. **I/O / decompression** — `gzip.open('rt')` uses stdlib zlib (slow) and pays a UTF-8 decode per line.
2. **Per-line parsing** — `extract_ids` (postprocess.py:10) allocates multiple Python objects per line: `strip()`, `split('?')`, `split('.')` per token, `filter()` with lambdas, `int()` conversions.

## Current State (after first round of optimisation)

`get_counts` now accepts `num_workers` and `chunk_size_bytes`. When `num_workers > 1`:

- `_iter_byte_chunks` (postprocess.py:46) decompresses via **isal.igzip** (Intel ISA-L, ~3× faster than zlib) and yields aligned byte chunks.
- `_count_bytes_chunk` (postprocess.py:26) parses chunks in pure Python bytes operations (no UTF-8 decode, no defaultdict).
- Chunks are dispatched to a `multiprocessing.Pool` via `imap_unordered`.

When `num_workers == 1` the old serial path (stdlib gzip + string-mode `extract_ids`) is still used — it does **not** benefit from isal.

### Benchmark results (48-core node, 5M synthetic reads, 78 MB compressed)

| Mode | Time | Speedup |
|---|---|---|
| serial-python (stdlib gzip + `extract_ids`) | 26.5 s | 1.0× |
| parallel-python (isal + Python `_count_bytes_chunk`) | 3.7 s | 7.1× |
| parallel-cython (isal + Cython `parse_chunk`) | 4.3 s | 6.2× |

Output verified identical across all 100 selections × 5M reads.

Benchmark script: `tests/benchmark_get_counts.py`

## Key Finding: Bottleneck is Decompression, Not Parsing

The Cython extension (`parse_chunk`) is **not faster** than the Python parallel version — it is marginally slower. This reveals that the bottleneck has shifted:

With 48 workers consuming chunks, the main process cannot supply them fast enough. All workers spend most of their time idle waiting for the next chunk from `_iter_byte_chunks`. Making each worker's parse step faster with C-level code just increases idle wait time — it does not reduce wall time.

**The real bottleneck is serial gzip decompression in the main process** (isal is already ~3× faster than stdlib zlib, but decompression is still single-threaded).

## Remaining Bottleneck

To go faster the architecture must change so that decompression can also be parallelised. Options:

## Potential Solutions

### Option A — Cython `parse_chunk` (implemented, no net gain)

Write `src/delt_hit/demultiplex/_parse.pyx` with a single `parse_chunk(bytes) -> dict` function. Uses `libc.string.memchr` for token boundary detection and manual digit parsing — no Python objects allocated during tokenisation. Drop-in replacement for `_count_bytes_chunk`.

**Expected gain:** 3–5× speedup in worker CPU time (parse dominates). Combined with existing parallelism: estimated **15–20× over original serial baseline**.

**Cost:** adds Cython build step (`setup.py`, `Cython>=3.0` in build-system). Python fallback kept for environments without a compiler.

### Option B — Re-encode IDs as fixed-width integers in the sequence tag

Instead of `?S0.2?B1.5?B2.3?B3.7?S1.1`, encode IDs as fixed-width zero-padded tokens so parsing is a direct offset read rather than a search. Requires changing the upstream tag encoding in `preprocess.py` / cutadapt adapter FASTA. Breaks compatibility with existing files.

### Option C — Pre-split input, run cutadapt on each part independently

Split the raw input FASTQ into N parts before cutadapt, producing N separate `reads_with_adapters_N.gz` files. Count each in parallel without multiprocessing within a single file. More I/O, harder orchestration. Avoids any API complexity.

### Option D — Polars / Arrow for counting

Load all header lines into a Polars DataFrame, split columns with native string expressions, group-by and count. Avoids Python per-line overhead. Requires buffering the full file or streaming batches. May be memory-intensive at 1B reads.

## Recommended Next Steps (open tasks)

See `issues/get_counts_performance_plan.md`.
