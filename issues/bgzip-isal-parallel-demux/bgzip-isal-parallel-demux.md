# Faster Demultiplexing: bgzip + isal + Multiprocessing

## Problem

`get_counts()` in `postprocess.py` iterates line-by-line over a gzipped FASTQ header file
(`reads_with_adapters.gz`) using Python's built-in `gzip` module. Two things are slow:

1. **Decompression** — stdlib `gzip` is single-threaded and slow.
2. **Parsing loop** — CPython has high per-line overhead; `extract_ids` runs once per read.

For a typical DEL experiment with tens of millions of reads this can take many minutes.

---

## What Was Done (branch `feat/bgzip-isal-parallel-demux`)

### Files changed

| File | Change |
|------|--------|
| `src/delt_hit/demultiplex/preprocess.py` | Final `zgrep` pipe now writes through `bgzip -@ {num_cores}` → `reads_with_adapters.bgz` |
| `src/delt_hit/demultiplex/postprocess.py` | `get_counts` replaced with isal + multiprocessing batch variant |
| `src/delt_hit/cli/demultiplex/api.py` | `input_path` updated to `.bgz` |
| `pyproject.toml` | `isal` added to dependencies |

### How it works now

```
bgzip -@ N  (parallel compression, written by the bash pipeline)
     ↓
isal.igzip.open  (fast single-threaded decompression stream)
     ↓
_batches()  (yields lists of 50k lines)
     ↓
multiprocessing.Pool.imap(_count_batch, ...)  (N workers parse in parallel)
     ↓
_merge()  (main process accumulates partial dicts)
```

**What this buys:**
- `isal` decompresses 3–5× faster than stdlib `gzip` (Intel ISA-L SIMD routines).
- Parsing (`extract_ids` + dict ops) runs across all CPU cores.
- `bgzip -@` compresses the output file faster (parallel threads).
- Memory is bounded: `imap` keeps only one batch in-flight per worker.

**What this does NOT do:**
- Decompression is still serial — one `isal` stream feeds all workers.
- If isal decompression is fast enough to saturate worker processes, this is fine.
  If decompression is the bottleneck, workers will sit idle waiting for batches.

### Prerequisite

`bgzip` must be available on the system (ships with `htslib`: `brew install htslib`).

---

## Profiling Decision Point

Before going further, profile whether decompression or parsing is the bottleneck:

```python
# rough timing split: time decompression alone vs. decompression + parsing
from isal import igzip
import time

with igzip.open('reads_with_adapters.bgz', 'rt') as f:
    t0 = time.perf_counter()
    for line in f:
        pass  # decompression only
    print("decompress only:", time.perf_counter() - t0)
```

Compare this against the full `get_counts` wall time with `num_workers=1`.

- **If decompress-only ≈ single-worker full time**: decompression is the bottleneck.
  Parallel decompression (see below) would help.
- **If decompress-only << single-worker full time**: parsing is the bottleneck.
  The current multiprocessing approach already addresses this; no need to change.

In practice, `isal` throughput on a modern CPU is ~1–3 GB/s decompressed.
`extract_ids` does a string split + two list comprehensions per line.
For short FASTQ headers (~100 bytes each), parsing ~10M reads/s is typical in CPython.
At 100 bytes/line that is 1 GB/s — roughly matching isal throughput, so decompression
may or may not be the limiter depending on the machine and read count.

---

## Full Parallel Decompression + Parsing (not yet implemented)

To decompress and parse in parallel, each worker must open the bgzf file independently
and seek to its own block range. This requires:

### 1. Get block offsets

BGZF files consist of independent compressed blocks (≤ 64 KB each). Block virtual offsets
can be listed with:

```bash
bgzip --list reads_with_adapters.bgz
```

Or in Python using `pysam`:

```python
import pysam
offsets = list(pysam.libcbgzf.get_block_offsets('reads_with_adapters.bgz'))
```

### 2. Split blocks across workers

Divide the list of `(start_offset, end_offset)` pairs into N roughly equal chunks.

### 3. Each worker seeks and decompresses independently

```python
import pysam

def _count_block_range(args):
    path, start_voffset, end_voffset = args
    with pysam.libcbgzf.BGZFile(path, 'rb') as f:
        f.seek(start_voffset)
        counts = {}
        while f.tell() < end_voffset:
            line = f.readline().decode()
            # ... extract_ids and count
        return counts
```

### Required changes

| File | Change |
|------|--------|
| `postprocess.py` | Replace `_batches` + `isal.igzip.open` with block-offset splitting + `pysam.libcbgzf.BGZFile` |
| `pyproject.toml` | Add `pysam` dependency (heavy: pulls in htslib C bindings) |

### Expected additional speedup

| Stage | Current | With block-parallel |
|-------|---------|---------------------|
| Decompression | ~1 thread (isal) | N threads (one per worker) |
| Parsing | N threads | N threads (same) |
| End-to-end | bottlenecked by whichever is slower | both scale with N |

On an 8-core machine decompression parallelism alone could give another 2–4× on top of
the current implementation, but **only if decompression is the measured bottleneck**.
If isal already feeds workers fast enough, the added complexity of `pysam` + virtual
offsets is not worth it.

### Recommendation

Profile first. If the `decompress-only` timing is within 2× of the full `get_counts`
time, then block-parallel decompression is worth implementing. Otherwise, the current
`isal` + `multiprocessing` approach is sufficient.

### Installation on Mac

```bash
brew install htslib
#brew install llvm
```

---

## Rust Extension (branch `feat/bgzip-isal-parallel-demux`)

A Rust PyO3 extension (`delt_hit._ext`) was added as an alternative to the Python
multiprocessing path. It uses `noodles-bgzf` for bgzf decoding and `rayon` for
parallel counting.

### Files changed

| File | Change |
|------|--------|
| `src/rust/lib.rs` | Rust implementation of `get_counts` |
| `Cargo.toml` | New crate: `delt-hit-ext`, lib name `_ext`, deps: pyo3 / noodles-bgzf / rayon / memchr |
| `pyproject.toml` | Build backend switched from `setuptools` to `maturin`; added `[tool.maturin]` with `python-source = "src"`, `module-name = "delt_hit._ext"` |
| `src/delt_hit/demultiplex/postprocess.py` | `get_counts` dispatches to Rust when `delt_hit._ext` is importable, falls back to Python path otherwise |
| `scripts/profile_get_counts.py` | Added `bench_get_counts_rust` step and Rust row in summary |

### How the Rust path works

```
bgzf::Reader (noodles-bgzf, serial decompression)
     ↓
Vec<Vec<u8>>  (all lines collected into memory)
     ↓
rayon::par_chunks  (parallel parsing across all CPU cores)
     ↓
Mutex<HashMap>  (merge partial counts under lock)
     ↓
PyDict  (returned to Python)
```

### Profiling results (5 M synthetic reads, Apple M-series, 11 cores)

| Step | Time | Throughput |
|------|------|------------|
| Decompress only (isal) | 2.28 s | 2.20 M lines/s |
| extract_ids only (Python, 1 thread) | 8.86 s | 0.56 M lines/s |
| get_counts Python w=1 | 9.53 s | 0.52 M reads/s |
| get_counts Python w=11 | 5.01 s | 1.00 M reads/s |
| **get_counts Rust** | **8.83 s** | **0.57 M reads/s** |

Decompression accounts for 24% of single-worker time → **parsing is the bottleneck**.

### Why the Rust path is not faster than Python w=11

The Rust implementation uses a two-phase approach:

1. **Read phase (serial):** all 5 M lines are decompressed and collected into a
   `Vec<Vec<u8>>` — this allocates 5 M heap buffers and takes roughly the same time
   as isal-only decompression plus allocation overhead.
2. **Parse phase (rayon parallel):** `par_chunks` distributes parsing across cores.

The dominant cost is the ~5–6 s of parsing work spread across all cores, which rayon
handles. However, the heavy upfront allocation in phase 1 adds latency that the Python
`multiprocessing` approach avoids by streaming batches directly into worker processes.
The net result is parity with Python w=1 rather than a speedup over Python w=11.

### Conclusion / next steps

The **Python w=11 path (isal + multiprocessing)** remains the fastest at 5.01 s.
The Rust extension does not improve on this for the current workload.

Options if further speedup is needed:

| Option | Expected gain | Complexity |
|--------|--------------|------------|
| Rust: avoid `Vec<Vec<u8>>` — parse lines directly during streaming read | reduces allocation, may close gap with Python w=11 | medium |
| Block-parallel bgzf (pysam or noodles seek) | decompression scales with N, worth ~1.5–2× extra | high (see section above) |
| Keep current Python w=11 path | no change needed | none |

Given that decompression is only 24% of total time, block-parallel decompression
alone is unlikely to move the needle much. Streaming parse in Rust (avoiding the
collect step) is the more promising Rust-side improvement.