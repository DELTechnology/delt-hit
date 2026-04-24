# `get_counts()` Optimization Notes

This note documents the changes made to `get_counts()` and the benchmark results collected afterward.

## What changed

### 1. Replace split-heavy parsing in the hot path

In [postprocess.py](/Users/adrianomartinelli/projects/delt-hit/src/delt_hit/demultiplex/postprocess.py), the counting path now uses a bytes parser with a small state machine instead of repeatedly doing:

- `split("?")`
- `split(".")`
- string filtering for `S` and `B`

The new hot-path parser is `_parse_labeled_line_bytes(...)`.

### 2. Replace nested partial maps with a flat counter key

Worker-local partial counts used to look like:

```python
{
  selection_ids: {
    barcodes: count
  }
}
```

Now workers count into a flat `Counter` keyed by a single tuple:

```python
(selection_0, selection_1, ..., -1, code_0, code_1, ...)
```

This makes worker-side counting and parent-side merging cheaper because the parent can use `Counter.update(...)` instead of nested Python loops.

### 3. Make chunk size configurable

`get_counts(...)` already accepted `chunk_size_bytes`, but DELT-Hit processing had it effectively fixed in practice. The processing path now reads:

- `experiment.chunk_size_bytes`

from the config and passes it through to `get_counts(...)`.

Relevant code:

- [postprocess.py](/Users/adrianomartinelli/projects/delt-hit/src/delt_hit/demultiplex/postprocess.py)
- [api.py](/Users/adrianomartinelli/projects/delt-hit/src/delt_hit/cli/demultiplex/api.py)

## Benchmarks

### Baseline before the change

The pre-change baseline below was supplied from the existing benchmark run on `5,000,000` reads:

```text
Generating 5,000,000 synthetic reads …
  done in 11.1s  (59 MB compressed)

Decompression only (isal) …
  1.35s  (3.71 M lines/s)

Materializing plain-text benchmark input …
  wrote 5,000,000 lines

Parsing only (extract_ids on preloaded lines) …
  7.74s  (0.65 M lines/s)

Serial (num_workers=1) …
  10.90s

Parallel (num_workers=11) …
  3.25s
  speedup: 3.4x
```

### Post-change benchmark on `5,000,000` reads

Command:

```bash
./.venv/bin/python tests/benchmark_get_counts.py --n-reads 5000000 --chunk-size-bytes 5000000
```

Results:

```text
Decompression only (isal): 1.43s
Parsing only (extract_ids on preloaded lines): 15.39s
Serial get_counts: 21.37s
Parallel get_counts (11 workers): 6.11s
Speedup: 3.5x
```

### Important interpretation

This result is mixed:

- `get_counts()` got slower after the refactor
- scaling stayed roughly similar (`3.4x` to `3.5x`)
- the standalone `extract_ids(...)` benchmark also got slower

That last point is expected because `extract_ids(...)` now goes through `line.encode()` to reuse the bytes parser. That makes the standalone parse-only benchmark harsher than before.

However, the bigger finding is that the full `get_counts()` path also regressed materially, so the new bytes/state-machine parser is not yet a win in Python.

## 25M-read chunk-size sweep

Commands were run with:

```bash
./.venv/bin/python tests/benchmark_get_counts.py --n-reads 25000000 --chunk-size-bytes <SIZE>
```

Three chunk sizes were tested:

- `5,000,000`
- `25,000,000`
- `100,000,000`

Saved raw outputs:

- [benchmark_5000000.txt](/Users/adrianomartinelli/projects/delt-hit/target/benchmarks/get_counts/benchmark_5000000.txt)
- [benchmark_25000000.txt](/Users/adrianomartinelli/projects/delt-hit/target/benchmarks/get_counts/benchmark_25000000.txt)
- [benchmark_100000000.txt](/Users/adrianomartinelli/projects/delt-hit/target/benchmarks/get_counts/benchmark_100000000.txt)

### Results

| Chunk size | Decompress only | Parse only | Serial `get_counts` | Parallel `get_counts` | Speedup |
|---|---:|---:|---:|---:|---:|
| `5,000,000` | `6.94s` | `77.40s` | `109.99s` | `32.62s` | `3.4x` |
| `25,000,000` | `7.10s` | `75.78s` | `106.42s` | `33.08s` | `3.2x` |
| `100,000,000` | `7.11s` | `77.87s` | `109.28s` | `41.17s` | `2.7x` |

## Findings

### 1. The new parser did not improve performance

The requested changes were implemented, but the benchmark results show the new version is slower than the old one.

That means:

- the bytes/state-machine parser is not currently beating the old split-based parser in pure Python
- flattening worker partials into a `Counter` did not offset the parser regression

### 2. Small-to-medium chunks are better than very large chunks

For the current implementation:

- `5MB` was best
- `25MB` was slightly worse
- `100MB` was clearly worse

So in this version, larger chunks do not improve scaling. They likely increase tail latency and reduce load balancing enough to outweigh lower scheduling overhead.

### 3. The main bottleneck is still parsing / Python work

At `25M` reads:

- decompression is around `7s`
- parse-only is around `76-78s`
- serial total is around `106-110s`

So the dominant cost remains parser/object overhead, not decompression.

### 4. Parallel scaling is still limited

With `11` workers, the best result here was about `3.4x`, not near-linear scaling.

The current evidence suggests the remaining limits are still mostly:

- Python parser cost inside each worker
- process coordination / IPC
- parent-side merge and reconstruction back into nested dicts

## What I would do next

The current refactor is informative, but it is not the optimization we want to keep.

The next promising options are:

1. Revert the hot-path parser change and benchmark a flatter merge strategy independently.
2. Keep the original parser, but profile these phases separately inside `get_counts()`:
   - chunk read/decompress
   - worker parse/count
   - parent merge
3. Benchmark a more compact flat key representation without reconstructing nested tuples in the hot path.
4. If we want a parser rewrite, move it to a compiled path:
   - Cython
   - Rust
   - Numba on a more array-like representation

## Bottom line

The requested optimization ideas were implemented and benchmarked.

The measurements say:

- the new bytes/state-machine parser is not faster in Python
- flat worker counters alone are not enough to recover that cost
- `5MB` chunks remain the best of the three tested sizes
- the real bottleneck is still Python parsing work, not decompression

So this was a useful experiment, but not yet a winning optimization.
