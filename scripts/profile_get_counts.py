"""Profile decompression vs. parsing bottleneck for get_counts().

Usage:
  uv run python scripts/profile_get_counts.py --num-reads 5000000

Steps:
1. Generate a synthetic reads_with_adapters.bgz file.
2. Time isal decompression only (no parsing).
3. Time extract_ids only (parsing, single-threaded, no pool overhead).
4. Time get_counts() with 1 worker (decompression + parsing, serial).
5. Time get_counts() with N workers (decompression + parallel parsing).

Line format (cutadapt info, '?'-delimited):
  <read_id>?...?S1.2?...?B1.4?B2.7?...
  - tokens containing 'S' carry selection IDs (last dot-separated field)
  - tokens containing 'B' carry barcode IDs (last dot-separated field, 1-indexed)
"""

import argparse
import random
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from multiprocessing import cpu_count
from pathlib import Path

# ---------------------------------------------------------------------------
# Synthetic data generation
# ---------------------------------------------------------------------------

def _make_line(rng: random.Random, num_selections: int, num_barcodes: int,
               sel_range: int, bc_range: int) -> str:
    sel_ids = [rng.randint(0, sel_range - 1) for _ in range(num_selections)]
    bc_ids  = [rng.randint(0, bc_range - 1)  for _ in range(num_barcodes)]
    sel_tokens = [f"X.S{i}.{s}" for i, s in enumerate(sel_ids, 1)]
    bc_tokens  = [f"X.B{i}.{b}" for i, b in enumerate(bc_ids, 1)]
    all_tokens = sel_tokens + bc_tokens
    rng.shuffle(all_tokens)
    return "read_" + str(rng.randint(0, 10**9)) + "?" + "?".join(all_tokens) + "\n"


def generate_bgz(path: Path, num_reads: int, seed: int = 42,
                 num_selections: int = 2, num_barcodes: int = 3,
                 sel_range: int = 4, bc_range: int = 100) -> None:
    rng = random.Random(seed)
    tmp_txt = path.with_suffix(".tmp.txt")
    print(f"  Writing {num_reads:,} synthetic lines to {tmp_txt} …", flush=True)
    with open(tmp_txt, "w") as fh:
        for _ in range(num_reads):
            fh.write(_make_line(rng, num_selections, num_barcodes, sel_range, bc_range))
    print(f"  Compressing with bgzip → {path} …", flush=True)
    subprocess.run(["bgzip", "-f", "-@", str(cpu_count()), str(tmp_txt)],
                   check=True)
    Path(str(tmp_txt) + ".gz").rename(path)  # bgzip appends .gz → rename to .bgz
    print(f"  Done. Size: {path.stat().st_size / 1e6:.1f} MB", flush=True)


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------

def bench_extract_ids_only(path: Path) -> tuple[float, int]:
    from isal import igzip
    from delt_hit.demultiplex.postprocess import extract_ids
    t0 = time.perf_counter()
    n = 0
    with igzip.open(path, "rt") as f:
        for line in f:
            extract_ids(line)
            n += 1
    return time.perf_counter() - t0, n


def bench_decompress_only(path: Path) -> float:
    from isal import igzip
    t0 = time.perf_counter()
    with igzip.open(path, "rt") as f:
        for _ in f:
            pass
    return time.perf_counter() - t0


def bench_plain_read_only(path: Path) -> float:
    t0 = time.perf_counter()
    with open(path) as f:
        for _ in f:
            pass
    return time.perf_counter() - t0


def bench_plain_extract_ids_only(path: Path) -> tuple[float, int]:
    from delt_hit.demultiplex.postprocess import extract_ids
    t0 = time.perf_counter()
    n = 0
    with open(path) as f:
        for line in f:
            extract_ids(line)
            n += 1
    return time.perf_counter() - t0, n


def bench_plain_get_counts(path: Path, num_workers: int,
                            batch_size: int = 150_000) -> float:
    import multiprocessing
    from delt_hit.demultiplex.postprocess import _merge

    counts: defaultdict = defaultdict(dict)

    def _batches(f):
        batch = []
        for line in f:
            batch.append(line)
            if len(batch) == batch_size:
                yield batch
                batch = []
        if batch:
            yield batch

    t0 = time.perf_counter()
    with open(path) as f, multiprocessing.Pool(num_workers) as pool:
        for partial in pool.imap(_count_batch_worker, _batches(f)):
            _merge(counts, partial)
    return time.perf_counter() - t0


def bench_plain_get_counts_offset(path: Path, num_workers: int,
                                   batch_size: int = 150_000) -> float:
    import multiprocessing, os
    from delt_hit.demultiplex.postprocess import _count_file_range, _merge

    size = os.path.getsize(path)
    chunk = size // num_workers
    ranges = [
        (path, i * chunk, (i + 1) * chunk if i < num_workers - 1 else size)
        for i in range(num_workers)
    ]

    counts: defaultdict = defaultdict(dict)
    t0 = time.perf_counter()
    with multiprocessing.Pool(num_workers) as pool:
        for partial in pool.map(_count_file_range, ranges):
            _merge(counts, partial)
    return time.perf_counter() - t0


def bench_decompress_then_offset(bgz_path: Path, num_workers: int) -> float:
    """New get_counts approach: decompress bgz to temp file, then seek-parallel."""
    import multiprocessing, os, tempfile
    from isal import igzip
    from delt_hit.demultiplex.postprocess import _count_file_range, _merge

    counts: defaultdict = defaultdict(dict)
    t0 = time.perf_counter()

    with tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp:
        tmp_path = tmp.name
    try:
        with igzip.open(bgz_path, 'rb') as src, open(tmp_path, 'wb') as dst:
            while chunk := src.read(1 << 20):
                dst.write(chunk)

        size = os.path.getsize(tmp_path)
        chunk = size // num_workers
        ranges = [(tmp_path, i * chunk, (i + 1) * chunk if i < num_workers - 1 else size)
                  for i in range(num_workers)]

        with multiprocessing.Pool(num_workers) as pool:
            for partial in pool.map(_count_file_range, ranges):
                _merge(counts, partial)
    finally:
        os.unlink(tmp_path)

    return time.perf_counter() - t0


def bench_get_counts_rust(path: Path) -> tuple[float, int]:
    from delt_hit._ext import get_counts as _rust
    t0 = time.perf_counter()
    counts = _rust(path)
    total = sum(sum(v.values()) for v in counts.values())
    return time.perf_counter() - t0, total


def _count_batch_worker(lines):
    from delt_hit.demultiplex.postprocess import extract_ids
    counts: defaultdict = defaultdict(dict)
    for line in lines:
        sel, bc = extract_ids(line)
        inner = counts[sel]
        inner[bc] = inner.get(bc, 0) + 1
    return dict(counts)


def bench_get_counts(path: Path, num_reads: int, num_workers: int,
                     batch_size: int = 150_000) -> float:
    from isal import igzip
    import multiprocessing
    from delt_hit.demultiplex.postprocess import _merge

    counts: defaultdict = defaultdict(dict)

    def _batches(f):
        batch = []
        for line in f:
            batch.append(line)
            if len(batch) == batch_size:
                yield batch
                batch = []
        if batch:
            yield batch

    t0 = time.perf_counter()
    with igzip.open(path, "rt") as f, \
            multiprocessing.Pool(num_workers) as pool:
        for partial in pool.imap(_count_batch_worker, _batches(f)):
            _merge(counts, partial)
    return time.perf_counter() - t0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--num-reads", type=int, default=5_000_000,
                        help="Number of synthetic reads (default: 5M)")
    parser.add_argument("--bgz-path", type=Path, default=None,
                        help="Re-use an existing .bgz file instead of generating one")
    parser.add_argument("--workers", type=int, nargs="+", default=None,
                        help="Worker counts to benchmark (default: 1 and cpu_count())")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    workers = args.workers or ([1, cpu_count()] if cpu_count() > 1 else [1])

    # --- generate or reuse synthetic file ---
    if args.bgz_path is not None:
        bgz_path = args.bgz_path
        if not bgz_path.exists():
            sys.exit(f"File not found: {bgz_path}")
        print(f"Re-using existing file: {bgz_path}")
        num_reads = args.num_reads
    else:
        tmp = tempfile.NamedTemporaryFile(suffix=".bgz", delete=False)
        bgz_path = Path(tmp.name)
        tmp.close()
        print(f"\n[1/4] Generating synthetic dataset ({args.num_reads:,} reads) …")
        generate_bgz(bgz_path, args.num_reads, seed=args.seed)
        num_reads = args.num_reads
        try:
            _run_benchmarks(bgz_path, num_reads, workers)
        finally:
            bgz_path.unlink(missing_ok=True)
        return

    _run_benchmarks(bgz_path, num_reads, workers)


def _run_benchmarks(bgz_path: Path, num_reads: int, workers: list[int]) -> None:
    print(f"\n[2/x] Benchmarking decompression only (isal, no parsing) …")
    decomp_time = bench_decompress_only(bgz_path)
    print(f"  decompress only: {decomp_time:.2f} s  "
          f"({num_reads / decomp_time / 1e6:.2f} M lines/s)")

    print(f"\n[3/x] Benchmarking extract_ids only (single-threaded, no pool) …")
    parse_time, n = bench_extract_ids_only(bgz_path)
    print(f"  extract_ids only: {parse_time:.2f} s  ({n / parse_time / 1e6:.2f} M lines/s)")

    results = {}
    for i, nw in enumerate(workers, start=4):
        label = f"get_counts (workers={nw})"
        print(f"\n[{i}/x] Benchmarking {label} …")
        t = bench_get_counts(bgz_path, num_reads, nw)
        results[nw] = t
        print(f"  {label}: {t:.2f} s  ({num_reads / t / 1e6:.2f} M reads/s)")

    decomp_offset_results = {}
    for i, nw in enumerate(workers, start=len(workers) + 4):
        label = f"decompress+offset get_counts (workers={nw})"
        print(f"\n[{i}/x] Benchmarking {label} …")
        t = bench_decompress_then_offset(bgz_path, nw)
        decomp_offset_results[nw] = t
        print(f"  {label}: {t:.2f} s  ({num_reads / t / 1e6:.2f} M reads/s)")

    next_step = len(workers) * 2 + 4
    print(f"\n[{next_step}/x] Benchmarking Rust get_counts …")
    rust_time, rust_n = bench_get_counts_rust(bgz_path)
    print(f"  Rust get_counts: {rust_time:.2f} s  ({rust_n / rust_time / 1e6:.2f} M reads/s)")

    # --- summary ---
    single = results.get(1)
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Reads:              {num_reads:,}")
    print(f"  File size:          {bgz_path.stat().st_size / 1e6:.1f} MB")
    print(f"  Decompress only:    {decomp_time:.2f} s")
    print(f"  Extract IDs only:   {parse_time:.2f} s  (decomp ratio {parse_time / decomp_time:.2f}x)")
    for nw, t in results.items():
        ratio = t / decomp_time
        print(f"  get_counts w={nw:<3}:  {t:.2f} s  (decomp ratio {ratio:.2f}x)")
    for nw, t in decomp_offset_results.items():
        ratio = t / decomp_time
        print(f"  decomp+offset w={nw:<3}: {t:.2f} s  (decomp ratio {ratio:.2f}x)")
    if single is not None:
        print(f"  Rust get_counts:    {rust_time:.2f} s  ({single / rust_time:.1f}x vs Python w=1)")
    else:
        print(f"  Rust get_counts:    {rust_time:.2f} s")

    if single is not None:
        ratio = decomp_time / single
        print()
        if ratio > 0.5:
            print("DIAGNOSIS: decompression is ≥50% of single-worker time.")
            print("  → block-parallel decompression (pysam BGZFile) would likely help.")
        else:
            print("DIAGNOSIS: decompression is <50% of single-worker time.")
            print("  → parsing is the bottleneck; current isal+multiprocessing is sufficient.")
    print("=" * 60)

    # --- plain-text section: decompress bgz → txt, rerun benches without decompression ---
    txt_path = bgz_path.with_suffix(".txt")
    print(f"\nDecompressing to plain text for no-decompress benchmarks …")
    subprocess.run(["bgzip", "-d", "-k", "-@", str(cpu_count()), "-c", str(bgz_path)],
                   stdout=open(txt_path, "wb"), check=True)
    print(f"  Plain text size: {txt_path.stat().st_size / 1e6:.1f} MB")

    try:
        print(f"\n[p1] Plain read-only (no parse, no decompress) …")
        plain_read = bench_plain_read_only(txt_path)
        print(f"  read-only: {plain_read:.2f} s  ({num_reads / plain_read / 1e6:.2f} M lines/s)")

        print(f"\n[p2] Plain extract_ids only (1 thread, no decompress) …")
        plain_parse, _ = bench_plain_extract_ids_only(txt_path)
        print(f"  extract_ids: {plain_parse:.2f} s  ({num_reads / plain_parse / 1e6:.2f} M lines/s)")

        plain_results = {}
        for i, nw in enumerate(workers, start=3):
            label = f"plain get_counts (workers={nw})"
            print(f"\n[p{i}] Benchmarking {label} …")
            t = bench_plain_get_counts(txt_path, nw)
            plain_results[nw] = t
            print(f"  {label}: {t:.2f} s  ({num_reads / t / 1e6:.2f} M reads/s)")

        offset_results = {}
        for i, nw in enumerate(workers, start=3 + len(workers)):
            label = f"plain get_counts offset (workers={nw})"
            print(f"\n[p{i}] Benchmarking {label} …")
            t = bench_plain_get_counts_offset(txt_path, nw)
            offset_results[nw] = t
            print(f"  {label}: {t:.2f} s  ({num_reads / t / 1e6:.2f} M reads/s)")

        plain_single = plain_results.get(1)
        print("\n" + "=" * 60)
        print("SUMMARY — plain text (no decompression)")
        print("=" * 60)
        print(f"  read-only:              {plain_read:.2f} s")
        print(f"  extract_ids only:       {plain_parse:.2f} s")
        for nw, t in plain_results.items():
            speedup = f"  ({plain_single / t:.1f}x vs w=1)" if plain_single and nw != 1 else ""
            print(f"  batched   w={nw:<3}:      {t:.2f} s  ({num_reads / t / 1e6:.2f} M reads/s){speedup}")
        offset_single = offset_results.get(1)
        for nw, t in offset_results.items():
            speedup = f"  ({offset_single / t:.1f}x vs w=1)" if offset_single and nw != 1 else ""
            print(f"  offset    w={nw:<3}:      {t:.2f} s  ({num_reads / t / 1e6:.2f} M reads/s){speedup}")
        if plain_single:
            print(f"\n  Decompress overhead (isal): +{single - plain_single:.2f} s vs plain batched w=1")
        print("=" * 60)
    finally:
        txt_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
