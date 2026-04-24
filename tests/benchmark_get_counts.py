"""Benchmark decompression, parsing, and end-to-end get_counts performance.

Synthetic format mirrors real cutadapt output:
  @read{i} x?S0.{s0}?B1.{b1}?B2.{b2}?B3.{b3}?S1.{s1}

Run with:
  uv run python tests/benchmark_get_counts.py
  uv run python tests/benchmark_get_counts.py --n-reads 100000
"""
import argparse
import gzip
import multiprocessing
import random
import tempfile
import time
from collections.abc import Iterable
from pathlib import Path

from isal import igzip

from delt_hit.demultiplex.postprocess import extract_ids, get_counts

DEFAULT_N_READS = 5_000_000
N_S = 10      # unique values per selection primer
N_BB = 100    # unique values per building-block region
SEED = 42


def generate_synthetic_file(path: Path, n_reads: int = DEFAULT_N_READS) -> None:
    rng = random.Random(SEED)
    with gzip.open(path, 'wb', compresslevel=1) as f:
        for i in range(n_reads):
            s0 = rng.randint(0, N_S - 1)
            s1 = rng.randint(0, N_S - 1)
            b1 = rng.randint(0, N_BB - 1)
            b2 = rng.randint(0, N_BB - 1)
            b3 = rng.randint(0, N_BB - 1)
            line = f'@read{i} x?S0.{s0}?B1.{b1}?B2.{b2}?B3.{b3}?S1.{s1}\n'
            f.write(line.encode())


def normalise(counts: dict) -> dict:
    """Convert nested defaultdicts / dicts to plain nested dicts for comparison."""
    return {k: dict(v) for k, v in counts.items()}


def bench_gzip_decompress_only(path: Path) -> tuple[float, int]:
    """Measure stdlib gzip decompression throughput without parsing."""
    t0 = time.perf_counter()
    n_lines = 0
    with gzip.open(path, "rt") as handle:
        for _ in handle:
            n_lines += 1
    return time.perf_counter() - t0, n_lines


def bench_isal_decompress_only(path: Path) -> tuple[float, int]:
    """Measure isal decompression throughput without parsing."""
    t0 = time.perf_counter()
    n_lines = 0
    with igzip.open(path, "rt") as handle:
        for _ in handle:
            n_lines += 1
    return time.perf_counter() - t0, n_lines


def bench_parse_only(lines: Iterable[str]) -> tuple[float, int]:
    """Measure parsing cost after decompression has already happened."""
    t0 = time.perf_counter()
    n_lines = 0
    for line in lines:
        extract_ids(line)
        n_lines += 1
    return time.perf_counter() - t0, n_lines


def decompress_to_plain_text(src_path: Path, dst_path: Path) -> int:
    """Materialize the compressed benchmark input as plain text once."""
    n_lines = 0
    with gzip.open(src_path, "rt") as src, open(dst_path, "wt") as dst:
        for line in src:
            dst.write(line)
            n_lines += 1
    return n_lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-reads", type=int, default=DEFAULT_N_READS)
    args = parser.parse_args()

    n_reads = args.n_reads
    n_workers = multiprocessing.cpu_count()

    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / 'reads_with_adapters.gz'
        plain_path = Path(tmpdir) / 'reads_with_adapters.txt'

        print(f"Generating {n_reads:,} synthetic reads …")
        t0 = time.perf_counter()
        generate_synthetic_file(path, n_reads)
        size_mb = path.stat().st_size / 1e6
        print(f"  done in {time.perf_counter() - t0:.1f}s  ({size_mb:.0f} MB compressed)\n")

        print("Decompression only (gzip) …")
        t_gzip_decompress, n_gzip_lines = bench_gzip_decompress_only(path)
        print(f"  {t_gzip_decompress:.2f}s  ({n_gzip_lines / t_gzip_decompress / 1e6:.2f} M lines/s)\n")

        print("Decompression only (isal) …")
        t_isal_decompress, n_isal_lines = bench_isal_decompress_only(path)
        print(f"  {t_isal_decompress:.2f}s  ({n_isal_lines / t_isal_decompress / 1e6:.2f} M lines/s)\n")

        print("Materializing plain-text benchmark input …")
        n_plain_lines = decompress_to_plain_text(path, plain_path)
        print(f"  wrote {n_plain_lines:,} lines\n")

        print("Parsing only (extract_ids on preloaded lines) …")
        with open(plain_path, "rt") as handle:
            t_parse_only, n_parse_lines = bench_parse_only(handle)
        print(f"  {t_parse_only:.2f}s  ({n_parse_lines / t_parse_only / 1e6:.2f} M lines/s)\n")

        print("Serial (num_workers=1) …")
        t0 = time.perf_counter()
        counts_serial = get_counts(input_path=path, num_reads=n_reads, num_workers=1)
        t_serial = time.perf_counter() - t0
        print(f"  {t_serial:.2f}s\n")

        print(f"Parallel (num_workers={n_workers}) …")
        t0 = time.perf_counter()
        counts_parallel = get_counts(input_path=path, num_reads=n_reads, num_workers=n_workers)
        t_parallel = time.perf_counter() - t0
        print(f"  {t_parallel:.2f}s")
        print(f"  speedup: {t_serial / t_parallel:.1f}x\n")

        print("Verifying output equality …")
        n_serial = normalise(counts_serial)
        n_parallel = normalise(counts_parallel)

        assert set(n_serial) == set(n_parallel), (
            f"Selection key mismatch: serial={set(n_serial)-set(n_parallel)} "
            f"parallel={set(n_parallel)-set(n_serial)}"
        )
        for sel in n_serial:
            assert n_serial[sel] == n_parallel[sel], (
                f"Count mismatch for selection {sel}"
            )

        total_reads = sum(sum(v.values()) for v in n_serial.values())
        print(f"  ✓ identical  ({total_reads:,} reads across {len(n_serial)} selections)")

        print("\nTiming summary")
        print(f"  gzip decompress only: {t_gzip_decompress:.2f}s")
        print(f"  isal decompress only: {t_isal_decompress:.2f}s")
        print(f"  parse only:           {t_parse_only:.2f}s")
        print(f"  serial get_counts:    {t_serial:.2f}s")
        print(f"  parallel get_counts:  {t_parallel:.2f}s")


if __name__ == '__main__':
    main()
