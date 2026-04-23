"""Benchmark serial vs parallel get_counts on a 5M synthetic reads file.

Synthetic format mirrors real cutadapt output:
  @read{i} x?S0.{s0}?B1.{b1}?B2.{b2}?B3.{b3}?S1.{s1}

Run with:
  uv run python tests/benchmark_get_counts.py
"""
import gzip
import multiprocessing
import random
import tempfile
import time
from pathlib import Path

from delt_hit.demultiplex.postprocess import get_counts

N_READS = 5_000_000
N_S = 10      # unique values per selection primer
N_BB = 100    # unique values per building-block region
SEED = 42


def generate_synthetic_file(path: Path, n_reads: int = N_READS) -> None:
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


def main():
    n_workers = multiprocessing.cpu_count()

    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / 'reads_with_adapters.gz'

        print(f"Generating {N_READS:,} synthetic reads …")
        t0 = time.perf_counter()
        generate_synthetic_file(path, N_READS)
        size_mb = path.stat().st_size / 1e6
        print(f"  done in {time.perf_counter() - t0:.1f}s  ({size_mb:.0f} MB compressed)\n")

        print("Serial (num_workers=1) …")
        t0 = time.perf_counter()
        counts_serial = get_counts(input_path=path, num_reads=N_READS, num_workers=1)
        t_serial = time.perf_counter() - t0
        print(f"  {t_serial:.2f}s\n")

        print(f"Parallel (num_workers={n_workers}) …")
        t0 = time.perf_counter()
        counts_parallel = get_counts(input_path=path, num_reads=N_READS, num_workers=n_workers)
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


if __name__ == '__main__':
    main()
