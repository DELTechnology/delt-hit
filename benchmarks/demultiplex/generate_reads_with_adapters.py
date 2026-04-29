#!/usr/bin/env python3
"""Generate synthetic reads_with_adapters.gz files for get_counts benchmarks."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from isal import igzip


DEFAULT_N_READS = 5_000_000
N_S = 10
N_BB = 100
SEED = 42


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-path", type=Path, required=True, help="Path to the output reads_with_adapters.gz file.")
    parser.add_argument("--n-reads", type=int, default=DEFAULT_N_READS, help="Number of synthetic reads to generate.")
    return parser.parse_args()


def generate_synthetic_file(path: Path, n_reads: int) -> None:
    rng = random.Random(SEED)
    path.parent.mkdir(parents=True, exist_ok=True)
    with igzip.open(path, "wb", compresslevel=1) as handle:
        for i in range(n_reads):
            s0 = rng.randint(0, N_S - 1)
            s1 = rng.randint(0, N_S - 1)
            b1 = rng.randint(0, N_BB - 1)
            b2 = rng.randint(0, N_BB - 1)
            b3 = rng.randint(0, N_BB - 1)
            line = f"@read{i} x?S0.{s0}?B1.{b1}?B2.{b2}?B3.{b3}?S1.{s1}\n"
            handle.write(line.encode())


def main() -> None:
    args = parse_args()
    generate_synthetic_file(args.output_path, args.n_reads)
    size_mb = args.output_path.stat().st_size / 1e6
    print(f"Wrote {args.n_reads:,} reads to {args.output_path} ({size_mb:.0f} MB compressed)")


if __name__ == "__main__":
    main()
