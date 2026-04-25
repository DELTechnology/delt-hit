"""Benchmark decompression, parsing, and end-to-end get_counts performance."""
import argparse
import multiprocessing
import tempfile
import time
from collections.abc import Iterable
from pathlib import Path

from isal import igzip

from delt_hit.demultiplex.postprocess import extract_ids, get_counts

DEFAULT_REPORT_DIR = Path("target/benchmarks/get_counts")

def normalise(counts: dict) -> dict:
    """Convert nested defaultdicts / dicts to plain nested dicts for comparison."""
    return {k: dict(v) for k, v in counts.items()}


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
    with igzip.open(src_path, "rt") as src, open(dst_path, "wt") as dst:
        for line in src:
            dst.write(line)
            n_lines += 1
    return n_lines


def report_path(*, input_path: Path, num_workers: int, chunk_size_bytes: int, report_dir: Path) -> Path:
    input_name = input_path.name.removesuffix(".gz")
    chunk_size_millions = chunk_size_bytes / 1_000_000
    filename = (
        f"input={input_name}-"
        f"num_workers={num_workers}-"
        f"chunk_size={chunk_size_millions:g}.txt"
    )
    return report_dir / filename


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-path", type=Path, required=True, help="Path to a pre-generated reads_with_adapters.gz input.")
    parser.add_argument("--num-workers", type=int, required=True, help="Worker count for the parallel get_counts benchmark.")
    parser.add_argument("--chunk-size-bytes", type=int, default=5_000_000, help="Chunk size passed to parallel get_counts.")
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR, help="Directory to write the short text report.")
    parser.add_argument("--run-decompress-only", action="store_true", help="Run the decompression-only benchmark.")
    parser.add_argument("--run-parse-only", action="store_true", help="Run the parse-only benchmark.")
    parser.add_argument("--run-serial", action="store_true", help="Run the serial get_counts benchmark.")
    args = parser.parse_args()

    n_workers = args.num_workers

    with tempfile.TemporaryDirectory() as tmpdir:
        path = args.input_path
        plain_path = Path(tmpdir) / 'reads_with_adapters.txt'
        size_mb = path.stat().st_size / 1e6
        print(f"Using pre-generated input {path} ({size_mb:.0f} MB compressed)\n")

        t_isal_decompress = None
        t_parse_only = None
        t_serial = None

        n_reads = sum(1 for _ in igzip.open(path, "rt"))

        if args.run_decompress_only:
            print("Decompression only (isal) …")
            t_isal_decompress, n_isal_lines = bench_isal_decompress_only(path)
            print(f"  {t_isal_decompress:.2f}s  ({n_isal_lines / t_isal_decompress / 1e6:.2f} M lines/s)\n")

        if args.run_parse_only:
            print("Materializing plain-text benchmark input …")
            n_plain_lines = decompress_to_plain_text(path, plain_path)
            print(f"  wrote {n_plain_lines:,} lines\n")

            print("Parsing only (extract_ids on preloaded lines) …")
            with open(plain_path, "rt") as handle:
                t_parse_only, n_parse_lines = bench_parse_only(handle)
            print(f"  {t_parse_only:.2f}s  ({n_parse_lines / t_parse_only / 1e6:.2f} M lines/s)\n")

        counts_serial = None
        if args.run_serial:
            print("Serial (num_workers=1) …")
            t0 = time.perf_counter()
            counts_serial = get_counts(input_path=path, num_reads=n_reads, num_workers=1)
            t_serial = time.perf_counter() - t0
            print(f"  {t_serial:.2f}s\n")

        print(f"Parallel (num_workers={n_workers}) …")
        t0 = time.perf_counter()
        counts_parallel = get_counts(
            input_path=path,
            num_reads=n_reads,
            num_workers=n_workers,
            chunk_size_bytes=args.chunk_size_bytes,
        )
        t_parallel = time.perf_counter() - t0
        print(f"  {t_parallel:.2f}s")
        if t_serial is not None:
            print(f"  speedup: {t_serial / t_parallel:.1f}x")
        print("")

        if counts_serial is not None:
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
            counts_match = True
            total_reads = sum(sum(v.values()) for v in n_serial.values())
            print(f"  ✓ identical  ({total_reads:,} reads across {len(n_serial)} selections)")
        else:
            counts_match = None
            total_reads = sum(sum(v.values()) for v in counts_parallel.values())

        print("\nTiming summary")
        if t_isal_decompress is not None:
            print(f"  isal decompress only: {t_isal_decompress:.2f}s")
        if t_parse_only is not None:
            print(f"  parse only:           {t_parse_only:.2f}s")
        if t_serial is not None:
            print(f"  serial get_counts:    {t_serial:.2f}s")
        print(f"  chunk size bytes:     {args.chunk_size_bytes}")
        print(f"  parallel get_counts:  {t_parallel:.2f}s")

        summary = [
            f"input_path={path}",
            f"input_size_mb={size_mb:.0f}",
            f"num_reads={n_reads}",
            f"num_workers={n_workers}",
            f"chunk_size_bytes={args.chunk_size_bytes}",
            f"parallel_get_counts_s={t_parallel:.2f}",
            f"total_reads={total_reads}",
        ]
        if t_isal_decompress is not None:
            summary.append(f"decompress_only_s={t_isal_decompress:.2f}")
        if t_parse_only is not None:
            summary.append(f"parse_only_s={t_parse_only:.2f}")
        if t_serial is not None:
            summary.append(f"serial_get_counts_s={t_serial:.2f}")
            summary.append(f"speedup={t_serial / t_parallel:.2f}")
        if counts_match is not None:
            summary.append(f"counts_match={counts_match}")
        args.report_dir.mkdir(parents=True, exist_ok=True)
        out_path = report_path(
            input_path=path,
            num_workers=n_workers,
            chunk_size_bytes=args.chunk_size_bytes,
            report_dir=args.report_dir,
        )
        out_path.write_text("\n".join(summary) + "\n")
        print(f"\nWrote report to {out_path}")


if __name__ == '__main__':
    main()
