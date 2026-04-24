#!/usr/bin/env python3
"""Generate a configurable synthetic FASTQ dataset plus truth artifacts.

The generator emits every combinatorial compound for ``n`` cycles with ``m``
building blocks per cycle, repeated ``q`` times per compound. Total reads are
``m ** n * q``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from itertools import product
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "target" / "synthetic"
DEFAULT_EXPERIMENT_NAME = "synthetic"
DEFAULT_NUM_CYCLES = 2
DEFAULT_BUILDING_BLOCKS_PER_CYCLE = 10
DEFAULT_NUM_READS_PER_COMPOUND = 1
DEFAULT_COMPRESSED = True

LIBRARY_TAG = "ACGTACGTAC"
CLOSING_TAG = "TTGGAACC"
TAG_LENGTH = 8
UMI_LENGTH = 8
QUALITY_CHAR = "I"


def int_to_dna(value: int, length: int = TAG_LENGTH) -> str:
    """Convert a non-negative integer to a fixed-width DNA tag."""
    alphabet = "ACGT"
    bases: list[str] = []
    for _ in range(length):
        bases.append(alphabet[value % 4])
        value //= 4
    return "".join(reversed(bases))


def expected_compounds(*, num_cycles: int, building_blocks_per_cycle: int) -> int:
    return building_blocks_per_cycle ** num_cycles


def expected_reads(
    *,
    num_cycles: int,
    building_blocks_per_cycle: int,
    num_reads_per_compound: int,
) -> int:
    return expected_compounds(
        num_cycles=num_cycles,
        building_blocks_per_cycle=building_blocks_per_cycle,
    ) * num_reads_per_compound


def make_building_blocks(*, num_cycles: int, building_blocks_per_cycle: int) -> list[list[dict[str, object]]]:
    building_blocks_by_cycle: list[list[dict[str, object]]] = []
    for cycle in range(1, num_cycles + 1):
        cycle_rows: list[dict[str, object]] = []
        offset = (cycle - 1) * 10_000
        for code in range(1, building_blocks_per_cycle + 1):
            cycle_rows.append(
                {
                    "cycle": cycle,
                    "code": code,
                    "building_block_id": f"BB{cycle}_{code:03d}",
                    "tag": int_to_dna(offset + code - 1),
                }
            )
        building_blocks_by_cycle.append(cycle_rows)
    return building_blocks_by_cycle


def write_building_blocks(path: Path, building_blocks_by_cycle: list[list[dict[str, object]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["cycle", "code", "building_block_id", "tag"],
            delimiter="\t",
        )
        writer.writeheader()
        for cycle_rows in building_blocks_by_cycle:
            writer.writerows(cycle_rows)


def write_expected_counts(
    path: Path,
    building_blocks_by_cycle: list[list[dict[str, object]]],
    *,
    num_reads_per_compound: int,
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    num_cycles = len(building_blocks_by_cycle)
    fieldnames = []
    for cycle in range(1, num_cycles + 1):
        fieldnames.append(f"code_{cycle}")
    for cycle in range(1, num_cycles + 1):
        fieldnames.append(f"building_block_id_{cycle}")
    fieldnames.append("count")

    row_count = 0
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for compound in product(*building_blocks_by_cycle):
            row_count += 1
            row = {f"code_{idx}": block["code"] for idx, block in enumerate(compound, start=1)}
            row.update(
                {
                    f"building_block_id_{idx}": block["building_block_id"]
                    for idx, block in enumerate(compound, start=1)
                }
            )
            row["count"] = num_reads_per_compound
            writer.writerow(row)
    return row_count


def write_fastq(
    path: Path,
    building_blocks_by_cycle: list[list[dict[str, object]]],
    *,
    num_reads_per_compound: int,
    experiment_name: str,
    compressed: bool,
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    read_count = 0
    open_func = gzip.open if compressed else open
    mode = "wt"
    with open_func(path, mode) as handle:
        for compound_idx, compound in enumerate(product(*building_blocks_by_cycle), start=1):
            tags = "".join(str(block["tag"]) for block in compound)
            ids = "_".join(str(block["building_block_id"]) for block in compound)
            for replicate_idx in range(1, num_reads_per_compound + 1):
                read_count += 1
                umi = int_to_dna(20_000 + read_count, UMI_LENGTH)
                sequence = f"{LIBRARY_TAG}{tags}{umi}{CLOSING_TAG}"
                quality = QUALITY_CHAR * len(sequence)
                handle.write(
                    f"@{experiment_name}_compound_{compound_idx:06d}_read_{replicate_idx:03d}_{ids}\n"
                    f"{sequence}\n"
                    "+\n"
                    f"{quality}\n"
                )
    return read_count


def write_manifest(
    path: Path,
    *,
    experiment_name: str,
    num_cycles: int,
    building_blocks_per_cycle: int,
    num_reads_per_compound: int,
    output_dir: Path,
    fastq_path: Path,
    expected_counts_path: Path,
    building_blocks_path: Path,
) -> None:
    manifest = {
        "experiment_name": experiment_name,
        "output_dir": str(output_dir.resolve()),
        "num_cycles": num_cycles,
        "building_blocks_per_cycle": building_blocks_per_cycle,
        "num_reads_per_compound": num_reads_per_compound,
        "expected_compounds": expected_compounds(
            num_cycles=num_cycles,
            building_blocks_per_cycle=building_blocks_per_cycle,
        ),
        "expected_reads": expected_reads(
            num_cycles=num_cycles,
            building_blocks_per_cycle=building_blocks_per_cycle,
            num_reads_per_compound=num_reads_per_compound,
        ),
        "library_tag": LIBRARY_TAG,
        "closing_tag": CLOSING_TAG,
        "tag_length": TAG_LENGTH,
        "umi_length": UMI_LENGTH,
        "files": {
            "fastq": str(fastq_path.resolve()),
            "expected_counts": str(expected_counts_path.resolve()),
            "building_blocks": str(building_blocks_path.resolve()),
        },
        "compressed_fastq": fastq_path.suffix == ".gz",
    }
    path.write_text(json.dumps(manifest, indent=2) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--num-cycles",
        type=int,
        default=DEFAULT_NUM_CYCLES,
        help="Number of barcode cycles to encode in each read.",
    )
    parser.add_argument(
        "--building-blocks-per-cycle",
        type=int,
        default=DEFAULT_BUILDING_BLOCKS_PER_CYCLE,
        help="Uniform number of building blocks to generate for every cycle.",
    )
    parser.add_argument(
        "--num-reads-per-compound",
        type=int,
        default=DEFAULT_NUM_READS_PER_COMPOUND,
        help="Replicate reads to emit for each combinatorial compound.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where the experiment folder is created.",
    )
    parser.add_argument(
        "--experiment-name",
        default=DEFAULT_EXPERIMENT_NAME,
        help="Name used for the experiment folder and FASTQ read identifiers.",
    )
    parser.add_argument(
        "--compressed",
        default=str(DEFAULT_COMPRESSED).lower(),
        choices=("true", "false"),
        help="Whether to write the FASTQ as gzip compressed output.",
    )
    args = parser.parse_args()
    if args.num_cycles < 1:
        parser.error("num_cycles must be at least 1")
    if args.building_blocks_per_cycle < 1:
        parser.error("building_blocks_per_cycle must be at least 1")
    if args.num_reads_per_compound < 1:
        parser.error("num_reads_per_compound must be at least 1")
    args.compressed = args.compressed == "true"
    return args


def main(
    *,
    num_cycles: int = DEFAULT_NUM_CYCLES,
    building_blocks_per_cycle: int = DEFAULT_BUILDING_BLOCKS_PER_CYCLE,
    num_reads_per_compound: int = DEFAULT_NUM_READS_PER_COMPOUND,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    experiment_name: str = DEFAULT_EXPERIMENT_NAME,
    compressed: bool = DEFAULT_COMPRESSED,
) -> None:
    if num_cycles < 1:
        raise ValueError("num_cycles must be at least 1")
    if building_blocks_per_cycle < 1:
        raise ValueError("building_blocks_per_cycle must be at least 1")
    if num_reads_per_compound < 1:
        raise ValueError("num_reads_per_compound must be at least 1")

    experiment_dir = output_dir / experiment_name
    experiment_dir.mkdir(parents=True, exist_ok=True)

    fastq_suffix = ".fastq.gz" if compressed else ".fastq"
    fastq_path = experiment_dir / f"{experiment_name}{fastq_suffix}"
    expected_counts_path = experiment_dir / "expected_counts.tsv"
    building_blocks_path = experiment_dir / "building_blocks.tsv"
    manifest_path = experiment_dir / "manifest.json"

    building_blocks_by_cycle = make_building_blocks(
        num_cycles=num_cycles,
        building_blocks_per_cycle=building_blocks_per_cycle,
    )
    write_building_blocks(building_blocks_path, building_blocks_by_cycle)
    compound_count = write_expected_counts(
        expected_counts_path,
        building_blocks_by_cycle,
        num_reads_per_compound=num_reads_per_compound,
    )
    read_count = write_fastq(
        fastq_path,
        building_blocks_by_cycle,
        num_reads_per_compound=num_reads_per_compound,
        experiment_name=experiment_name,
        compressed=compressed,
    )
    write_manifest(
        manifest_path,
        experiment_name=experiment_name,
        num_cycles=num_cycles,
        building_blocks_per_cycle=building_blocks_per_cycle,
        num_reads_per_compound=num_reads_per_compound,
        output_dir=experiment_dir,
        fastq_path=fastq_path,
        expected_counts_path=expected_counts_path,
        building_blocks_path=building_blocks_path,
    )

    print(f"Wrote FASTQ with {read_count} reads to {fastq_path}")
    print(f"Wrote expected counts for {compound_count} compounds to {expected_counts_path}")
    print(f"Wrote building block metadata to {building_blocks_path}")
    print(f"Wrote manifest to {manifest_path}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(
        num_cycles=cli_args.num_cycles,
        building_blocks_per_cycle=cli_args.building_blocks_per_cycle,
        num_reads_per_compound=cli_args.num_reads_per_compound,
        output_dir=cli_args.output_dir,
        experiment_name=cli_args.experiment_name,
        compressed=cli_args.compressed,
    )
