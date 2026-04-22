#!/usr/bin/env python3
"""Generate small synthetic DELi decoding experiments.

The fixtures contain two-, three-, and four-cycle DEL libraries with:
- 100 cycle-1 building blocks
- 10 building blocks for each additional cycle
- configurable FASTQ read depth per encoded compound combination
"""

from __future__ import annotations

import argparse
import csv
import json
from itertools import product
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "other_tools" / "DELi"
BARCODE_COUNTS = (2, 3, 4)
BUILDING_BLOCK_COUNTS_BY_CYCLE = {
    1: 100,
    2: 10,
    3: 10,
    4: 10,
}
DEFAULT_NUM_READS_PER_COMPOUND = 1
TAG_LENGTH = 8
UMI_LENGTH = 8

LIBRARY_TAG = "ACGTACGTAC"
CLOSING_TAG = "TTGGAACC"
QUALITY_CHAR = "I"


def int_to_dna(value: int, length: int = TAG_LENGTH) -> str:
    """Convert a non-negative integer to a fixed-width DNA tag."""
    alphabet = "ACGT"
    bases: list[str] = []
    for _ in range(length):
        bases.append(alphabet[value % 4])
        value //= 4
    return "".join(reversed(bases))


def make_building_blocks(prefix: str, count: int, offset: int) -> list[dict[str, str]]:
    """Create deterministic, unique tagged building blocks."""
    return [
        {
            "id": f"{prefix}{idx + 1:03d}",
            "tag": int_to_dna(offset + idx),
        }
        for idx in range(count)
    ]


def write_building_block_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "tag"])
        writer.writeheader()
        writer.writerows(rows)


def compound_count(building_block_rows: list[list[dict[str, str]]]) -> int:
    """Return the number of encoded compounds in a combinatorial library."""
    total = 1
    for rows in building_block_rows:
        total *= len(rows)
    return total


def write_fastq(
    path: Path,
    building_block_rows: list[list[dict[str, str]]],
    *,
    num_reads_per_compound: int,
) -> int:
    """Write one or more FASTQ reads for every building-block combination."""
    path.parent.mkdir(parents=True, exist_ok=True)
    read_idx = 0
    with path.open("w") as handle:
        for compound_idx, blocks in enumerate(product(*building_block_rows), start=1):
            block_ids = "_".join(block["id"] for block in blocks)
            tags = "".join(block["tag"] for block in blocks)
            for replicate_idx in range(1, num_reads_per_compound + 1):
                read_idx += 1
                umi = int_to_dna(20_000 + read_idx, UMI_LENGTH)
                sequence = f"{LIBRARY_TAG}{tags}{umi}{CLOSING_TAG}"
                quality = QUALITY_CHAR * len(sequence)
                handle.write(
                    f"@synthetic_combo_{compound_idx:04d}_read_{replicate_idx:03d}_{block_ids}\n"
                    f"{sequence}\n"
                    "+\n"
                    f"{quality}\n"
                )
    return read_idx


def write_deli_config(path: Path, data_dir: Path) -> None:
    """Write a local DELi config file compatible with DELi 0.2.x."""
    path.write_text(
        "\n".join(
            [
                "[DEFAULT]",
                "",
                "[deli.data]",
                f"deli_data_dir = {data_dir.resolve()}",
                "",
                "[deli.hamming]",
                "nuc_2_int = A:0,T:1,C:2,G:3",
                "",
                "[deli.hamming.8_4]",
                "hamming_order = p0,p1,p2,d3,p4,d5,d6,d7",
                "custom_order = p0,p1,p2,d3,p4,d5,d6,d7",
                "",
                "[deli.hamming.16_5]",
                "hamming_order = p0,p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15",
                "custom_order = p0,p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15",
                "",
                "[deli.hamming.7_3]",
                "hamming_order = p1,p2,d3,p4,d5,d6,d7",
                "custom_order = p1,p2,d3,p4,d5,d6,d7",
                "",
                "[deli.hamming.15_4]",
                "hamming_order = p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15",
                "custom_order = p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15",
                "",
                "[deli.buildingblocks]",
                "BB_MASK = ###",
                "",
            ]
        )
    )


def make_library(
    *,
    barcode_count: int,
    building_block_csvs: list[Path],
) -> dict:
    return {
        "dna_barcode_on": None,
        "bb_sets": [
            {"cycle": cycle, "bb_set_path": str(path)}
            for cycle, path in enumerate(building_block_csvs[:barcode_count], start=1)
        ],
        "barcode_schema": {
            "library": {"tag": LIBRARY_TAG},
            **{f"bb{cycle}": {"tag": "N" * TAG_LENGTH} for cycle in range(1, barcode_count + 1)},
            "umi": {"tag": "N" * UMI_LENGTH},
            "closing": {"tag": CLOSING_TAG},
        },
    }


def make_decode_config(
    *,
    barcode_count: int,
    fastq_path: Path,
    library_json: Path,
    min_read_length: int,
    num_reads_per_compound: int,
) -> dict:
    compound_reads = f"{num_reads_per_compound} read"
    if num_reads_per_compound != 1:
        compound_reads += "s"

    return {
        "selection_id": f"synthetic_{barcode_count}cycle_selection",
        "target_id": "synthetic_target",
        "selection_condition": "synthetic_smoke_test",
        "additional_info": f"Every {barcode_count}-barcode compound combination appears with {compound_reads}.",
        "data_ran": "2026-04-20T00:00:00",
        "sequence_files": [str(fastq_path)],
        "libraries": [str(library_json)],
        "decode_settings": {
            "library_error_tolerance": 0.0,
            "min_library_overlap": len(LIBRARY_TAG),
            "alignment_algorithm": "semi",
            "bb_calling_approach": "bio",
            "revcomp": False,
            "max_read_length": 128,
            "min_read_length": min_read_length,
            "read_type": "single",
            "disable_error_correction": True,
            "umi_clustering": False,
            "umi_min_distance": 2,
        },
    }


def write_decode_settings_v02(path: Path, *, min_read_length: int) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "ignore_tool_compounds": True,
                "demultiplexer_algorithm": "regex",
                "demultiplexer_mode": "library",
                "realign": False,
                "library_error_tolerance": 0,
                "library_error_correction_mode_str": "disable",
                "min_library_overlap": len(LIBRARY_TAG),
                "revcomp": False,
                "library_wiggle": False,
                "wiggle": False,
                "decode_matching_approach": "first_perfect",
                "max_read_length": 128,
                "min_read_length": min_read_length,
                "default_error_correction_mode_str": "disable",
            },
            sort_keys=False,
        )
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "settings",
        nargs="*",
        help="Optional key=value settings. Supported: num_reads_per_compound=<int>.",
    )
    parser.add_argument(
        "--num-reads-per-compound",
        type=int,
        default=None,
        help="FASTQ reads to emit for each building-block combination.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Directory where other_tools/DELi-style experiment folders are written.",
    )
    args = parser.parse_args()

    for setting in args.settings:
        key, sep, value = setting.partition("=")
        if sep != "=" or key != "num_reads_per_compound":
            parser.error(f"unsupported setting {setting!r}; use num_reads_per_compound=<int>")
        if args.num_reads_per_compound is not None:
            parser.error("set num_reads_per_compound only once")
        try:
            args.num_reads_per_compound = int(value)
        except ValueError:
            parser.error("num_reads_per_compound must be an integer")

    if args.num_reads_per_compound is None:
        args.num_reads_per_compound = DEFAULT_NUM_READS_PER_COMPOUND
    if args.num_reads_per_compound < 1:
        parser.error("num_reads_per_compound must be at least 1")
    return args


def write_experiment(
    *,
    output_root: Path,
    barcode_count: int,
    building_blocks_by_cycle: dict[int, list[dict[str, str]]],
    num_reads_per_compound: int,
) -> int:
    experiment_name = f"synthetic_{barcode_count}cycle"
    experiment_dir = output_root / experiment_name
    bb_dir = experiment_dir / "data" / "building_blocks"
    lib_dir = experiment_dir / "data" / "libraries"
    input_dir = experiment_dir / "input"
    output_dir = experiment_dir / "output_v02"

    for path in [
        bb_dir,
        lib_dir,
        experiment_dir / "data" / "reactions",
        experiment_dir / "data" / "tool_compounds",
        input_dir,
        output_dir,
    ]:
        path.mkdir(parents=True, exist_ok=True)

    building_block_csvs = []
    for cycle in range(1, barcode_count + 1):
        csv_path = bb_dir / f"synthetic_bb{cycle}.csv"
        write_building_block_csv(csv_path, building_blocks_by_cycle[cycle])
        building_block_csvs.append(csv_path)

    fastq_path = input_dir / f"{experiment_name}.fastq"
    library_json = lib_dir / f"{experiment_name}_library.json"
    min_read_length = len(LIBRARY_TAG) + barcode_count * TAG_LENGTH + UMI_LENGTH

    cycle_rows = [building_blocks_by_cycle[cycle] for cycle in range(1, barcode_count + 1)]
    read_count = write_fastq(
        fastq_path,
        cycle_rows,
        num_reads_per_compound=num_reads_per_compound,
    )

    library = make_library(barcode_count=barcode_count, building_block_csvs=building_block_csvs)
    library_json.write_text(json.dumps(library, indent=2) + "\n")

    decode_config = make_decode_config(
        barcode_count=barcode_count,
        fastq_path=fastq_path,
        library_json=library_json,
        min_read_length=min_read_length,
        num_reads_per_compound=num_reads_per_compound,
    )
    (experiment_dir / "decode_synthetic.yaml").write_text(yaml.safe_dump(decode_config, sort_keys=False))
    write_decode_settings_v02(experiment_dir / "decode_settings_v02.yaml", min_read_length=min_read_length)
    write_deli_config(experiment_dir / "deli_config", experiment_dir / "data")
    return read_count


def main(
    *,
    num_reads_per_compound: int = DEFAULT_NUM_READS_PER_COMPOUND,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
) -> None:
    if num_reads_per_compound < 1:
        raise ValueError("num_reads_per_compound must be at least 1")

    building_blocks_by_cycle = {
        cycle: make_building_blocks(f"BB{cycle}_", count, offset=(cycle - 1) * 10_000)
        for cycle, count in BUILDING_BLOCK_COUNTS_BY_CYCLE.items()
    }

    generated_read_counts = {}
    output_root.mkdir(parents=True, exist_ok=True)
    for barcode_count in BARCODE_COUNTS:
        generated_read_counts[barcode_count] = write_experiment(
            output_root=output_root,
            barcode_count=barcode_count,
            building_blocks_by_cycle=building_blocks_by_cycle,
            num_reads_per_compound=num_reads_per_compound,
        )

    for cycle, rows in building_blocks_by_cycle.items():
        print(f"Prepared {len(rows)} cycle-{cycle} building blocks")
    for barcode_count in BARCODE_COUNTS:
        experiment_name = f"synthetic_{barcode_count}cycle"
        compounds = compound_count([building_blocks_by_cycle[cycle] for cycle in range(1, barcode_count + 1)])
        print(
            f"Wrote {generated_read_counts[barcode_count]} FASTQ reads "
            f"({compounds} compounds x {num_reads_per_compound} reads) "
            f"to {output_root / experiment_name / 'input' / f'{experiment_name}.fastq'}"
        )
    print(f"Wrote DELi experiments to {output_root}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(
        num_reads_per_compound=cli_args.num_reads_per_compound,
        output_root=cli_args.output_root,
    )
