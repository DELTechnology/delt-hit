#!/usr/bin/env python3
"""Convert a neutral synthetic dataset into a runnable DELi benchmark sandbox."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[3]
BENCHMARKS_ROOT = PROJECT_ROOT / "benchmarks" / "demultiplex"
DATA_ROOT = BENCHMARKS_ROOT / "data"
TOOLS_ROOT = BENCHMARKS_ROOT / "tools"
TOOL_ROOT = TOOLS_ROOT / "deli"
DEFAULT_DATASET_NAME = "synthetic_2cycle_10bbpc_1m"


def read_manifest(dataset_dir: Path) -> dict:
    return json.loads((dataset_dir / "manifest.json").read_text())


def read_building_blocks(dataset_dir: Path) -> list[dict[str, str]]:
    with (dataset_dir / "building_blocks.tsv").open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_deli_config(path: Path, data_dir: Path) -> None:
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-name",
        default=DEFAULT_DATASET_NAME,
        help="Dataset directory name under benchmarks/demultiplex/data/.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Scratch workspace root for one dataset. Uses <data-dir>/data/<dataset>, caches benchmark config under benchmarks/demultiplex/tools/deli/<dataset>, and writes generated DELi input assets under <data-dir>/tools/deli/<dataset>.",
    )
    return parser.parse_args()


def main(*, dataset_name: str = DEFAULT_DATASET_NAME, data_dir: Path | None = None) -> None:
    config_dir = (TOOL_ROOT / dataset_name).resolve()
    if data_dir is None:
        dataset_dir = (DATA_ROOT / dataset_name).resolve()
        runtime_dir = config_dir
    else:
        data_dir = data_dir.resolve()
        dataset_dir = data_dir / "data" / dataset_name
        runtime_dir = data_dir / "tools" / "deli" / dataset_name

    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset directory not found: {dataset_dir}")

    manifest = read_manifest(dataset_dir)
    building_blocks = read_building_blocks(dataset_dir)
    dataset_name = manifest.get("dataset_name", manifest.get("experiment_name"))
    if dataset_name is None:
        raise KeyError("Manifest must contain dataset_name")
    num_cycles = manifest["num_cycles"]
    num_errors = int(manifest.get("num_errors", 0))
    if num_errors not in {0, 1}:
        raise ValueError(f"Unsupported num_errors={num_errors}; only 0 and 1 are supported")

    data_dir = runtime_dir / "data"
    bb_dir = data_dir / "building_blocks"
    lib_dir = data_dir / "libraries"
    reaction_dir = data_dir / "reactions"
    tool_dir = data_dir / "tool_compounds"

    for path in [bb_dir, lib_dir, reaction_dir, tool_dir]:
        path.mkdir(parents=True, exist_ok=True)

    for cycle in range(1, num_cycles + 1):
        cycle_rows = [
            {"id": row["building_block_id"], "tag": row["tag"]}
            for row in building_blocks
            if int(row["cycle"]) == cycle
        ]
        with (bb_dir / f"{dataset_name}_bb{cycle}.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["id", "tag"])
            writer.writeheader()
            writer.writerows(cycle_rows)

    library_json = lib_dir / f"{dataset_name}_library.json"
    library_data = {
        "dna_barcode_on": None,
        "bb_sets": [
            {
                "cycle": cycle,
                "bb_set_path": str((bb_dir / f"{dataset_name}_bb{cycle}.csv").resolve()),
            }
            for cycle in range(1, num_cycles + 1)
        ],
        "barcode_schema": {
            "umi": {"tag": "N" * manifest["umi_length"]},
            "library": {"tag": manifest["library_tag"]},
            **{
                f"bb{cycle}": {
                    "tag": "N" * manifest["tag_length"],
                    **({"error_correction": "hamming_dist:1"} if num_errors == 1 else {}),
                }
                for cycle in range(1, num_cycles + 1)
            },
            "closing": {"tag": manifest["closing_tag"]},
        },
    }
    library_json.write_text(json.dumps(library_data, indent=2) + "\n")

    config_dir.mkdir(parents=True, exist_ok=True)

    selection_yaml = config_dir / "decode_synthetic.yaml"
    selection_data = {
        "selection_id": f"{dataset_name}_selection",
        "target_id": "synthetic_target",
        "selection_condition": "synthetic_benchmark",
        "additional_info": (
            f"Synthetic neutral dataset with {manifest['expected_reads']} reads "
            f"across {manifest['expected_compounds']} encoded compounds."
        ),
        "data_ran": "2026-04-24T00:00:00",
        "sequence_files": [str(Path(manifest["files"]["fastq"]).resolve())],
        "libraries": [str(library_json.resolve())],
        "decode_settings": {
            "library_error_tolerance": num_errors,
            "min_library_overlap": len(manifest["library_tag"]),
            "alignment_algorithm": "semi",
            "bb_calling_approach": "bio",
            "revcomp": False,
            "max_read_length": 128,
            "min_read_length": (
                len(manifest["library_tag"]) + num_cycles * manifest["tag_length"] + manifest["umi_length"]
            ),
            "read_type": "single",
            "disable_error_correction": num_errors == 0,
            "umi_clustering": False,
            "umi_min_distance": 2,
        },
    }
    selection_yaml.write_text(yaml.safe_dump(selection_data, sort_keys=False))

    decode_settings_path = config_dir / "decode_settings_v02.yaml"
    decode_settings_path.write_text(
        yaml.safe_dump(
            {
                "ignore_tool_compounds": True,
                "demultiplexer_algorithm": "regex",
                "demultiplexer_mode": "library",
                "realign": False,
                "library_error_tolerance": num_errors,
                "library_error_correction_mode_str": "hamming_dist:1" if num_errors == 1 else "disable",
                "min_library_overlap": len(manifest["library_tag"]),
                "revcomp": False,
                "library_wiggle": False,
                "wiggle": False,
                "decode_matching_approach": "first_perfect",
                "max_read_length": 128,
                "min_read_length": (
                    len(manifest["library_tag"]) + num_cycles * manifest["tag_length"] + manifest["umi_length"]
                ),
                "default_error_correction_mode_str": "hamming_dist:1" if num_errors == 1 else "disable",
            },
            sort_keys=False,
        )
    )

    deli_config = config_dir / "deli_config"
    write_deli_config(deli_config, data_dir)

    print(f"Wrote DELi selection file to {selection_yaml}")
    print(f"Wrote DELi decode settings to {decode_settings_path}")
    print(f"Wrote DELi config to {deli_config}")
    print(f"Wrote DELi library JSON to {library_json}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(dataset_name=cli_args.dataset_name, data_dir=cli_args.data_dir)
