#!/usr/bin/env python3
"""Convert neutral synthetic FASTQ artifacts into a DELT-Hit config.yaml.

This converter targets demultiplex/count benchmarking only. It emits the
whitelist/config structure DELT-Hit needs, but does not add chemical
enumeration metadata.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from delt_hit.utils import write_yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATASET_DIR = PROJECT_ROOT / "benchmarks" / "data" / "synthetic_2cycle_1m"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "benchmarks" / "tools" / "delt"


def read_manifest(dataset_dir: Path) -> dict:
    return json.loads((dataset_dir / "manifest.json").read_text())


def read_building_blocks(dataset_dir: Path) -> list[dict[str, str]]:
    with (dataset_dir / "building_blocks.tsv").open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=DEFAULT_DATASET_DIR,
        help="Directory containing manifest.json, building_blocks.tsv, and the FASTQ file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where DELT-Hit inputs are written. Defaults to benchmarks/tools/delt/<dataset-name>.",
    )
    parser.add_argument(
        "--num-cores",
        type=int,
        default=1,
        help="Value to write into experiment.num_cores in config.yaml.",
    )
    return parser.parse_args()

def main(
    *,
    dataset_dir: Path = DEFAULT_DATASET_DIR,
    output_dir: Path | None = None,
    num_cores: int = 1,
) -> None:
    dataset_dir = dataset_dir.expanduser().resolve()
    manifest = read_manifest(dataset_dir)
    building_blocks = read_building_blocks(dataset_dir)
    experiment_name = manifest["experiment_name"]
    num_cycles = manifest["num_cycles"]

    output_dir = (output_dir or (DEFAULT_OUTPUT_DIR / dataset_dir.name)).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    fastq_path = Path(manifest["files"]["fastq"]).resolve()

    whitelists = {
        "C0": [{"codon": manifest["closing_tag"]}],
        "S0": [{"name": "synthetic_selection", "codon": manifest["library_tag"]}],
        "S1": [{"name": "synthetic_selection", "codon": "N" * manifest["umi_length"] + manifest["closing_tag"]}],
    }
    for cycle in range(1, num_cycles + 1):
        sheet_name = f"B{cycle - 1}"
        cycle_rows = [
            {
                "index": int(row["code"]) - 1,
                "codon": row["tag"],
                "building_block_id": row["building_block_id"],
            }
            for row in building_blocks
            if int(row["cycle"]) == cycle
        ]
        whitelists[sheet_name] = cycle_rows

    config = {
        "experiment": {
            "name": experiment_name,
            "fastq_path": str(fastq_path),
            "save_dir": str(output_dir.resolve()),
            "num_cores": num_cores,
        },
        "selections": {
            "synthetic_selection": {
                "operator": "Codex",
                "date": "2026-04-24",
                "group": "synthetic",
                "target": "synthetic_target",
                "S0": manifest["library_tag"],
                "S1": "N" * manifest["umi_length"] + manifest["closing_tag"],
                "ids": [0, 0],
            }
        },
        "structure": [
            {"name": "S0", "type": "selection", "max_error_rate": 0, "indels": False},
            *[
                {"name": f"B{cycle - 1}", "type": "building_block", "max_error_rate": 0, "indels": False}
                for cycle in range(1, num_cycles + 1)
            ],
            {"name": "S1", "type": "selection", "max_error_rate": 0, "indels": False},
        ],
        "whitelists": whitelists,
        "library": {
            "products": [],
            "educts": [],
            "reactions": [],
            "bb_edges": [],
            "other_edges": [],
            "building_blocks": [f"B{cycle - 1}" for cycle in range(1, num_cycles + 1)],
        },
        "catalog": {
            "compounds": {},
            "reactions": {},
        },
    }

    config_path = output_dir / "config.yaml"
    write_yaml(config, config_path)

    print(f"Wrote DELT-Hit config to {config_path}")
    print(f"Referenced FASTQ at {fastq_path}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(
        dataset_dir=cli_args.dataset_dir,
        output_dir=cli_args.output_dir,
        num_cores=cli_args.num_cores,
    )
