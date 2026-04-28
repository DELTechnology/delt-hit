#!/usr/bin/env python3
"""Convert a neutral synthetic dataset into a runnable DELT-Hit benchmark sandbox."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BENCHMARKS_ROOT = PROJECT_ROOT / "benchmarks"
DATA_ROOT = BENCHMARKS_ROOT / "data"
TOOLS_ROOT = BENCHMARKS_ROOT / "tools"
TOOL_ROOT = TOOLS_ROOT / "delt"
DEFAULT_DATASET_NAME = "synthetic_2cycle_1m"
DELT_HIT_PYTHON = PROJECT_ROOT / ".venv" / "bin" / "python"


def read_manifest(dataset_dir: Path) -> dict:
    return json.loads((dataset_dir / "manifest.json").read_text())


def read_building_blocks(dataset_dir: Path) -> list[dict[str, str]]:
    with (dataset_dir / "building_blocks.tsv").open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_command(cmd: list[str | Path], *, env: dict[str, str] | None = None) -> None:
    subprocess.run([str(item) for item in cmd], cwd=PROJECT_ROOT, env=env, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-name",
        default=DEFAULT_DATASET_NAME,
        help="Dataset directory name under benchmarks/data/.",
    )
    parser.add_argument(
        "--num-cores",
        type=int,
        default=1,
        help="Value to write into experiment.num_cores in config.yaml.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Scratch workspace root for one dataset. Uses <data-dir>/data/<dataset> and writes tool inputs under <data-dir>/tools/delt/<dataset>.",
    )
    return parser.parse_args()


def main(*, dataset_name: str = DEFAULT_DATASET_NAME, num_cores: int = 1, data_dir: Path | None = None) -> None:
    if data_dir is None:
        dataset_dir = (DATA_ROOT / dataset_name).resolve()
        output_dir = (TOOL_ROOT / dataset_name).resolve()
    else:
        data_dir = data_dir.resolve()
        dataset_dir = data_dir / "data" / dataset_name
        output_dir = data_dir / "tools" / "delt" / dataset_name

    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset directory not found: {dataset_dir}")

    manifest = read_manifest(dataset_dir)
    building_blocks = read_building_blocks(dataset_dir)
    experiment_name = manifest["experiment_name"]
    num_cycles = manifest["num_cycles"]
    num_errors = int(manifest.get("num_errors", 0))
    if num_errors not in {0, 1}:
        raise ValueError(f"Unsupported num_errors={num_errors}; only 0 and 1 are supported")

    output_dir.mkdir(parents=True, exist_ok=True)
    fastq_path = Path(manifest["files"]["fastq"]).resolve()
    library_error_rate = num_errors / len(manifest["library_tag"])
    barcode_error_rate = num_errors / int(manifest["tag_length"])
    umi_length = int(manifest["umi_length"])
    # Absorb UMI bases into C0 with Ns — cutadapt cannot demux all-N adapters,
    # so we lump the UMI skip into the closing-tag step.
    # Reads start with: [UMI][library_tag][bb0]...[bbN][closing_tag]
    # cutadapt anchors all adapters at 5', so S0 must absorb the UMI prefix
    # with Ns, then the barcodes strip naturally, and C0 absorbs any remaining
    # UMI bases (none, since S0 already consumed them).
    selection_codon = "N" * umi_length + manifest["library_tag"]
    selection_error_rate = num_errors / len(manifest["library_tag"])
    closing_codon = manifest["closing_tag"]
    closing_error_rate = num_errors / len(closing_codon)
    whitelists = {
        "C0": [{"codon": closing_codon}],
        "S0": [{"name": "synthetic_selection", "codon": selection_codon}],
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
                "ids": [0],
            }
        },
        "structure": [
            {"name": "S0", "type": "selection", "max_error_rate": selection_error_rate, "indels": False},
            *[
                {"name": f"B{cycle - 1}", "type": "building_block", "max_error_rate": barcode_error_rate, "indels": False}
                for cycle in range(1, num_cycles + 1)
            ],
            {"name": "C0", "type": "constant", "max_error_rate": closing_error_rate, "indels": False},
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
    with config_path.open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)

    env = os.environ.copy()
    env["PATH"] = f"{PROJECT_ROOT / '.venv' / 'bin'}:{env['PATH']}"
    run_command(
        [DELT_HIT_PYTHON, "-m", "delt_hit.cli.main", "demultiplex", "prepare", "--config_path", config_path],
        env=env,
    )

    print(f"Wrote DELT-Hit config to {config_path}")
    print(f"Prepared demultiplex inputs under {output_dir / experiment_name / 'demultiplex'}")


if __name__ == "__main__":
    cli_args = parse_args()
    main(dataset_name=cli_args.dataset_name, num_cores=cli_args.num_cores, data_dir=cli_args.data_dir)
