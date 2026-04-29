from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path
import sys


def load_generator_module():
    script_path = Path(__file__).resolve().parents[1] / "benchmarks" / "demultiplex" / "generate_synthetic_fastq.py"
    spec = importlib.util.spec_from_file_location("generate_synthetic_fastq", script_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def count_fastq_reads(path: Path) -> int:
    open_func = gzip.open if path.suffix == ".gz" else open
    with open_func(path, "rt") as handle:
        return sum(1 for _ in handle) // 4


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_small_smoke_case(tmp_path):
    module = load_generator_module()

    module.main(
        num_cycles=2,
        building_blocks_per_cycle=3,
        num_reads_per_compound=2,
        output_dir=tmp_path,
        experiment_name="smoke",
    )

    experiment_dir = tmp_path / "smoke_err=0"
    manifest = json.loads((experiment_dir / "manifest.json").read_text())
    expected_counts = read_tsv(experiment_dir / "expected_counts.tsv")

    assert manifest["expected_compounds"] == 9
    assert manifest["expected_reads"] == 18
    assert manifest["compressed_fastq"] is True
    assert count_fastq_reads(experiment_dir / "smoke_err=0.fastq.gz") == 18
    assert len(expected_counts) == 9
    assert {int(row["count"]) for row in expected_counts} == {2}
    assert {int(row["code_1"]) for row in expected_counts} == {1, 2, 3}
    assert {int(row["code_2"]) for row in expected_counts} == {1, 2, 3}


def test_higher_cycle_schema_and_building_blocks(tmp_path):
    module = load_generator_module()

    module.main(
        num_cycles=4,
        building_blocks_per_cycle=2,
        num_reads_per_compound=1,
        output_dir=tmp_path,
        experiment_name="shape",
    )

    experiment_dir = tmp_path / "shape_err=0"
    manifest = json.loads((experiment_dir / "manifest.json").read_text())
    expected_counts = read_tsv(experiment_dir / "expected_counts.tsv")
    building_blocks = read_tsv(experiment_dir / "building_blocks.tsv")

    assert manifest["num_cycles"] == 4
    assert manifest["building_blocks_per_cycle"] == 2
    assert manifest["expected_compounds"] == 16
    assert manifest["expected_reads"] == 16
    assert manifest["compressed_fastq"] is True
    assert count_fastq_reads(experiment_dir / "shape_err=0.fastq.gz") == 16
    assert len(expected_counts) == 16
    assert len(building_blocks) == 8
    assert set(expected_counts[0]) == {
        "code_1",
        "code_2",
        "code_3",
        "code_4",
        "building_block_id_1",
        "building_block_id_2",
        "building_block_id_3",
        "building_block_id_4",
        "count",
    }
