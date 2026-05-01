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
        dataset_name="smoke",
    )

    experiment_dir = tmp_path / "smoke"
    manifest = json.loads((experiment_dir / "manifest.json").read_text())
    expected_counts = read_tsv(experiment_dir / "expected_counts.tsv")

    assert manifest["expected_compounds"] == 9
    assert manifest["expected_reads"] == 18
    assert manifest["compressed_fastq"] is True
    assert count_fastq_reads(experiment_dir / "smoke.fastq.gz") == 18
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
        dataset_name="shape",
    )

    experiment_dir = tmp_path / "shape"
    manifest = json.loads((experiment_dir / "manifest.json").read_text())
    expected_counts = read_tsv(experiment_dir / "expected_counts.tsv")
    building_blocks = read_tsv(experiment_dir / "building_blocks.tsv")

    assert manifest["num_cycles"] == 4
    assert manifest["building_blocks_per_cycle"] == 2
    assert manifest["expected_compounds"] == 16
    assert manifest["expected_reads"] == 16
    assert manifest["compressed_fastq"] is True
    assert count_fastq_reads(experiment_dir / "shape.fastq.gz") == 16
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


def test_building_block_generation_only_enforces_hamming_distance_when_errors_enabled():
    module = load_generator_module()

    no_error_blocks = module.make_building_blocks(
        num_cycles=1,
        building_blocks_per_cycle=4,
        num_errors=0,
    )
    error_tolerant_blocks = module.make_building_blocks(
        num_cycles=1,
        building_blocks_per_cycle=4,
        num_errors=1,
    )

    no_error_tags = [row["tag"] for row in no_error_blocks[0]]
    error_tolerant_tags = [row["tag"] for row in error_tolerant_blocks[0]]

    assert no_error_tags == [module.int_to_dna(idx, module.TAG_LENGTH) for idx in range(4)]
    assert module.hamming_distance(no_error_tags[0], no_error_tags[1]) == 1
    assert all(
        module.hamming_distance(left, right) >= module.DEFAULT_BARCODE_MIN_HAMMING_DISTANCE
        for idx, left in enumerate(error_tolerant_tags)
        for right in error_tolerant_tags[idx + 1 :]
    )
