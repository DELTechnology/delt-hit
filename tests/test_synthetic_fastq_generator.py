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


def load_matrix_module():
    script_path = Path(__file__).resolve().parents[1] / "benchmarks" / "demultiplex" / "generate_synthetic_benchmark_matrix.py"
    spec = importlib.util.spec_from_file_location("generate_synthetic_benchmark_matrix", script_path)
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


def test_two_cycle_large_library_profile_has_expected_presets():
    module = load_matrix_module()

    specs = module.PROFILES["two_cycle_large_libraries"]
    by_name = {spec.experiment_name: spec for spec in specs}

    assert set(by_name) == {
        "synthetic_2cycle_100bbpc_1m",
        "synthetic_2cycle_100bbpc_10m",
        "synthetic_2cycle_100bbpc_100m",
        "synthetic_2cycle_100bbpc_1000m",
        "synthetic_2cycle_1000bbpc_1m",
        "synthetic_2cycle_1000bbpc_10m",
        "synthetic_2cycle_1000bbpc_100m",
        "synthetic_2cycle_1000bbpc_1000m",
    }
    assert by_name["synthetic_2cycle_100bbpc_1m"].num_reads_per_compound == 100
    assert by_name["synthetic_2cycle_100bbpc_1000m"].num_reads_per_compound == 100_000
    assert by_name["synthetic_2cycle_1000bbpc_1m"].num_reads_per_compound == 1
    assert by_name["synthetic_2cycle_1000bbpc_1000m"].num_reads_per_compound == 1_000


def test_matrix_profile_generation_writes_new_two_cycle_dataset(tmp_path):
    module = load_matrix_module()
    spec = module.make_dataset_spec(
        num_cycles=2,
        building_blocks_per_cycle=100,
        depth_label="1m",
    )

    module.generate_dataset(
        num_cycles=spec.num_cycles,
        building_blocks_per_cycle=spec.building_blocks_per_cycle,
        num_reads_per_compound=spec.num_reads_per_compound,
        output_dir=tmp_path,
        experiment_name=spec.experiment_name,
        compressed=True,
    )

    experiment_dir = tmp_path / f"{spec.experiment_name}_err=0"
    manifest = json.loads((experiment_dir / "manifest.json").read_text())

    assert manifest["experiment_name"] == "synthetic_2cycle_100bbpc_1m_err=0"
    assert manifest["building_blocks_per_cycle"] == 100
    assert manifest["num_cycles"] == 2
    assert manifest["num_reads_per_compound"] == 100
    assert manifest["expected_compounds"] == 10_000
    assert manifest["expected_reads"] == 1_000_000
