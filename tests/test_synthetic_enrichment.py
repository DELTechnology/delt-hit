from __future__ import annotations

import importlib.util
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest


def load_module(module_name: str, relative_path: str):
    module_path = Path(__file__).resolve().parents[1] / relative_path
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_generator_module():
    return load_module("generate_synthetic_enrichment", "benchmarks/generate_synthetic_enrichment.py")


def load_compare_module():
    return load_module("compare_synthetic_enrichment", "benchmarks/compare_synthetic_enrichment.py")


def r_with_required_packages_available() -> bool:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False
    cmd = [
        rscript,
        "--vanilla",
        "-e",
        (
            "packages <- c('tidyverse', 'GGally', 'edgeR', 'limma'); "
            "status <- if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) 0 else 1; "
            "quit(status = status)"
        ),
    ]
    return subprocess.run(cmd, capture_output=True, check=False).returncode == 0


def test_generator_smoke_case(tmp_path):
    module = load_generator_module()

    dataset_dir = module.main(
        building_blocks_per_cycle=4,
        positive_low=3,
        positive_medium=2,
        positive_high=1,
        output_dir=tmp_path,
        experiment_name="smoke",
        seed=7,
    )

    assert dataset_dir == tmp_path / "smoke"
    assert (dataset_dir / "config.yaml").exists()
    assert (dataset_dir / "truth.tsv").exists()
    assert (dataset_dir / "manifest.json").exists()
    for selection_name in module.selection_names():
        assert (dataset_dir / "selections" / selection_name / "counts.txt").exists()


def test_generated_count_tables_match_analysis_schema(tmp_path):
    module = load_generator_module()
    dataset_dir = module.main(
        building_blocks_per_cycle=5,
        positive_low=2,
        positive_medium=2,
        positive_high=1,
        output_dir=tmp_path,
        experiment_name="schema",
        seed=11,
    )

    for selection_name in module.selection_names():
        table = pd.read_csv(dataset_dir / "selections" / selection_name / "counts.txt", sep="\t")
        assert list(table.columns) == ["code_0", "code_1", "count", "id"]
        assert table["code_0"].min() == 0
        assert table["code_1"].min() == 0
        assert table["code_0"].max() == 4
        assert table["code_1"].max() == 4


def test_truth_panel_counts_are_exact_and_non_overlapping(tmp_path):
    module = load_generator_module()
    dataset_dir = module.main(
        building_blocks_per_cycle=6,
        positive_low=4,
        positive_medium=3,
        positive_high=2,
        output_dir=tmp_path,
        experiment_name="truth",
        seed=5,
    )

    metadata, truth = module.read_truth_table(dataset_dir / "truth.tsv")

    assert metadata["positive_low"] == 4
    assert metadata["positive_medium"] == 3
    assert metadata["positive_high"] == 2
    assert (truth["truth_tier"] == "positive_low").sum() == 4
    assert (truth["truth_tier"] == "positive_medium").sum() == 3
    assert (truth["truth_tier"] == "positive_high").sum() == 2
    assert (truth["truth_tier"] == "negative").sum() == 36 - 9
    assert truth["id"].is_unique


def test_negatives_have_no_built_in_group_shift(tmp_path):
    module = load_generator_module()
    dataset_dir = module.main(
        building_blocks_per_cycle=6,
        positive_low=3,
        positive_medium=2,
        positive_high=1,
        output_dir=tmp_path,
        experiment_name="nulls",
        seed=17,
    )

    _, truth = module.read_truth_table(dataset_dir / "truth.tsv")
    negatives = truth[truth["truth_tier"] == "negative"]

    assert (negatives["target_fold_change"] == 1.0).all()
    assert negatives["protein_mean"].equals(negatives["no_protein_mean"])


def test_same_seed_is_reproducible(tmp_path):
    module = load_generator_module()
    first_dir = module.main(
        building_blocks_per_cycle=8,
        positive_low=5,
        positive_medium=3,
        positive_high=2,
        output_dir=tmp_path,
        experiment_name="seed_a",
        seed=23,
    )
    second_dir = module.main(
        building_blocks_per_cycle=8,
        positive_low=5,
        positive_medium=3,
        positive_high=2,
        output_dir=tmp_path,
        experiment_name="seed_b",
        seed=23,
    )

    _, truth_first = module.read_truth_table(first_dir / "truth.tsv")
    _, truth_second = module.read_truth_table(second_dir / "truth.tsv")
    pd.testing.assert_frame_equal(truth_first, truth_second)

    for selection_name in module.selection_names():
        counts_first = pd.read_csv(first_dir / "selections" / selection_name / "counts.txt", sep="\t")
        counts_second = pd.read_csv(second_dir / "selections" / selection_name / "counts.txt", sep="\t")
        pd.testing.assert_frame_equal(counts_first, counts_second)


def test_higher_dispersion_increases_replicate_variance(tmp_path):
    module = load_generator_module()
    low_noise_dir = module.main(
        building_blocks_per_cycle=20,
        positive_low=10,
        positive_medium=5,
        positive_high=3,
        dispersion=0.01,
        output_dir=tmp_path,
        experiment_name="low_noise",
        seed=31,
    )
    high_noise_dir = module.main(
        building_blocks_per_cycle=20,
        positive_low=10,
        positive_medium=5,
        positive_high=3,
        dispersion=0.60,
        output_dir=tmp_path,
        experiment_name="high_noise",
        seed=31,
    )

    _, low_truth = module.read_truth_table(low_noise_dir / "truth.tsv")
    _, high_truth = module.read_truth_table(high_noise_dir / "truth.tsv")
    pd.testing.assert_frame_equal(low_truth, high_truth)

    def average_replicate_variance(dataset_dir: Path) -> float:
        tables = []
        for selection_name in module.selection_names():
            table = pd.read_csv(dataset_dir / "selections" / selection_name / "counts.txt", sep="\t")
            tables.append(table.rename(columns={"count": selection_name})[["id", selection_name]])
        merged = tables[0]
        for table in tables[1:]:
            merged = merged.merge(table, on="id")
        return float(merged.drop(columns=["id"]).var(axis=1).mean())

    assert average_replicate_variance(high_noise_dir) > average_replicate_variance(low_noise_dir)


def test_runner_scoring_prefers_strong_positive_controls():
    compare_module = load_compare_module()
    truth = pd.DataFrame(
        [
            {"id": "0_0", "code_0": 0, "code_1": 0, "truth_tier": "positive_high"},
            {"id": "0_1", "code_0": 0, "code_1": 1, "truth_tier": "positive_medium"},
            {"id": "1_0", "code_0": 1, "code_1": 0, "truth_tier": "positive_low"},
            {"id": "1_1", "code_0": 1, "code_1": 1, "truth_tier": "negative"},
        ]
    )
    counts_results = pd.DataFrame(
        [
            {"code_0": 0, "code_1": 0, "counts_score": 120.0},
            {"code_0": 0, "code_1": 1, "counts_score": 60.0},
            {"code_0": 1, "code_1": 0, "counts_score": 20.0},
            {"code_0": 1, "code_1": 1, "counts_score": 1.0},
        ]
    )
    edger_results = pd.DataFrame(
        [
            {"code_0": 0, "code_1": 0, "edgeR_logFC": 4.0, "edgeR_FDR": 0.001},
            {"code_0": 0, "code_1": 1, "edgeR_logFC": 2.5, "edgeR_FDR": 0.010},
            {"code_0": 1, "code_1": 0, "edgeR_logFC": 1.2, "edgeR_FDR": 0.050},
            {"code_0": 1, "code_1": 1, "edgeR_logFC": 0.1, "edgeR_FDR": 0.900},
        ]
    )

    comparison = compare_module.build_comparison_table(truth, counts_results, edger_results)
    counts_metrics = compare_module.compute_method_metrics(
        comparison,
        score_column="counts_score",
        rank_column="counts_rank",
        top_k=3,
    )
    edger_metrics = compare_module.compute_method_metrics(
        comparison,
        score_column="edgeR_logFC",
        rank_column="edgeR_rank",
        top_k=3,
    )

    assert counts_metrics["precision_at_k"] == pytest.approx(1.0)
    assert edger_metrics["precision_at_k"] == pytest.approx(1.0)
    assert counts_metrics["tier_recall_at_k"]["positive_high"] == pytest.approx(1.0)
    assert edger_metrics["tier_recall_at_k"]["positive_high"] == pytest.approx(1.0)


@pytest.mark.skipif(not r_with_required_packages_available(), reason="R packages for counts/edgeR analysis are unavailable")
def test_generated_config_runs_through_counts_and_edger(tmp_path):
    generator = load_generator_module()
    compare_module = load_compare_module()

    dataset_dir = generator.main(
        building_blocks_per_cycle=5,
        positive_low=2,
        positive_medium=2,
        positive_high=1,
        dispersion=0.05,
        output_dir=tmp_path,
        experiment_name="r_smoke",
        seed=41,
    )

    output_dir = compare_module.main(dataset_dir=dataset_dir, top_k=5)

    assert (output_dir / "comparison.tsv").exists()
    assert (output_dir / "comparison.json").exists()
    assert (output_dir / "summary.md").exists()
