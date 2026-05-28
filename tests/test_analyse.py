from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from delt_hit.cli.analyse.api import Analyse
from delt_hit.cli.analyse.zscore import compute_zscore_stats


def write_yaml(path: Path, data: dict) -> Path:
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def write_counts(path: Path, rows: list[dict]) -> Path:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def make_analysis_config(tmp_path: Path) -> tuple[Path, Path]:
    selection_a = write_counts(
        tmp_path / "protein_counts.txt",
        [
            {"code_0": 0, "code_1": 0, "count": 10, "id": "0_0"},
            {"code_0": 0, "code_1": 1, "count": 5, "id": "0_1"},
        ],
    )
    selection_b = write_counts(
        tmp_path / "no_protein_counts.txt",
        [
            {"code_0": 0, "code_1": 0, "count": 3, "id": "0_0"},
            {"code_0": 0, "code_1": 1, "count": 4, "id": "0_1"},
        ],
    )
    analysis_config = write_yaml(
        tmp_path / "analysis.yaml",
        {
            "experiments": [
                {
                    "name": "protein_vs_no_protein",
                    "save_dir": str(tmp_path / "analysis_output"),
                    "selections": [
                        {"name": "protein_rep1", "group": "protein", "counts_path": str(selection_a)},
                        {"name": "no_protein_rep1", "group": "no_protein", "counts_path": str(selection_b)},
                    ],
                }
            ]
        },
    )
    return analysis_config, tmp_path / "analysis_output" / "protein_vs_no_protein"


def make_project_config(tmp_path: Path) -> tuple[Path, Path]:
    config_path = write_yaml(
        tmp_path / "config.yaml",
        {
            "experiment": {
                "name": "demo_project",
                "save_dir": str(tmp_path / "project_output"),
            },
            "library": {
                "building_blocks": ["B0", "B1"],
            },
            "whitelists": {
                "B0": [{"codon": "AAA"}, {"codon": "CCC"}],
                "B1": [{"codon": "GGG"}, {"codon": "TTT"}, {"codon": "ATA"}],
            },
        },
    )
    counts_path = write_counts(
        tmp_path / "counts.txt",
        [
            {"code_0": 0, "code_1": 0, "count": 6, "id": "0_0"},
            {"code_0": 0, "code_1": 1, "count": 1, "id": "0_1"},
            {"code_0": 1, "code_1": 0, "count": 0, "id": "1_0"},
        ],
    )
    return config_path, counts_path


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"method": "counts", "name": "x"}, "`analysis_config` is required"),
        ({"method": "counts", "analysis_config": Path("analysis.yaml")}, "`name` is required"),
        ({"method": "counts", "analysis_config": Path("analysis.yaml"), "name": "x", "counts": Path("counts.txt")}, "`counts` is not supported"),
        ({"method": "edgeR", "analysis_config": Path("analysis.yaml"), "name": "x", "config_path": Path("config.yaml")}, "`config_path` is not supported"),
        ({"method": "z_score", "counts": Path("counts.txt")}, "`config_path` is required"),
        ({"method": "z_score", "config_path": Path("config.yaml")}, "`counts` is required"),
        ({"method": "z_score", "config_path": Path("config.yaml"), "counts": Path("counts.txt"), "analysis_config": Path("analysis.yaml")}, "`analysis_config` is not supported"),
        ({"method": "z_score", "config_path": Path("config.yaml"), "counts": Path("counts.txt"), "name": "x"}, "`name` is not supported"),
    ],
)
def test_enrichment_argument_validation(kwargs, message):
    analyse = Analyse()
    with pytest.raises(ValueError, match=message):
        analyse.enrichment(**kwargs)


def test_deseq2_raises_not_implemented(tmp_path):
    analysis_config, _ = make_analysis_config(tmp_path)
    analyse = Analyse()

    with pytest.raises(NotImplementedError, match="DESeq2 analysis is not implemented"):
        analyse.enrichment(
            method="DESeq2",
            analysis_config=analysis_config,
            name="protein_vs_no_protein",
        )


def test_counts_analysis_still_generates_merged_inputs_and_script(tmp_path):
    analysis_config, output_root = make_analysis_config(tmp_path)
    analyse = Analyse()

    analyse.enrichment(
        method="counts",
        analysis_config=analysis_config,
        name="protein_vs_no_protein",
    )

    data = pd.read_csv(output_root / "data.csv")
    samples = pd.read_csv(output_root / "samples.csv")
    script_path = output_root / "counts" / "enrichment_counts.R"

    assert list(samples.columns) == ["name", "group"]
    assert set(data["name"]) == {"protein_rep1", "no_protein_rep1"}
    assert script_path.exists()


def test_edger_analysis_generates_script(tmp_path):
    analysis_config, output_root = make_analysis_config(tmp_path)
    analyse = Analyse()

    analyse.enrichment(
        method="edgeR",
        analysis_config=analysis_config,
        name="protein_vs_no_protein",
    )

    assert (output_root / "edgeR" / "enrichment_edgeR.R").exists()


def test_zscore_analysis_writes_stats_hits_and_script(tmp_path):
    config_path, counts_path = make_project_config(tmp_path)
    analyse = Analyse()

    analyse.enrichment(
        method="z_score",
        config_path=config_path,
        counts=counts_path,
    )

    output_dir = tmp_path / "project_output" / "demo_project" / "z_score"
    stats = pd.read_csv(output_dir / "stats.csv")
    hits = pd.read_csv(output_dir / "hits.csv")

    assert (output_dir / "enrichment_z_score.R").exists()
    assert "expected_count" in stats.columns
    assert "sigma" in stats.columns
    assert "z_score" in stats.columns
    assert "norm_z_score" in stats.columns
    assert hits.iloc[0]["norm_z_score"] == pytest.approx(stats["norm_z_score"].max())


def test_compute_zscore_stats_matches_expected_values():
    counts = pd.DataFrame(
        [
            {"code_0": 0, "code_1": 0, "count": 6},
            {"code_0": 0, "code_1": 1, "count": 1},
            {"code_0": 1, "code_1": 0, "count": 0},
        ]
    )
    stats = compute_zscore_stats(counts=counts, library_size=6)

    total_count = 7.0
    expected_count = total_count / 6.0
    sigma = (total_count * (1.0 / 6.0) * (5.0 / 6.0)) ** 0.5

    assert stats["expected_count"].tolist() == pytest.approx([expected_count] * len(stats))
    assert stats["sigma"].tolist() == pytest.approx([sigma] * len(stats))
    assert stats.loc[0, "z_score"] == pytest.approx((6.0 - expected_count) / sigma)
    assert stats.loc[0, "norm_z_score"] == pytest.approx(((6.0 - expected_count) / sigma) / (total_count**0.5))


def test_zscore_requires_derivable_library_size(tmp_path):
    bad_config = write_yaml(
        tmp_path / "bad_config.yaml",
        {
            "experiment": {
                "name": "demo_project",
                "save_dir": str(tmp_path / "project_output"),
            },
            "library": {"building_blocks": ["B0"]},
            "whitelists": {},
        },
    )
    counts_path = write_counts(tmp_path / "counts.txt", [{"code_0": 0, "count": 1}])
    analyse = Analyse()

    with pytest.raises(ValueError, match="missing whitelist entries"):
        analyse.enrichment(method="z_score", config_path=bad_config, counts=counts_path)


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([{"count": 1}], "at least one `code_\\*` column"),
        ([{"code_0": 0}], "must contain a `count` column"),
        ([{"code_0": 0, "count": 0}], "positive total count"),
    ],
)
def test_zscore_validates_counts_schema(tmp_path, rows, message):
    config_path = write_yaml(
        tmp_path / "config.yaml",
        {
            "experiment": {
                "name": "demo_project",
                "save_dir": str(tmp_path / "project_output"),
            },
            "library": {"building_blocks": ["B0"]},
            "whitelists": {"B0": [{"codon": "AAA"}]},
        },
    )
    counts_path = write_counts(tmp_path / "counts.txt", rows)
    analyse = Analyse()

    with pytest.raises(ValueError, match=message):
        analyse.enrichment(method="z_score", config_path=config_path, counts=counts_path)
