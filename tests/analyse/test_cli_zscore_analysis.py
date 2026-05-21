from pathlib import Path

import pandas as pd
import pytest

from delt_hit.cli.analyse.api import Analyse
from delt_hit.cli.analyse.zscore import compute_zscore, zscore_analysis


def _write_counts(path: Path, rows: list[tuple[str, str, int]]) -> None:
    frame = pd.DataFrame(rows, columns=["code_1", "code_2", "count"])
    frame["id"] = [f"cmpd_{idx}" for idx in range(len(frame))]
    frame.to_csv(path, sep="\t", index=False)


def test_zscore_analysis_writes_protocol_outputs_and_zero_fills(tmp_path: Path):
    data = pd.DataFrame(
        [
            {"code_1": "A", "code_2": "A", "name": "protein_1", "count": 90},
            {"code_1": "A", "code_2": "B", "name": "protein_1", "count": 10},
            {"code_1": "A", "code_2": "A", "name": "control_1", "count": 10},
            {"code_1": "A", "code_2": "B", "name": "control_1", "count": 70},
            {"code_1": "A", "code_2": "C", "name": "control_1", "count": 20},
        ]
    )
    samples = pd.DataFrame(
        [
            {"name": "protein_1", "group": "protein"},
            {"name": "control_1", "group": "no_protein"},
        ]
    )

    save_dir = tmp_path / "zscore"
    zscore_analysis(data=data, samples=samples, library_size=4, save_dir=save_dir)

    assert (save_dir / "stats.csv").exists()
    assert (save_dir / "hits.csv").exists()
    assert (save_dir / "protein_1.csv").exists()
    assert (save_dir / "control_1.csv").exists()

    protein = pd.read_csv(save_dir / "protein_1.csv")
    zero_filled = protein.loc[(protein["code_1"] == "A") & (protein["code_2"] == "C")]
    assert zero_filled["count"].item() == 0
    assert protein.loc[protein["code_2"] == "A", "expected_fraction"].item() == 0.25
    assert protein.loc[protein["code_2"] == "A", "expected_count"].item() == 25.0
    assert protein.loc[protein["code_2"] == "A", "fold_enrichment"].item() == 3.6

    stats = pd.read_csv(save_dir / "stats.csv")
    assert list(stats.columns) == ["code_1", "code_2", "protein", "no_protein", "enrichment"]
    top_hit = pd.read_csv(save_dir / "hits.csv").iloc[0]
    assert top_hit["code_1"] == "A"
    assert top_hit["code_2"] == "A"
    assert top_hit["enrichment"] > 0


def test_zscore_analysis_supports_three_cycle_code_columns(tmp_path: Path):
    data = pd.DataFrame(
        [
            {"code_1": "A", "code_2": "B", "code_3": "C", "name": "sel_1", "count": 8},
            {"code_1": "A", "code_2": "B", "code_3": "D", "name": "sel_1", "count": 2},
            {"code_1": "A", "code_2": "B", "code_3": "C", "name": "sel_2", "count": 1},
            {"code_1": "A", "code_2": "B", "code_3": "D", "name": "sel_2", "count": 9},
        ]
    )
    samples = pd.DataFrame(
        [
            {"name": "sel_1", "group": "protein"},
            {"name": "sel_2", "group": "no_protein"},
        ]
    )

    zscore_analysis(data=data, samples=samples, library_size=8, save_dir=tmp_path)

    stats = pd.read_csv(tmp_path / "stats.csv")
    assert list(stats.columns[:3]) == ["code_1", "code_2", "code_3"]
    assert len(stats) == 2


def test_compute_zscore_rejects_zero_total_selection():
    data = pd.DataFrame(
        [
            {"code_1": "A", "name": "empty", "group": "protein", "count": 0},
            {"code_1": "B", "name": "empty", "group": "protein", "count": 0},
        ]
    )

    with pytest.raises(ValueError, match="zero total"):
        compute_zscore(data=data, library_size=2)


def test_zscore_analysis_rejects_missing_code_columns(tmp_path: Path):
    data = pd.DataFrame([{"name": "sel_1", "count": 10}])
    samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

    with pytest.raises(ValueError, match="code_\\*"):
        zscore_analysis(data=data, samples=samples, library_size=10, save_dir=tmp_path)


def test_zscore_analysis_rejects_missing_code_values(tmp_path: Path):
    data = pd.DataFrame(
        [
            {"code_1": "A", "name": "sel_1", "count": 10},
            {"code_1": None, "name": "sel_1", "count": 5},
        ]
    )
    samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

    with pytest.raises(ValueError, match="missing code_\\* values"):
        zscore_analysis(data=data, samples=samples, library_size=10, save_dir=tmp_path)


def test_zscore_analysis_rejects_library_size_smaller_than_observed(tmp_path: Path):
    data = pd.DataFrame(
        [
            {"code_1": "A", "name": "sel_1", "count": 10},
            {"code_1": "B", "name": "sel_1", "count": 5},
            {"code_1": "C", "name": "sel_1", "count": 1},
        ]
    )
    samples = pd.DataFrame([{"name": "sel_1", "group": "protein"}])

    with pytest.raises(ValueError, match="smaller than the observed"):
        zscore_analysis(data=data, samples=samples, library_size=2, save_dir=tmp_path)


def test_analyse_enrichment_zscore_wires_cli_method(tmp_path: Path):
    protein_counts = tmp_path / "protein.tsv"
    control_counts = tmp_path / "control.tsv"
    _write_counts(protein_counts, [("A", "A", 90), ("A", "B", 10)])
    _write_counts(control_counts, [("A", "A", 10), ("A", "B", 90)])

    config = tmp_path / "analysis.yaml"
    output_root = tmp_path / "out"
    config.write_text(
        f"""
experiments:
  - name: demo
    save_dir: {output_root}
    library_size: 4
    selections:
      - name: protein_1
        counts_path: {protein_counts}
        group: protein
      - name: control_1
        counts_path: {control_counts}
        group: no_protein
""".strip()
    )

    Analyse().enrichment(config_path=config, name="demo", method="zscore")

    experiment_dir = output_root / "demo"
    assert (experiment_dir / "data.csv").exists()
    assert (experiment_dir / "samples.csv").exists()
    assert (experiment_dir / "zscore" / "stats.csv").exists()
    assert (experiment_dir / "zscore" / "hits.csv").exists()


def test_analyse_enrichment_rejects_unknown_method(tmp_path: Path):
    counts = tmp_path / "counts.tsv"
    _write_counts(counts, [("A", "A", 10), ("A", "B", 5)])
    config = tmp_path / "analysis.yaml"
    config.write_text(
        f"""
experiments:
  - name: demo
    save_dir: {tmp_path / "out"}
    selections:
      - name: protein_1
        counts_path: {counts}
        group: protein
""".strip()
    )

    with pytest.raises(ValueError, match="Unknown enrichment method"):
        Analyse().enrichment(config_path=config, name="demo", method="unknown")
