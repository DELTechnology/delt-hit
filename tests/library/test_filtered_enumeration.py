from pathlib import Path

import pandas as pd
import pytest

from delt_hit.cli.library.api import Library, get_combinations
from delt_hit.utils import write_yaml


def make_test_config(tmp_path: Path) -> tuple[dict, Path]:
    config = {
        "experiment": {
            "name": "mini",
            "save_dir": str(tmp_path),
        },
        "library": {
            "bb_edges": [["B0", "rxn_join"], ["B1", "rxn_join"], ["rxn_join", "product_1"]],
            "other_edges": [],
            "building_blocks": ["B0", "B1"],
            "products": ["product_1"],
        },
        "catalog": {
            "reactions": {
                "rxn_join": {"smirks": "[C:1].[O:2]>>[C:1][O:2]"},
            },
            "compounds": {},
        },
        "whitelists": {
            "B0": [
                {"index": 0, "smiles": "C", "reaction": "rxn_join", "product": "product_1", "educt": "B1"},
                {"index": 1, "smiles": "CC", "reaction": "rxn_join", "product": "product_1", "educt": "B1"},
            ],
            "B1": [
                {"index": 0, "smiles": "O", "reaction": "rxn_join", "product": "product_1", "educt": "B0"},
                {"index": 1, "smiles": "CO", "reaction": "rxn_join", "product": "product_1", "educt": "B0"},
            ],
        },
    }
    config_path = tmp_path / "config.yaml"
    write_yaml(config, config_path)
    return config, config_path


def test_get_combinations_full_mode_returns_cartesian_product(tmp_path):
    cfg, _ = make_test_config(tmp_path)

    combs = get_combinations(cfg=cfg, building_block_names=["B0", "B1"])

    assert len(combs) == 4
    assert combs[0][0]["index"] == 0
    assert combs[0][1]["index"] == 0
    assert combs[-1][0]["index"] == 1
    assert combs[-1][1]["index"] == 1


def test_get_combinations_filtered_mode_accepts_code_columns(tmp_path):
    cfg, _ = make_test_config(tmp_path)
    counts_path = tmp_path / "observed.tsv"
    counts_path.write_text("code_0\tcode_1\n1\t0\n0\t1\n")

    combs = get_combinations(
        cfg=cfg,
        building_block_names=["B0", "B1"],
        counts_path=counts_path,
    )

    assert [(comb[0]["index"], comb[1]["index"]) for comb in combs] == [(1, 0), (0, 1)]


def test_get_combinations_filtered_mode_respects_top_n(tmp_path):
    cfg, _ = make_test_config(tmp_path)
    counts_path = tmp_path / "observed.tsv"
    counts_path.write_text("code_0\tcode_1\n1\t1\n0\t0\n")

    combs = get_combinations(
        cfg=cfg,
        building_block_names=["B0", "B1"],
        counts_path=counts_path,
        top_n=1,
    )

    assert [(comb[0]["index"], comb[1]["index"]) for comb in combs] == [(1, 1)]


def test_get_combinations_filtered_mode_rejects_missing_columns(tmp_path):
    cfg, _ = make_test_config(tmp_path)
    counts_path = tmp_path / "bad.tsv"
    counts_path.write_text("B0\n0\n")

    with pytest.raises(AssertionError, match="expects the same `code_\\*` columns"):
        get_combinations(
            cfg=cfg,
            building_block_names=["B0", "B1"],
            counts_path=counts_path,
        )


def test_get_combinations_filtered_mode_rejects_out_of_range_code(tmp_path):
    cfg, _ = make_test_config(tmp_path)
    counts_path = tmp_path / "bad.tsv"
    counts_path.write_text("code_0\tcode_1\n2\t0\n")

    with pytest.raises(AssertionError, match="outside whitelist range"):
        get_combinations(
            cfg=cfg,
            building_block_names=["B0", "B1"],
            counts_path=counts_path,
        )


def test_enumerate_filtered_mode_writes_named_library_and_selected_rows(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)
    counts_path = tmp_path / "observed.tsv"
    counts_path.write_text("code_0\tcode_1\n1\t0\n0\t1\n")

    class DummyFigure:
        def savefig(self, *_args, **_kwargs):
            return None

    class DummyAxes:
        figure = DummyFigure()

    monkeypatch.setattr("delt_hit.cli.library.api.visualize_reaction_graph", lambda _g: DummyAxes())

    library = Library()
    library.enumerate(
        config_path=config_path,
        counts_path=counts_path,
        top_n=1,
        library_name="observed_hits",
    )

    library_path = tmp_path / "mini" / "library" / "observed_hits.parquet"
    assert library_path.exists()

    df = pd.read_parquet(library_path)
    assert df[["code_0", "code_1"]].to_dict("records") == [{"code_0": 1, "code_1": 0}]


def test_enumerate_filtered_mode_requires_library_name(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)
    counts_path = tmp_path / "observed.tsv"
    counts_path.write_text("code_0\tcode_1\n1\t0\n")

    monkeypatch.setattr("delt_hit.cli.library.api.visualize_reaction_graph", lambda _g: None)

    with pytest.raises(AssertionError, match="library_name"):
        Library().enumerate(config_path=config_path, counts_path=counts_path)
