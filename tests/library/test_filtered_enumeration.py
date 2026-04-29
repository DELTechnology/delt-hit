from pathlib import Path

import pandas as pd
import pytest

from delt_hit.cli.library.api import Library, get_combinations
from delt_hit.cli.visualize.api import Visualize
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
            "compounds": {
                "scaffold": {"smiles": "CC(=O)O"},
            },
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

    monkeypatch.setattr(
        "delt_hit.cli.library.api.visualize_reaction_graph",
        lambda _g, include_smirks=False: DummyAxes(),
    )

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


def test_visualize_enumerate_writes_expected_input_bundle(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)

    class DummyFigure:
        def savefig(self, path, *_args, **_kwargs):
            Path(path).write_text("figure")

        def tight_layout(self):
            return None

    class DummyAxes:
        figure = DummyFigure()

    monkeypatch.setattr("delt_hit.cli.visualize.api.save_graph_visualizations", lambda **_kwargs: None)
    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_reaction_schemes", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_smiles", lambda *args, **kwargs: DummyAxes())

    Visualize().enumerate(
        config_path=config_path,
    )

    viz_dir = tmp_path / "mini" / "library" / "visualization"
    assert (viz_dir / "building_blocks_B0.png").exists()
    assert (viz_dir / "building_blocks_B1.png").exists()
    assert (viz_dir / "compounds" / "scaffold.png").exists()


def test_visualize_enumerate_flags_limit_outputs(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)

    class DummyFigure:
        def savefig(self, path, *_args, **_kwargs):
            Path(path).write_text("figure")

        def tight_layout(self):
            return None

    class DummyAxes:
        figure = DummyFigure()

    monkeypatch.setattr("delt_hit.cli.visualize.api.save_graph_visualizations", lambda **_kwargs: None)
    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_reaction_schemes", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_smiles", lambda *args, **kwargs: DummyAxes())

    Visualize().enumerate(
        config_path=config_path,
        compounds=True,
    )

    viz_dir = tmp_path / "mini" / "library" / "visualization"
    assert (viz_dir / "compounds" / "scaffold.png").exists()
    assert not (viz_dir / "building_blocks_B0.png").exists()
    assert not (viz_dir / "building_blocks_B1.png").exists()


def test_visualize_enumerate_passes_tile_size_to_structure_rendering(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)
    captured_sizes = []

    class DummyFigure:
        def savefig(self, path, *_args, **_kwargs):
            Path(path).write_text("figure")

        def tight_layout(self):
            return None

    class DummyAxes:
        figure = DummyFigure()

    def capture_visualize_smiles(*args, **kwargs):
        captured_sizes.append(kwargs.get("sub_img_size"))
        return DummyAxes()

    monkeypatch.setattr("delt_hit.cli.visualize.api.save_graph_visualizations", lambda **_kwargs: None)
    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_reaction_schemes", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_smiles", capture_visualize_smiles)

    Visualize().enumerate(
        config_path=config_path,
        building_blocks=True,
        compounds=True,
        tile_size=420,
    )

    assert captured_sizes
    assert set(captured_sizes) == {(420, 420)}


def test_properties_default_mode_writes_properties_parquet(tmp_path):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "library.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles": ["CO", "CCO"]}).to_parquet(library_path, index=False)

    Library().properties(config_path=config_path)

    properties_path = tmp_path / "mini" / "library" / "properties" / "properties.parquet"
    assert properties_path.exists()

    df = pd.read_parquet(properties_path)
    assert list(df["smiles"]) == ["CO", "CCO"]
    assert "prop_mw" in df.columns


def test_properties_named_library_writes_named_properties_parquet(tmp_path):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "observed_hits.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles": ["CO"]}).to_parquet(library_path, index=False)

    Library().properties(config_path=config_path, library_name="observed_hits")

    properties_path = tmp_path / "mini" / "library" / "properties" / "observed_hits.parquet"
    assert properties_path.exists()

    df = pd.read_parquet(properties_path)
    assert list(df["smiles"]) == ["CO"]
    assert "prop_mw" in df.columns


def test_properties_library_path_overrides_library_name(tmp_path):
    _, config_path = make_test_config(tmp_path)
    library_dir = tmp_path / "mini" / "library"
    library_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame({"smiles": ["CC"]}).to_parquet(library_dir / "observed_hits.parquet", index=False)
    explicit_path = tmp_path / "custom_library.parquet"
    pd.DataFrame({"smiles": ["CO"]}).to_parquet(explicit_path, index=False)

    Library().properties(
        config_path=config_path,
        library_name="observed_hits",
        library_path=explicit_path,
    )

    properties_path = tmp_path / "properties" / "custom_library.parquet"
    assert properties_path.exists()

    df = pd.read_parquet(properties_path)
    assert list(df["smiles"]) == ["CO"]


def test_properties_named_library_requires_existing_parquet(tmp_path):
    _, config_path = make_test_config(tmp_path)

    with pytest.raises(AssertionError, match="Library file not found"):
        Library().properties(config_path=config_path, library_name="observed_hits")
