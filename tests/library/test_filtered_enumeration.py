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
        "structure": [
            {"name": "B0", "type": "building_block", "strand": None},
            {"name": "B1", "type": "building_block", "strand": None},
        ],
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


def make_dual_display_test_config(tmp_path: Path) -> tuple[dict, Path]:
    config = {
        "experiment": {
            "name": "dual",
            "save_dir": str(tmp_path),
        },
        "library": {
            "bb_edges": [
                ["B0", "rxn_a"], ["B1", "rxn_a"], ["rxn_a", "product_a"],
                ["B2", "rxn_b"], ["B3", "rxn_b"], ["rxn_b", "product_b"],
            ],
            "other_edges": [],
            "building_blocks": ["B0", "B1", "B2", "B3"],
            "products": ["product_a", "product_b"],
            "reactions": ["rxn_a", "rxn_b"],
        },
        "catalog": {
            "reactions": {
                "rxn_a": {"smirks": "[C:1].[O:2]>>[C:1][O:2]"},
                "rxn_b": {"smirks": "[C:1].[O:2]>>[C:1][O:2]"},
            },
            "compounds": {},
        },
        "structure": [
            {"name": "B0", "type": "building_block", "strand": "A"},
            {"name": "B1", "type": "building_block", "strand": "A"},
            {"name": "B2", "type": "building_block", "strand": "B"},
            {"name": "B3", "type": "building_block", "strand": "B"},
        ],
        "whitelists": {
            "B0": [
                {"index": 0, "smiles": "C", "reaction": "rxn_a", "product": "product_a", "educt": "B1"},
                {"index": 1, "smiles": "CC", "reaction": "rxn_a", "product": "product_a", "educt": "B1"},
            ],
            "B1": [
                {"index": 0, "smiles": "O", "reaction": "rxn_a", "product": "product_a", "educt": "B0"},
            ],
            "B2": [
                {"index": 0, "smiles": "C", "reaction": "rxn_b", "product": "product_b", "educt": "B3"},
            ],
            "B3": [
                {"index": 0, "smiles": "O", "reaction": "rxn_b", "product": "product_b", "educt": "B2"},
                {"index": 1, "smiles": "CO", "reaction": "rxn_b", "product": "product_b", "educt": "B2"},
            ],
        },
    }
    config_path = tmp_path / "dual_config.yaml"
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


def test_get_combinations_dual_display_config_keeps_existing_product_order(tmp_path):
    cfg, _ = make_dual_display_test_config(tmp_path)

    combs = get_combinations(cfg=cfg, building_block_names=["B0", "B1", "B2", "B3"])

    assert len(combs) == 4
    assert [(comb[0]["index"], comb[1]["index"], comb[2]["index"], comb[3]["index"]) for comb in combs] == [
        (0, 0, 0, 0),
        (0, 0, 0, 1),
        (1, 0, 0, 0),
        (1, 0, 0, 1),
    ]


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


def test_enumerate_parallel_workers_match_serial_output(tmp_path):
    _, config_path = make_test_config(tmp_path)
    counts_path = tmp_path / "observed.tsv"
    counts_path.write_text("code_0\tcode_1\n1\t0\n0\t1\n1\t1\n0\t0\n")

    Library().enumerate(
        config_path=config_path,
        counts_path=counts_path,
        library_name="serial",
        num_workers=1,
    )
    Library().enumerate(
        config_path=config_path,
        counts_path=counts_path,
        library_name="parallel",
        num_workers=2,
    )

    serial = pd.read_parquet(tmp_path / "mini" / "library" / "serial.parquet")
    parallel = pd.read_parquet(tmp_path / "mini" / "library" / "parallel.parquet")
    pd.testing.assert_frame_equal(parallel, serial)


def test_enumerate_rejects_nonpositive_worker_count(tmp_path):
    _, config_path = make_test_config(tmp_path)

    with pytest.raises(AssertionError, match="positive integer"):
        Library().enumerate(config_path=config_path, num_workers=0)


def test_enumerate_writes_reaction_graphs_under_visualization_dir(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)
    captured_save_dir = None

    def capture_save_graph_visualizations(**kwargs):
        nonlocal captured_save_dir
        captured_save_dir = kwargs["save_dir"]

    monkeypatch.setattr(
        "delt_hit.cli.library.api.save_graph_visualizations",
        capture_save_graph_visualizations,
    )

    Library().enumerate(config_path=config_path)

    assert captured_save_dir == tmp_path / "mini" / "library" / "visualization"


def test_enumerate_dual_display_writes_paired_smiles_columns(tmp_path, monkeypatch):
    _, config_path = make_dual_display_test_config(tmp_path)

    class DummyFigure:
        def savefig(self, *_args, **_kwargs):
            return None

    class DummyAxes:
        figure = DummyFigure()

    monkeypatch.setattr("delt_hit.cli.library.api.visualize_reaction_graph", lambda _g: DummyAxes())

    Library().enumerate(config_path=config_path)

    library_path = tmp_path / "dual" / "library" / "library.parquet"
    df = pd.read_parquet(library_path)

    assert list(df.columns) == ["code_0", "code_1", "code_2", "code_3", "smiles_a", "smiles_b"]
    assert df[["code_0", "code_1", "code_2", "code_3"]].to_dict("records") == [
        {"code_0": 0, "code_1": 0, "code_2": 0, "code_3": 0},
        {"code_0": 0, "code_1": 0, "code_2": 0, "code_3": 1},
        {"code_0": 1, "code_1": 0, "code_2": 0, "code_3": 0},
        {"code_0": 1, "code_1": 0, "code_2": 0, "code_3": 1},
    ]
    assert list(df["smiles_a"]) == ["CO", "CO", "CCO", "CCO"]
    assert list(df["smiles_b"]) == ["CO", "COC", "CO", "COC"]


def test_enumerate_single_strand_subset_keeps_global_code_names(tmp_path, monkeypatch):
    cfg, _ = make_dual_display_test_config(tmp_path)
    strand_b_cfg = {
        **cfg,
        "library": {
            **cfg["library"],
            "bb_edges": [["B2", "rxn_b"], ["B3", "rxn_b"], ["rxn_b", "product_b"]],
            "other_edges": [],
            "building_blocks": ["B2", "B3"],
            "products": ["product_b"],
            "reactions": ["rxn_b"],
        },
        "catalog": {
            **cfg["catalog"],
            "reactions": {"rxn_b": cfg["catalog"]["reactions"]["rxn_b"]},
        },
        "structure": [entry for entry in cfg["structure"] if entry["name"] in {"B2", "B3"}],
        "whitelists": {name: cfg["whitelists"][name] for name in ["B2", "B3"]},
    }

    class DummyFigure:
        def savefig(self, *_args, **_kwargs):
            return None

    class DummyAxes:
        figure = DummyFigure()

    monkeypatch.setattr("delt_hit.cli.library.api.visualize_reaction_graph", lambda _g: DummyAxes())

    from delt_hit.cli.library.api import enumerate_single_strand, prepare_graph_bundle

    df = enumerate_single_strand(
        cfg=strand_b_cfg,
        graph_bundle=prepare_graph_bundle(cfg=strand_b_cfg),
        building_block_names=["B2", "B3"],
    )

    assert list(df.columns) == ["code_2", "code_3", "smiles"]


def test_enumerate_dual_display_filtered_mode_uses_full_code_rows(tmp_path, monkeypatch):
    _, config_path = make_dual_display_test_config(tmp_path)
    counts_path = tmp_path / "observed.tsv"
    counts_path.write_text("code_0\tcode_1\tcode_2\tcode_3\n1\t0\t0\t1\n0\t0\t0\t0\n")

    class DummyFigure:
        def savefig(self, *_args, **_kwargs):
            return None

    class DummyAxes:
        figure = DummyFigure()

    monkeypatch.setattr("delt_hit.cli.library.api.visualize_reaction_graph", lambda _g: DummyAxes())

    Library().enumerate(
        config_path=config_path,
        counts_path=counts_path,
        library_name="observed_hits",
    )

    library_path = tmp_path / "dual" / "library" / "observed_hits.parquet"
    df = pd.read_parquet(library_path)

    assert df[["code_0", "code_1", "code_2", "code_3"]].to_dict("records") == [
        {"code_0": 0, "code_1": 0, "code_2": 0, "code_3": 0},
        {"code_0": 0, "code_1": 0, "code_2": 0, "code_3": 1},
        {"code_0": 1, "code_1": 0, "code_2": 0, "code_3": 0},
        {"code_0": 1, "code_1": 0, "code_2": 0, "code_3": 1},
    ]


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
    assert (viz_dir / "building_blocks" / "B0" / "0.png").exists()
    assert (viz_dir / "building_blocks" / "B0" / "1.png").exists()
    assert (viz_dir / "building_blocks" / "B1" / "0.png").exists()
    assert (viz_dir / "building_blocks" / "B1" / "1.png").exists()
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
    assert not (viz_dir / "building_blocks" / "B0").exists()
    assert not (viz_dir / "building_blocks" / "B1").exists()


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
    assert len(captured_sizes) == 5


def test_visualize_enumerate_omits_structure_titles_and_legends(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)
    captured_titles = []
    captured_legends = []

    class DummyFigure:
        def savefig(self, path, *_args, **_kwargs):
            Path(path).write_text("figure")

        def tight_layout(self):
            return None

    class DummyAxes:
        figure = DummyFigure()

    def capture_visualize_smiles(*args, **kwargs):
        captured_titles.append(kwargs.get("title"))
        captured_legends.append(kwargs.get("legends"))
        return DummyAxes()

    monkeypatch.setattr("delt_hit.cli.visualize.api.save_graph_visualizations", lambda **_kwargs: None)
    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_reaction_schemes", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_smiles", capture_visualize_smiles)

    Visualize().enumerate(
        config_path=config_path,
        building_blocks=True,
        compounds=True,
    )

    assert captured_titles
    assert set(captured_titles) == {""}
    assert set(captured_legends) == {None}


def test_visualize_library_writes_named_library_panels(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "observed_hits.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "smiles": ["CO", "CCO"],
            "code_0": [1, 0],
            "code_1": [0, 1],
        }
    ).to_parquet(library_path, index=False)

    class DummyFigure:
        def savefig(self, path, *_args, **_kwargs):
            Path(path).write_text("figure")

        def tight_layout(self):
            return None

    class DummyAxes:
        figure = DummyFigure()

    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_smiles", lambda *args, **kwargs: DummyAxes())

    Visualize().library(
        config_path=config_path,
        library_name="observed_hits",
    )

    viz_dir = tmp_path / "mini" / "library" / "visualization"
    assert (viz_dir / "libraries" / "observed_hits" / "B0=1-B1=0.png").exists()
    assert (viz_dir / "libraries" / "observed_hits" / "B0=0-B1=1.png").exists()


def test_visualize_library_passes_clean_rendering_options(tmp_path, monkeypatch):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "observed_hits.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "smiles": ["CO", "CCO"],
            "code_0": [1, 0],
            "code_1": [0, 1],
        }
    ).to_parquet(library_path, index=False)
    captured_kwargs = []

    class DummyFigure:
        def savefig(self, path, *_args, **_kwargs):
            Path(path).write_text("figure")

        def tight_layout(self):
            return None

    class DummyAxes:
        figure = DummyFigure()

    def capture_visualize_smiles(*args, **kwargs):
        captured_kwargs.append(kwargs)
        return DummyAxes()

    monkeypatch.setattr("delt_hit.cli.visualize.api.visualize_smiles", capture_visualize_smiles)

    Visualize().library(
        config_path=config_path,
        library_name="observed_hits",
        tile_size=420,
    )

    assert {entry["legends"] for entry in captured_kwargs} == {None}
    assert {entry["title"] for entry in captured_kwargs} == {""}
    assert {entry["sub_img_size"] for entry in captured_kwargs} == {(420, 420)}
    assert {entry["nrow"] for entry in captured_kwargs} == {1}


def test_visualize_library_requires_existing_named_library(tmp_path):
    _, config_path = make_test_config(tmp_path)

    with pytest.raises(AssertionError, match="Library file not found"):
        Visualize().library(config_path=config_path, library_name="observed_hits")


def test_visualize_library_rejects_dual_display_outputs(tmp_path):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "observed_hits.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles_a": ["CO"], "smiles_b": ["NCl"]}).to_parquet(library_path, index=False)

    with pytest.raises(AssertionError, match="Dual-display libraries"):
        Visualize().library(config_path=config_path, library_name="observed_hits")


def test_properties_default_mode_writes_properties_parquet(tmp_path):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "library.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles": ["CO", "CCO"]}).to_parquet(library_path, index=False)

    Library().properties(config_path=config_path)

    properties_dir = tmp_path / "mini" / "library" / "properties" / "library"
    properties_path = properties_dir / "properties.parquet"
    assert properties_path.exists()
    assert (properties_dir / "prop_mw.png").exists()

    df = pd.read_parquet(properties_path)
    assert list(df["smiles"]) == ["CO", "CCO"]
    assert "prop_mw" in df.columns


def test_properties_named_library_writes_named_properties_parquet(tmp_path):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "observed_hits.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles": ["CO"]}).to_parquet(library_path, index=False)

    Library().properties(config_path=config_path, library_name="observed_hits")

    properties_dir = tmp_path / "mini" / "library" / "properties" / "observed_hits"
    properties_path = properties_dir / "properties.parquet"
    assert properties_path.exists()
    assert (properties_dir / "prop_mw.png").exists()

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

    properties_dir = tmp_path / "properties" / "custom_library"
    properties_path = properties_dir / "properties.parquet"
    assert properties_path.exists()
    assert (properties_dir / "prop_mw.png").exists()

    df = pd.read_parquet(properties_path)
    assert list(df["smiles"]) == ["CO"]


def test_properties_named_library_requires_existing_parquet(tmp_path):
    _, config_path = make_test_config(tmp_path)

    with pytest.raises(AssertionError, match="Library file not found"):
        Library().properties(config_path=config_path, library_name="observed_hits")


def test_properties_reject_dual_display_outputs(tmp_path):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "library.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles_a": ["CO"], "smiles_b": ["NCl"]}).to_parquet(library_path, index=False)

    with pytest.raises(AssertionError, match="Dual-display libraries"):
        Library().properties(config_path=config_path)


def test_represent_rejects_dual_display_outputs(tmp_path):
    _, config_path = make_test_config(tmp_path)
    library_path = tmp_path / "mini" / "library" / "library.parquet"
    library_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles_a": ["CO"], "smiles_b": ["NCl"]}).to_parquet(library_path, index=False)

    with pytest.raises(AssertionError, match="Dual-display libraries"):
        Library().represent(config_path=config_path)


def test_complete_reaction_graph_logging_does_not_raise_unboundlocalerror(tmp_path):
    from delt_hit.cli.library.api import complete_reaction_graph, get_reaction_graph

    graph = get_reaction_graph(
        steps=[("broken_rxn", "product_1")],
        reactions={"broken_rxn": {"smirks": "[C:1].[O:2]>>[C:1][O:2]"}},
        compounds={},
        products={"product_1": {"smiles": None}},
    )

    with pytest.raises(SystemExit):
        complete_reaction_graph(graph, errors="raise")
