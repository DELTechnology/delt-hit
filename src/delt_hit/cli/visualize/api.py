from pathlib import Path

import pandas as pd
from tqdm import tqdm

from delt_hit.cli.helper import get_library_dir, get_named_library_path
from delt_hit.cli.library.api import (
    close_figure,
    prepare_graph_bundle,
    save_figure_outputs,
    visualize_reaction_graph,
    visualize_reaction_schemes,
    visualize_smiles,
)
from delt_hit.utils import read_yaml


def save_graph_visualizations(*, graph_bundle: dict, save_dir: Path, dpi: int = 300) -> None:
    """Write the standard reaction graph visualizations to disk."""
    for graph, filename in [
        (graph_bundle['bb_G'], 'reaction_graph_building_blocks'),
        (graph_bundle['add_G'], 'reaction_graph_additional'),
        (graph_bundle['G'], 'reaction_graph'),
    ]:
        ax = visualize_reaction_graph(graph)
        save_figure_outputs(ax.figure, save_dir / filename, dpi=dpi)
        close_figure(ax.figure)


class Visualize:
    def enumerate(
        self,
        *,
        config_path: Path,
        building_block_ids: list[str] | None = None,
        dpi: int = 300,
        tile_size: int = 300,
        graph: bool = False,
        reactions: bool = False,
        building_blocks: bool = False,
        compounds: bool = False,
    ):
        """Generate visualization panels for enumeration inputs.

        Args:
            config_path: Path to the YAML config file.
            building_block_ids: Optional subset of building block IDs to consider.
            dpi: Raster DPI used when exporting PNGs.
            tile_size: Pixel width and height used by RDKit for each molecule tile.
            graph: Whether to save reaction graph visualizations.
            reactions: Whether to save reaction scheme panels from SMIRKS.
            building_blocks: Whether to save building block structure panels.
            compounds: Whether to save configured compound structure grids.
        """
        if not any([graph, reactions, building_blocks, compounds]):
            graph = True
            reactions = True
            building_blocks = True
            compounds = True

        cfg = read_yaml(config_path)
        lib_dir = get_library_dir(config_path)
        lib_dir.mkdir(parents=True, exist_ok=True)
        visualization_dir = lib_dir / "visualization"
        visualization_dir.mkdir(parents=True, exist_ok=True)

        graph_bundle = prepare_graph_bundle(cfg=cfg)

        if graph:
            save_graph_visualizations(graph_bundle=graph_bundle, save_dir=visualization_dir, dpi=dpi)

        if reactions:
            reactions_dir = visualization_dir / 'reactions'
            reactions_dir.mkdir(parents=True, exist_ok=True)
            visualize_reaction_schemes(cfg['catalog']['reactions'], save_dir=reactions_dir, dpi=dpi)

        if building_blocks:
            building_block_names = sorted(graph_bundle['building_blocks'])
            if building_block_ids:
                building_block_names = [name for name in building_block_names if name in building_block_ids]

            building_blocks_dir = visualization_dir / "building_blocks"
            building_blocks_dir.mkdir(parents=True, exist_ok=True)

            for bb_name in tqdm(building_block_names, desc="Building block families"):
                whitelist = cfg['whitelists'][bb_name]
                bb_dir = building_blocks_dir / bb_name
                bb_dir.mkdir(parents=True, exist_ok=True)

                for entry in tqdm(whitelist, desc=f"{bb_name} panels", leave=False):
                    smiles = entry['smiles']
                    if pd.isna(smiles):
                        continue

                    ax = visualize_smiles(
                        smiles=[smiles],
                        legends=[f"{bb_name}:{entry['index']}"],
                        title=f'{bb_name} Building Blocks',
                        nrow=1,
                        sub_img_size=(tile_size, tile_size),
                    )
                    save_figure_outputs(ax.figure, bb_dir / str(entry['index']), dpi=dpi)
                    close_figure(ax.figure)

        if compounds:
            compound_entries = cfg['catalog'].get('compounds', {})
            compounds_dir = visualization_dir / "compounds"
            compounds_dir.mkdir(parents=True, exist_ok=True)

            for name, entry in tqdm(compound_entries.items(), desc="Compound panels"):
                smiles = entry.get('smiles')
                if pd.isna(smiles):
                    continue
                compound_ax = visualize_smiles(
                    smiles=[smiles],
                    legends=[name],
                    title=name,
                    nrow=1,
                    sub_img_size=(tile_size, tile_size),
                )
                save_figure_outputs(compound_ax.figure, compounds_dir / name, dpi=dpi)
                close_figure(compound_ax.figure)

    def library(
        self,
        *,
        config_path: Path,
        library_name: str,
        dpi: int = 300,
        tile_size: int = 300,
    ):
        """Generate one structure panel per entry for a named library parquet.

        Args:
            config_path: Path to the YAML config file.
            library_name: Named library parquet stem inside the experiment library dir.
            dpi: Raster DPI used when exporting PNGs.
            tile_size: Pixel width and height used by RDKit for each molecule tile.
        """
        lib_path = get_named_library_path(config_path, library_name)
        assert lib_path.exists(), f"Library file not found at {lib_path}"

        df = pd.read_parquet(lib_path)
        assert "smiles" in df.columns, "Dual-display libraries with `smiles_a`/`smiles_b` are not supported by `visualize library`"
        legend_df = df.loc[df["smiles"].notna()]
        legends = build_code_legends(legend_df)
        filenames = build_code_filenames(legend_df)

        visualization_dir = get_library_dir(config_path) / "visualization"
        visualization_dir.mkdir(parents=True, exist_ok=True)
        library_dir = visualization_dir / "libraries" / library_name
        library_dir.mkdir(parents=True, exist_ok=True)

        smiles = legend_df["smiles"].tolist()
        for i, smiles_value in tqdm(enumerate(smiles), total=len(smiles), desc=f"{library_name} panels"):
            legend = legends[i] if legends is not None else str(i)
            ax = visualize_smiles(
                smiles=[smiles_value],
                legends=[legend],
                title=library_name,
                nrow=1,
                sub_img_size=(tile_size, tile_size),
            )
            filename = filenames[i] if filenames is not None else str(i)
            save_figure_outputs(ax.figure, library_dir / filename, dpi=dpi)
            close_figure(ax.figure)


def build_code_legends(data: pd.DataFrame) -> list[str] | None:
    """Build ``code_*`` legends for a library dataframe when available."""
    code_columns = [col for col in data.columns if col.startswith("code_")]
    code_columns = sorted(code_columns, key=lambda name: int(name.split("_", maxsplit=1)[1]))
    if not code_columns:
        return None
    return data[code_columns].astype(str).agg(":".join, axis=1).tolist()


def build_code_filenames(data: pd.DataFrame) -> list[str] | None:
    """Build per-entry filenames from ``code_*`` columns when available."""
    code_columns = [col for col in data.columns if col.startswith("code_")]
    code_columns = sorted(code_columns, key=lambda name: int(name.split("_", maxsplit=1)[1]))
    if not code_columns:
        return None

    filenames = []
    for _, row in data[code_columns].iterrows():
        parts = [f"B{col.split('_', maxsplit=1)[1]}={row[col]}" for col in code_columns]
        filenames.append("-".join(parts))
    return filenames
