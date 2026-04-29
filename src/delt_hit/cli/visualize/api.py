from pathlib import Path

import pandas as pd

from delt_hit.cli.library.api import (
    LibraryPaths,
    close_figure,
    prepare_graph_bundle,
    save_figure_outputs,
    save_graph_visualizations,
    visualize_reaction_schemes,
    visualize_smiles,
)
from delt_hit.utils import read_yaml


class Visualize(LibraryPaths):
    def enumerate(
        self,
        *,
        config_path: Path,
        building_block_ids: list[str] | None = None,
        nrow: int = 25,
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
            nrow: Number of molecules per row in structure grids.
            dpi: Raster DPI used when exporting PNGs.
            tile_size: Pixel width and height used by RDKit for each molecule tile.
            graph: Whether to save reaction graph visualizations.
            reactions: Whether to save reaction scheme panels from SMIRKS.
            building_blocks: Whether to save building block structure grids.
            compounds: Whether to save configured compound structure grids.
        """
        if not any([graph, reactions, building_blocks, compounds]):
            graph = True
            reactions = True
            building_blocks = True
            compounds = True

        cfg = read_yaml(config_path)
        lib_dir = self.get_library_dir(config_path=config_path)
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

            for bb_name in building_block_names:
                whitelist = cfg['whitelists'][bb_name]
                smiles = [entry['smiles'] for entry in whitelist if not pd.isna(entry['smiles'])]
                legends = [f"{bb_name}:{entry['index']}" for entry in whitelist if not pd.isna(entry['smiles'])]

                if not smiles:
                    continue

                ax = visualize_smiles(
                    smiles=smiles,
                    legends=legends,
                    title=f'{bb_name} Building Blocks',
                    nrow=nrow,
                    sub_img_size=(tile_size, tile_size),
                )
                save_figure_outputs(ax.figure, visualization_dir / f"building_blocks_{bb_name}", dpi=dpi)
                close_figure(ax.figure)

        if compounds:
            compound_entries = cfg['catalog'].get('compounds', {})
            compounds_dir = visualization_dir / "compounds"
            compounds_dir.mkdir(parents=True, exist_ok=True)

            for name, entry in compound_entries.items():
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
