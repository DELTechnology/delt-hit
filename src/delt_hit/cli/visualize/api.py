from pathlib import Path

import pandas as pd

from delt_hit.cli.library.api import (
    LibraryPaths,
    close_figure,
    prepare_graph_bundle,
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
        nrow: int = 10,
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
            graph: Whether to save reaction graph visualizations.
            reactions: Whether to save reaction scheme panels from SMIRKS.
            building_blocks: Whether to save building block structure grids.
            compounds: Whether to save configured compound structure grids.
        """
        output_name: str = "visualization"

        if not any([graph, reactions, building_blocks, compounds]):
            graph = True
            reactions = True
            building_blocks = True
            compounds = True

        cfg = read_yaml(config_path)
        lib_dir = self.get_library_dir(config_path=config_path)
        lib_dir.mkdir(parents=True, exist_ok=True)

        graph_bundle = prepare_graph_bundle(cfg=cfg)

        if graph:
            save_graph_visualizations(graph_bundle=graph_bundle, save_dir=lib_dir)

        if reactions:
            reactions_dir = lib_dir / output_name / 'reactions'
            reactions_dir.mkdir(parents=True, exist_ok=True)
            visualize_reaction_schemes(cfg['catalog']['reactions'], save_dir=reactions_dir)

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
                )
                ax.figure.savefig(lib_dir / output_name / f"{bb_name}.png", dpi=300)
                close_figure(ax.figure)

        if compounds:
            compound_entries = cfg['catalog'].get('compounds', {})
            compound_smiles = []
            compound_legends = []
            for name, entry in compound_entries.items():
                smiles = entry.get('smiles')
                if pd.isna(smiles):
                    continue
                compound_smiles.append(smiles)
                compound_legends.append(name)

            if compound_smiles:
                compound_ax = visualize_smiles(
                    smiles=compound_smiles,
                    legends=compound_legends,
                    title='Compounds',
                    nrow=nrow,
                )
                compound_ax.figure.savefig(lib_dir / output_name / f"compounds.png", dpi=300)
                close_figure(compound_ax.figure)
