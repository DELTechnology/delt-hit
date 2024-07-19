from delt_core.quality_control.compare_output import compare_counts_with_legacy
from delt_core.cli.demultiplex.cmds import create_lists
from pathlib import Path
from delt_core.demultiplex.utils import Config

config_path = Path('/Users/adrianomartinelli/polybox/decl-data/db/experiments/OST1-e0-2024-07-12-20-29-10/config.yml')
config_path = Path('/Users/adrianomartinelli/polybox/decl-data/db/experiments/BAF1-2024-07-12-21-15-02/config.yml')
config = Config.from_yaml(config_path).model_dump()
legacy_results_dir = Path('/Users/adrianomartinelli/polybox/decl-data/raw-files/Evaluation_2301_704504_NF2-AG_yOST')
legacy_results_dir = Path('/Users/adrianomartinelli/polybox/decl-data/raw-files/Evaluation_2301_701501_NF2-AG_TRIM28_dBAF')

# create_lists(config_path, output_dir=None)

compare_counts_with_legacy(config, legacy_results_dir)