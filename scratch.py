from delt_core.quality_control.compare_output import compare_counts_with_legacy
from pathlib import Path
from delt_core.demultiplex.utils import Config

config = Config.from_yaml(Path('/Users/adrianomartinelli/polybox/decl-data/db/experiments/OST1-e0-2024-07-12-19-15-53/config.yml')).model_dump()
legacy_results_dir = Path('/Users/adrianomartinelli/polybox/decl-data/raw-files/Evaluation_2301_704504_NF2-AG_yOST')
compare_counts_with_legacy(config, legacy_results_dir)
