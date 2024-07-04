from delt_core.quality_control.compare_output import compare_counts_with_legacy
from delt_core.demultiplex.preprocess import read_yaml
from pathlib import Path

# from delt_core.cli.demultiplex.cmds import init
#
# from delt_core.demultiplex.preprocess import convert_struct_file
# convert_struct_file('/Users/adrianomartinelli/polybox/decl-data/raw-files/structureNF2.txt')
#
# init('', '', '', '', '', '')

config = read_yaml(Path('/Users/adrianomartinelli/polybox/decl-data/data_andy/experiments/no-errors/config.yml'))
config = read_yaml(Path('/Users/adrianomartinelli/polybox/decl-data/data_andy/experiments/ignore-C1/config.yml'))
legacy_results_dir = Path('/Users/adrianomartinelli/polybox/decl-data/data_andy/evaluation_cpp')
compare_counts_with_legacy(config, legacy_results_dir)
