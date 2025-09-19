from delt_core.demultiplex.preprocess import generate_input_files
from pathlib import Path

config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/templates/config.yaml')
generate_input_files(config_path=config_path)