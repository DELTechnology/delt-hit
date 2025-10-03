from pathlib import Path
import yaml

from delt_core.demultiplex import experiment_from_excel, whitelists_from_excel, catalog_from_excel, config_from_excel, structure_from_excel

path = Path('/Users/adrianomartinelli/projects/delt/delt-core/templates/library.xlsx')
save_config = config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/templates/config.yaml')

experiment_from_excel(path)
catalog_from_excel(path)
whitelists_from_excel(path)
structure_from_excel(path)
config = config_from_excel(path)

with open(save_config, 'w') as f:
    yaml.dump(config, f, sort_keys=False, indent=2)