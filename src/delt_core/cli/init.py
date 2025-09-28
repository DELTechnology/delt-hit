from pathlib import Path

from loguru import logger

from delt_core.demultiplex.parser import config_from_excel
from delt_core.utils import write_yaml

excel_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/NF.xlsx')

def init(*, excel_path: Path):
    config = config_from_excel(excel_path)
    save_dir = Path(config['experiment']['save_dir']).expanduser().resolve()
    name = config['experiment']['name']

    save_path = save_dir / name / 'config.yaml'
    save_path.parent.mkdir(parents=True, exist_ok=True)
    write_yaml(config, save_path)
    logger.info(f'Configuration created at {save_path}')
