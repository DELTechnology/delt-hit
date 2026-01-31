from pathlib import Path

from loguru import logger

from delt_hit.demultiplex.parser import config_from_excel
from delt_hit.utils import write_yaml

def init(*, excel_path: Path):
    config = config_from_excel(excel_path.expanduser().resolve())
    save_dir = Path(config['experiment']['save_dir']).expanduser().resolve()
    name = config['experiment']['name']

    save_path = save_dir / name / 'config.yaml'
    save_path.parent.mkdir(parents=True, exist_ok=True)
    write_yaml(config, save_path)
    logger.info(f'Configuration created at {save_path}')
