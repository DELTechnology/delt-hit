from pathlib import Path
from importlib import resources
import subprocess

from delt_core.utils import hash_dict
from delt_core.demultiplex.validation import Config


def run(
        config_file: Path,
        target = list,
        control = list,
) -> None:
    print('Not yet implemented.')

