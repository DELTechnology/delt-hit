from delt_core.cli.library import Library
from pathlib import Path

lib = Library()

path = Path('/Volumes/T7/experiments/experiment-GB/config.yaml')
libraries = lib.enumerate(config_path=path)
