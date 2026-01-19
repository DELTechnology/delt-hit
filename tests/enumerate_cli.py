from delt_core.cli.library import Library
from pathlib import Path

lib = Library()

path = Path('/Volumes/T7/experiments/experiment-GB/config.yaml')
# path = Path('/Volumes/T7/experiments/experiment-NF2-with-precursors/config.yaml')
i = 103564
libraries = lib.enumerate(config_path=path, overwrite=True)
