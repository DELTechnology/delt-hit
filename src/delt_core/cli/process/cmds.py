from pathlib import Path
import typing as tp

from ... import compute as c


def compute_smiles_cli(
        input_path: tp.Tuple,
        output_path: str = None,
):
    if not output_path:
        output_path = Path(input_path[0]).parent / 'smiles.txt'

    libraries = []
    for library in input_path:
        libraries += [c.load_data(library)]
    
    c.compute_smiles(libraries, Path(output_path))


def compute_counts_cli():
    c.compute_counts()


def evaluate_cli():
    pass

