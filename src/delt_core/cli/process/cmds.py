from pathlib import Path
import typing as tp

from ... import process as p


def compute_smiles_cli(
        input_path: tp.Tuple,
        output_path: str = None,
):
    if not output_path:
        output_path = Path(input_path[0]).parent / 'smiles.txt'

    libraries = []
    for library in input_path:
        libraries += [p.load_data(library)]
    
    p.compute_smiles(libraries, Path(output_path))


def compute_counts_cli():
    p.compute_counts()


def evaluate_cli():
    pass

