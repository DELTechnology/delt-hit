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


def compute_counts_cli(
        input_file: str,
        struct_file: str,
        output_file: str = None,
):
    input_file = Path(input_file)
    structure = c.read_json(struct_file)
    if not output_file:
        output_file = input_file.parent / 'counts.json'

    counts = c.compute_counts(structure, input_file)
    c.write_json(counts, output_file)
    print(sum(counts.values()))


def evaluate_cli():
    pass

