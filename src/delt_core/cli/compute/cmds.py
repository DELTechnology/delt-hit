from pathlib import Path

from ... import compute as c


def compute_smiles(
        input_path: tuple,
        output_path: str = None,
) -> None:
    if not output_path:
        output_path = Path(input_path[0]).parent / 'smiles.txt.gz'

    libraries = []
    for library in input_path:
        libraries += [c.load_data(library)]
    
    c.compute_smiles(libraries, Path(output_path))


def compute_properties(
        input_path: Path,
) -> None:
    c.compute_properties(Path(input_path))

