import gzip
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
    
    c.compute_smiles(libraries, input_path, Path(output_path))


def merge_libraries(
        input_path: tuple,
) -> None:
    libraries = []
    for library in input_path:
        libraries += [c.load_data(library)]
    c.merge_excel_files(libraries, input_path)


def compute_properties(
        input_file: Path,
) -> None:
    input_file = Path(input_file)
    output_dir = input_file.parent.parent / 'properties'
    output_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(input_file, 'rt') as file:
        header = file.readline().split('\t')
    indices = [header.index(i) for i in header if i.split('_')[0] == 'Product']
    for i, index in enumerate(indices, 1):
        output_file = output_dir / f'properties_L{i}.txt.gz'
        c.compute_properties(input_file, index, output_file)


def plot(
        input_file: Path,
) -> None:
    input_file = Path(input_file)
    output_dir = input_file.parent / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    c.plot_properties(input_file, output_dir)

