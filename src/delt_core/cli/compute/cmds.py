import gzip
from pathlib import Path

from ... import compute as c


def compute_smiles(
        input_files: tuple[Path, Path],
) -> None:
    if len(input_files) == 1:
        output_file = f'{Path(input_files[0]).stem}_smiles.txt.gz'
    else:
        output_file = f'{Path(input_files[0]).stem}-{Path(input_files[1]).stem}_smiles.txt.gz'

    output_dir = Path(input_files[0]).parent / 'smiles'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_file

    libraries = []
    for library in input_files:
        libraries += [c.load_data(library)]
    
    c.compute_smiles(libraries, output_path)


def merge_libraries(
        input_files: tuple[Path, Path],
) -> None:
    libraries = []
    for library in input_files:
        libraries += [c.load_data(library)]
    output_dir = Path(input_files[0]).parent
    name = f'{Path(input_files[0]).stem}-{Path(input_files[1]).stem}'
    output_file = output_dir / f'{name}.xlsx'
    c.merge_excel_files(libraries, output_file)


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
        library = input_file.stem.split('_smiles')[0]
        output_file = output_dir / f'{library}_properties_L{i}.txt.gz'
        c.compute_properties(input_file, index, output_file)


def plot(
        input_file: Path,
) -> None:
    input_file = Path(input_file)
    output_dir = input_file.parent / 'plots'
    output_dir.mkdir(parents=True, exist_ok=True)
    c.plot_properties(input_file, output_dir)

