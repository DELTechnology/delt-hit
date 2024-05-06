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
        output_file = input_file.parent / 'counts.txt'

    counts = c.compute_counts(structure, input_file)

    with open(output_file, 'w') as file:
        num_codes = len(next(iter(counts.keys())))
        header = ['Count', *[f'Code{i}' for i in range(1, num_codes + 1)], '\n']
        file.write('\t'.join(header))
        for codes, count in counts.items():
            row = '\t'.join(str(code) for code in codes)
            file.write(f'{count}\t{row}\n')
    # c.write_json(counts, output_file)
    print(sum(counts.values()))


def evaluate_cli():
    pass

