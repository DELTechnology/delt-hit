# %%
import json
from pathlib import Path
import typing as tp

import click
import numpy as np
import pandas as pd

# %%
BASES = ['A', 'T', 'C', 'G']
rng = np.random.default_rng()


def read_struct_file(path_to_struct_file):
    path_to_struct_file = Path(path_to_struct_file)
    with open(path_to_struct_file, 'r') as f:
        lines = f.readlines()

    file_name = lines[0].strip()
    elements = list(map(lambda x: x.strip().split('\t'), lines[2:]))
    elements = pd.DataFrame.from_records(elements, columns=['start', 'end', 'region_type', 'path_to_barcode_list'])
    elements['start'] = elements['start'].astype(int) - 1
    elements['end'] = elements['end'].astype(int) - 1
    elements['path_to_barcode_list'] = [path_to_struct_file.parent / i for i in elements['path_to_barcode_list']]

    barcodes = [open(i, 'r').readlines() for i in elements['path_to_barcode_list']]
    barcodes = map(lambda x: list(map(lambda y: y.strip(), x)), barcodes)
    barcodes = map(lambda x: filter(lambda y: len(y) > 0, x), barcodes)
    barcodes = map(lambda x: list(map(lambda y: list(y), x)), barcodes)
    barcodes = list(barcodes)
    elements['barcodes'] = barcodes

    elements = elements.sort_values('start')
    return dict(file_name=file_name,
                elements=elements)


def generate_reads(config):
    number_of_reads = config['number_of_reads']
    elements = config['elements']
    regions = elements[['start', 'end', 'barcodes']].to_records(index=False)

    read_length = elements.end.max() + 1
    reads = []
    indices = np.full((len(regions), number_of_reads), -1)

    for i in range(number_of_reads):
        read = np.array(['X'] * read_length)

        for j, region in enumerate(regions):
            start, end, barcodes = region
            idx = rng.choice(range(0, len(barcodes)))
            barcode = barcodes[idx]
            indices[j, i] = idx
            assert end - start + 1 == len(barcode)
            read[start:(end + 1)] = barcode
        
        reads.append(read)

    return reads, indices


def create_all_barcode_regions_single_position_error(reads, elements, relative_position):
    # NOTE: this only works as long as all sequences have the same length and there were no insertions or deletions
    reads_stacked = np.stack(reads)
    regions = elements[elements.region_type.isin(['B'])]
    positions = regions.start + relative_position
    base = rng.choice(BASES)
    print(base)
    reads_stacked[:, positions.values] = base
    return [i for i in reads_stacked]


def create_insertion_at_read_start_error(reads, elements, insertion_length, min_length=None, max_length=None, p=None):
    if insertion_length == 'random':
        inserts = [np.array([rng.choice(BASES)] * l) for l in rng.choice(range(min_length, max_length+1), size=len(reads), p=p)]
        reads = [np.insert(read, 0, insert) for read, insert in zip(reads, inserts)]
    else:
        insert = np.array([rng.choice(BASES)] * insertion_length)
        reads = list(map(lambda x: np.insert(x, 0, insert), reads))
    return reads


def create_errors(config, reads):
    errors = config['errors']
    elements = config['elements']
    for error in errors:
        if error['type'] == 'all_barcode_regions_single_position':
            reads = create_all_barcode_regions_single_position_error(reads, elements, **error['kwargs'])
        elif error['type'] == 'insertion_at_read_start':
            pass
            reads = create_insertion_at_read_start_error(reads, elements, **error['kwargs'])
        else:
            raise NotImplementedError()

    return reads


def generate_fastq_file(config, reads):
    path_to_output = config['path_to_output']
    content = []
    for i, read in enumerate(reads):
        phred_scores = '~' * len(read)
        content.append(f'@simulated-read-{i}\n')
        content.append(''.join(read) + '\n')
        content.append('+\n')
        content.append(phred_scores + '\n')

    with open(path_to_output, 'w') as f:
        f.writelines(content)


def compute_counts(
        structure: tp.Dict,
        indices: np.ndarray,
):
    code_indices = [i for i, key in enumerate(structure.keys()) if key[0] == 'B']
    counts = {}
    for index in indices.T:
        if -1 in index:
            continue
        code = tuple(index[code_indices])
        counts[code] = counts.get(code, 0) + 1
    return counts


def save_counts(
        counts: tp.Dict,
        path: str,
):
    with open(path, 'w') as file:
        num_codes = len(next(iter(counts.keys())))
        header = ['Count', *[f'Code{i}' for i in range(1, num_codes + 1)], '\n']
        file.write('\t'.join(header))
        for codes, count in counts.items():
            row = '\t'.join(str(code) for code in codes)
            file.write(f'{count}\t{row}\n')





# %%
@click.command()
@click.argument('path_to_config', type=click.Path(exists=True))
def main(path_to_config):
    path_to_config = Path(path_to_config)
    with open(path_to_config, 'r') as f:
        config = json.load(f)
    path_to_struct_file = config['path_to_struct_file']
    struct_dict = read_struct_file(path_to_struct_file)
    config.update(struct_dict)
    reads, indices = generate_reads(config)
    reads = create_errors(config, reads)
    generate_fastq_file(config, reads)

    # Compute counts.
    path_to_struct_file = Path(path_to_struct_file).parent.parent / 'structure.json'
    with open(path_to_struct_file, 'r') as file:
        structure = json.load(file)
    counts = compute_counts(structure, indices)
    path_to_counts = Path(config['path_to_output']).parent / 'counts_true.txt'
    save_counts(counts, path_to_counts)


if __name__ == '__main__':
    main()
