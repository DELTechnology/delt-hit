import json
from pathlib import Path

import click
import numpy as np
import pandas as pd

BASES = {'A', 'T', 'C', 'G'}
rng = np.random.default_rng()


def read_struct_file(path_to_struct_file):
    path_to_struct_file = Path(path_to_struct_file)
    with open(path_to_struct_file, 'r') as f:
        lines = f.readlines()

    file_name = lines[0].strip()
    path_to_output = path_to_struct_file.parent / file_name
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
                path_to_output=path_to_output,
                elements=elements)


def generate_reads(config):
    number_of_reads = config['number_of_reads']
    elements = config['elements']
    regions = elements[['start', 'end', 'barcodes']].to_records(index=False)

    read_length = elements.end.max() + 1
    n = 0
    reads = []
    while n < number_of_reads:
        read = np.array(['X'] * read_length)
        for region in regions:
            start, end, barcodes = region
            barcode = rng.choice(barcodes)
            assert end - start + 1 == len(barcode)
            read[start:(end + 1)] = barcode
        reads.append(read)
        n += 1
    return reads


def create_all_barcode_regions_single_position_error(reads, elements, relative_position):
    # NOTE: this only works as long as all sequences have the same length and there were no insertions or deletions
    reads_stacked = np.stack(reads)
    regions = elements[elements.region_type.isin(['B'])]
    positions = regions.start + relative_position
    reads_stacked[:, positions.values] = 'X'
    return [i for i in reads_stacked]


def create_insertion_at_read_start_error(reads, elements, insertion_length, min_length=None, max_length=None, p=None):
    if insertion_length == 'random':
        inserts = [np.array(['X'] * l) for l in rng.choice(range(min_length, max_length+1), size=len(reads), p=p)]
        reads = [np.insert(read, 0, insert) for read, insert in zip(reads, inserts)]
    else:
        insert = np.array(['X'] * insertion_length)
        reads = list(map(lambda x: np.insert(x, 0, insert), reads))
    return reads


def create_errors(config, reads):
    errors = config['errors']
    elements = config['elements']
    for error in errors:
        if error['type'] == 'all_barcode_regions_single_position':
            reads = create_all_barcode_regions_single_position_error(reads, elements, **error['kwargs'])
        elif error['type'] == 'insertion_at_read_start':
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


@click.command()
@click.argument('path_to_config', type=click.Path(exists=True))
def main(path_to_config):
    path_to_config = Path(path_to_config)
    with open(path_to_config, 'r') as f:
        config = json.load(f)
    path_to_struct_file, number_of_reads = config['path_to_struct_file'], config['number_of_reads']
    struct_dict = read_struct_file(path_to_struct_file)
    config.update(struct_dict)
    reads = generate_reads(config)
    reads = create_errors(config, reads)
    generate_fastq_file(config, reads)


if __name__ == '__main__':
    main()
