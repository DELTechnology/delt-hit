from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from .utils import read_txt, write_txt, read_yaml
from delt_core.cli.demultiplex.cmds import create_lists
from delt_core.demultiplex.postprocess import extract_ids, get_selection_ids, save_counts
from delt_core.demultiplex.preprocess import read_yaml


BASES = ['A', 'T', 'C', 'G']
rng = np.random.default_rng()


def read_structure(
        config_file: Path,
) -> dict:
    config = read_yaml(config_file)
    root = Path(config['Root'])
    output_dir = root / 'codon_lists'
    structure = create_lists(config_file, output_dir=output_dir)
    elements = []
    start = 0
    for element, values in structure.items():
        sequences = read_txt(values['Path'])
        sequences = list(map(lambda x: [*x.strip()], sequences))
        length = len(sequences[0])
        end = start + length - 1
        elements += [[element, element[0], start, end, sequences]]
        start += length
    columns = ['region', 'region_type', 'start', 'end', 'barcodes']
    elements = pd.DataFrame.from_records(elements, columns=columns)
    return (dict(elements=elements))


def generate_reads(config):
    num_reads = config['Simulation']['NumReads']
    elements = config['elements']
    regions = elements[['region', 'start', 'end', 'barcodes']].to_records(index=False)
    read_length = elements.end.max() + 1
    reads, reads_info = [], []
    for i in range(num_reads):
        read = np.array(['X'] * read_length)
        read_info = f'@simulated-read-{i} '
        for j, region in enumerate(regions):
            region, start, end, barcodes = region
            idx = rng.choice(range(0, len(barcodes)))
            barcode = barcodes[idx]
            assert end - start + 1 == len(barcode)
            read[start:(end + 1)] = barcode
            read_info += f'?{j}-{region}.{idx}'
        reads.append(read)
        reads_info.append(read_info)
    return reads, reads_info


def create_all_barcode_regions_single_position_error(reads, elements, relative_position):
    # NOTE: this only works as long as all sequences have the same length and there were no insertions or deletions
    reads_stacked = np.stack(reads)
    regions = elements[elements.region_type.isin(['B'])]
    positions = regions.start + relative_position
    base = rng.choice(BASES)
    reads_stacked[:, positions.values] = base
    return [i for i in reads_stacked]


def create_insertion_at_read_start_error(reads, insertion_length, min_length=None, max_length=None, p=None):
    if insertion_length == 'random':
        inserts = [np.array([rng.choice(BASES)] * l) for l in rng.choice(range(min_length, max_length+1), size=len(reads), p=p)]
        reads = [np.insert(read, 0, insert) for read, insert in zip(reads, inserts)]
    else:
        insert = np.array([rng.choice(BASES)] * insertion_length)
        reads = list(map(lambda x: np.insert(x, 0, insert), reads))
    return reads


def create_errors(config, reads):
    errors = config['Simulation']['Errors']
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
    path_to_output = Path(config['Simulation']['OutputFile'])
    content = []
    for i, read in enumerate(reads):
        phred_scores = '~' * len(read)
        content.append(f'@simulated-read-{i}\n')
        content.append(''.join(read) + '\n')
        content.append('+\n')
        content.append(phred_scores + '\n')
    write_txt(content, path_to_output)


def run_simulation(
        config_file: Path,
) -> None:
    
    # Generate FASTQ file.
    struct_dict = read_structure(config_file)
    config = read_yaml(config_file)
    config.update(struct_dict)
    reads, reads_info = generate_reads(config)
    reads = create_errors(config, reads)
    generate_fastq_file(config, reads)

    # Generate table for counts.
    root = Path(config['Root'])
    output_dir = root / 'evaluations_true'
    output_dir.mkdir(parents=True, exist_ok=True)
    counts = defaultdict(lambda: defaultdict(int))
    for line in reads_info:
        ids = extract_ids(line)
        counts[ids['selection_ids']][ids['barcodes']] += 1
    list_of_selection_primer_ids = list(counts.keys())
    selection_ids = get_selection_ids(list_of_selection_primer_ids, config)
    counts = {selection_id: val for selection_id, val in zip(selection_ids, counts.values())}
    save_counts(counts, output_dir, config)


if __name__ == '__main__':
    
    config_file = 'config/config.yml'
    run_simulation(config_file)

