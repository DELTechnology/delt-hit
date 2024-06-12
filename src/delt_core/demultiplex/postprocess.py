from collections import defaultdict
import gzip
import json
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def perform_selection(
        reads: list,
        num_reads: int,
) -> None:
    selections = defaultdict(lambda: defaultdict(int))
    for line in tqdm(reads, total=num_reads, ncols=100):
        _, *adapters = line.strip().split('?')
        selection_id = '_'.join([i.split('.')[-1] for i in filter(lambda x: 'S' in x, adapters)])
        barcodes = tuple(i.split('.')[-1] for i in filter(lambda x: 'B' in x, adapters))
        selections[selection_id][barcodes] += 1
    return selections


def save_counts(
        selections: dict,
        output_dir: Path,
) -> None:
    for selection_id, counts in tqdm(selections.items(), ncols=100):
        counts = [(j, *i) for i, j in zip(counts.keys(), counts.values())]
        df = pd.DataFrame.from_records(counts, columns=['count', 'barcode1', 'barcode2'])
        df = df.astype(int)
        df.sort_values(['barcode1', 'barcode2'], inplace=True)
        path_counts = output_dir / f'{selection_id}.tsv'
        df.to_csv(path_counts, index=False, sep='\t')


def compute_counts(
        input_file: str,
        output_dir: str,
) -> None:
    input_dir = input_file.parent
    num_reads = json.load(open(
        sorted((input_dir).glob('*.cutadapt.json'))[-1]
    ))['read_counts']['output']
    with gzip.open(input_file, 'rt') as f:
        selections = perform_selection(f, num_reads)
    save_counts(selections, output_dir)

