import gzip
import json

import pandas as pd
from tqdm import tqdm


def compute_counts(
        input_file: str,
        output_dir: str,
) -> None:

    # Path to the gzip file
    input_dir = input_file.parent
    number_of_reads = json.load(open(
        sorted((input_dir).glob('*.cutadapt.json'))[-1]
    ))['read_counts']['output']

    # Open the gzip file
    from collections import defaultdict, Counter
    selections = defaultdict(lambda: defaultdict(int))
    with gzip.open(input_file, 'rt') as f:
        # Read the file line by line
        for line in tqdm(f, total=number_of_reads, ncols=100):
            # Split the line by '-'
            _, *adapters = line.strip().split('?')
            selection_id = '_'.join([i.split('.')[-1] for i in filter(lambda x: 'S' in x, adapters)])
            barcodes = tuple(i.split('.')[-1] for i in filter(lambda x: 'B' in x, adapters))

            selections[selection_id][barcodes] += 1

    for selection_id, counts in tqdm(selections.items(), ncols=100):
        # counts = [(*i, j) for i,j in zip(counts.keys(), counts.values())]
        counts = [(j, *i) for i,j in zip(counts.keys(), counts.values())]
        # df = pd.DataFrame.from_records(counts, columns=['barcode1', 'barcode2', 'count'])
        df = pd.DataFrame.from_records(counts, columns=['count', 'barcode1', 'barcode2'])
        df = df.astype(int)
        df.sort_values(['barcode1', 'barcode2'], inplace=True)

        path_counts = output_dir / f'{selection_id}.tsv'
        df.to_csv(path_counts, index=False, sep='\t')

