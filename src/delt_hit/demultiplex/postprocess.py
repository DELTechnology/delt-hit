from collections import defaultdict
import gzip
from pathlib import Path

import pandas as pd
from tqdm import tqdm

def extract_ids(line: str):
    _, *adapters = line.strip().split('?')
    selection_ids = [i.split('.')[-1] for i in filter(lambda x: 'S' in x, adapters)]
    selection_ids = tuple(map(int, selection_ids))
    barcodes = tuple(int(i.split('.')[-1]) + 1 for i in filter(lambda x: 'B' in x, adapters))
    return {'selection_ids': selection_ids, 'barcodes': barcodes}


def save_counts(counts: dict, output_dir: Path, ids_to_name: dict = None,
                as_files: bool = True, sort_by_counts: bool = True) -> None:

    num_codes = len(list(list(counts.values())[0].keys())[0])
    codon_cols = [f'code_{i}' for i in range(1, num_codes + 1)]
    columns = codon_cols + ['count']

    sort_by_cols = 'count' if sort_by_counts else codon_cols

    for selection_ids, count in tqdm(counts.items(), ncols=100):
        rows = [(*k, v, "_".join(map(str, k))) for k, v in count.items()]
        df = pd.DataFrame.from_records(rows, columns=[*columns, 'id'])
        df = df.astype({k: int for k in columns})
        df.sort_values(sort_by_cols, ascending=False, inplace=True)

        if ids_to_name is None:
            name = '-'.join(map(str, selection_ids))
        else:
            name = ids_to_name[selection_ids]

        if as_files:
            output_file = output_dir / f'{name}_counts.txt'
            output_file.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_file, index=False, sep='\t')
        else:
            selection_dir = output_dir / name
            selection_dir.mkdir(parents=True, exist_ok=True)
            output_file = selection_dir / f'counts.txt'
            df.to_csv(output_file, index=False, sep='\t')


def get_counts(*, input_path: Path, num_reads: int) -> dict:
    with gzip.open(input_path, 'rt') as f:
        counts = defaultdict(lambda: defaultdict(int))
        for line in tqdm(f, total=num_reads, ncols=100):
            ids = extract_ids(line)
            counts[ids['selection_ids']][ids['barcodes']] += 1
    return counts


