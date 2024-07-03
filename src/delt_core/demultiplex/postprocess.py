from collections import defaultdict
import gzip
import json
from pathlib import Path
import shutil

import pandas as pd
from tqdm import tqdm

from .preprocess import read_yaml, get_selections, hash_dict


def get_selection_id(
        ids: list,
        config_file: Path,
) -> int:
    config = read_yaml(config_file)
    root = Path(config['Root'])
    keys = config['Structure'].keys()
    keys_s = [key for key in keys if key.startswith('S')]
    primers = []
    for id, key in zip(ids, keys_s):
        primer_file = root / 'codon_lists' / f'{key}.txt'
        with open(primer_file, 'r') as f:
            primer_list = [primer.strip() for primer in f.readlines()]
        primers += [primer_list[id]]
    selections = get_selections(config)
    p1 = (selections['FwdPrimer'] == primers[0])
    p2 = (selections['RevPrimer'] == primers[1])
    return selections[p1 & p2]['SelectionID'].squeeze()


def perform_selection(
        reads: list,
        num_reads: int,
        config_file: Path,
) -> dict:
    counts = defaultdict(lambda: defaultdict(int))
    for line in tqdm(reads, total=num_reads, ncols=100):
        _, *adapters = line.strip().split('?')
        primer_ids = [i.split('.')[-1] for i in filter(lambda x: 'S' in x, adapters)]
        primer_ids = list(map(int, primer_ids))
        selection_id = get_selection_id(primer_ids, config_file)
        barcodes = tuple(i.split('.')[-1] for i in filter(lambda x: 'B' in x, adapters))
        counts[selection_id][barcodes] += 1
    return counts


def save_counts(
        counts: dict,
        output_dir: Path,
        config_file: Path,
) -> None:
    for selection_id, count in tqdm(counts.items(), ncols=100):
        count = [(j, *i) for i, j in zip(count.keys(), count.values())]
        df = pd.DataFrame.from_records(count, columns=['Count', 'Code1', 'Code2'])
        df = df.astype(int)
        df.sort_values(['Code1', 'Code2'], inplace=True)
        selection_dir = output_dir / f'selection-{selection_id}'
        selection_dir.mkdir(parents=True, exist_ok=True)
        config = read_yaml(config_file)
        hash_value = hash_dict(config['Structure'])
        output_file = selection_dir / f'{hash_value}.txt'
        df.to_csv(output_file, index=False, sep='\t')
        shutil.copy(config_file, output_file.with_suffix('.yml'))


def compute_counts(
        config_file: Path,
        input_file: str,
        output_dir: str,
) -> None:
    input_dir = input_file.parent
    num_reads = json.load(open(
        sorted((input_dir).glob('*.cutadapt.json'))[-1]
    ))['read_counts']['output']
    with gzip.open(input_file, 'rt') as f:
        counts = perform_selection(f, num_reads, config_file)
    save_counts(counts, output_dir, config_file)

