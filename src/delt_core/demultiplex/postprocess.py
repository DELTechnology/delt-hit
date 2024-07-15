from collections import defaultdict
import gzip
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from .preprocess import get_selections
from .utils import write_yaml, hash_dict, Config


def get_selection_ids(
        list_of_selection_primer_ids: list[list],
        config: dict,
) -> list[int]:
    root = Path(config['Root'])
    experiment_name = config['Experiment']['Name']

    # NOTE: this requires to respect the order
    structure_keys = filter(lambda x: x.startswith('S'), config['Structure'].keys())
    primer_lists = []
    for key in structure_keys:
        primer_file = root / 'experiments' / experiment_name / 'codon_lists' / f'{key}.txt'
        with open(primer_file, 'r') as f:
            primer_lists.append([primer.strip() for primer in f.readlines()])

    selections = get_selections(config)
    selection_ids = []
    for selection_primer_ids in list_of_selection_primer_ids:
        # NOTE: we would need to specify the column_name for the primers to avoid implicitly assume the order
        fwd_primer = primer_lists[0][selection_primer_ids[0]]
        rev_primer = primer_lists[1][selection_primer_ids[1]]
        p1 = (selections['FwdPrimer'] == fwd_primer)
        p2 = (selections['RevPrimer'] == rev_primer)
        selection_id = selections[p1 & p2]['SelectionID'].squeeze()
        selection_ids.append(selection_id)
    return selection_ids


def extract_ids(line: str):
    _, *adapters = line.strip().split('?')
    selection_ids = [i.split('.')[-1] for i in filter(lambda x: 'S' in x, adapters)]
    selection_ids = tuple(map(int, selection_ids))
    barcodes = tuple(i.split('.')[-1] for i in filter(lambda x: 'B' in x, adapters))
    return {'selection_ids': selection_ids, 'barcodes': barcodes}


def save_counts(
        counts: dict,
        output_dir: Path,
        config: dict,
) -> None:
    if not counts:
        return
    hash_value = hash_dict(config['Structure'])
    num_codes = len(list(list(counts.values())[0].keys())[0])
    columns = [f'Code{i}' for i in range(1, num_codes + 1)]
    columns.insert(0, 'Count')
    
    for selection_id, count in tqdm(counts.items(), ncols=100):
        count = [(j, *i) for i, j in zip(count.keys(), count.values())]
        df = pd.DataFrame.from_records(count, columns=columns)
        df = df.astype(int)
        df.sort_values(columns[1:], inplace=True)
        selection_dir = output_dir / f'selection-{selection_id}'
        selection_dir.mkdir(parents=True, exist_ok=True)
        output_file = selection_dir / f'{hash_value}.txt'
        df.to_csv(output_file, index=False, sep='\t')
        Config(**config).write_yaml(output_file)


def compute_counts(
        *,
        config: dict,
        input_file: Path,
        output_dir: Path,
        num_reads: int,
) -> None:
    with gzip.open(input_file, 'rt') as f:
        counts = defaultdict(lambda: defaultdict(int))
        for line in tqdm(f, total=num_reads, ncols=100):
            ids = extract_ids(line)
            counts[ids['selection_ids']][ids['barcodes']] += 1

    list_of_selection_primer_ids = list(counts.keys())
    selection_ids = get_selection_ids(list_of_selection_primer_ids, config)
    counts = {selection_id: val for selection_id, val in zip(selection_ids, counts.values())}
    save_counts(counts, output_dir, config)

