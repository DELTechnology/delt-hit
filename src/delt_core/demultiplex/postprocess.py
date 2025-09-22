from collections import defaultdict
import gzip
from pathlib import Path

import pandas as pd
from tqdm import tqdm

# from .preprocess import get_selections
# from ..utils import hash_dict
# from .validation import Config


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
        p1 = (selections['FwdPrimer1'] == fwd_primer)
        p2 = (selections['RevPrimer1'] == rev_primer)
        selection_id = selections[p1 & p2]['SelectionID'].squeeze()
        selection_ids.append(selection_id)
    return selection_ids


def extract_ids(line: str):
    _, *adapters = line.strip().split('?')
    selection_ids = [i.split('.')[-1] for i in filter(lambda x: 'S' in x, adapters)]
    selection_ids = tuple(map(int, selection_ids))
    barcodes = tuple(int(i.split('.')[-1]) + 1 for i in filter(lambda x: 'B' in x, adapters))
    return {'selection_ids': selection_ids, 'barcodes': barcodes}


def save_counts(counts: dict, output_dir: Path, ids_to_name: dict = None) -> None:

    num_codes = len(list(list(counts.values())[0].keys())[0])
    columns = [f'code_{i}' for i in range(1, num_codes + 1)] + ['count']

    for selection_ids, count in tqdm(counts.items(), ncols=100):
        count = [(*k, v) for k, v in count.items()]
        df = pd.DataFrame.from_records(count, columns=columns)
        df = df.astype(int)
        df.sort_values(columns[:2], inplace=True)
        # df.sort_values('count', inplace=True, ascending=False)

        if ids_to_name is None:
            name = '-'.join(map(str, selection_ids))
        else:
            name = ids_to_name[selection_ids]
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

def compute_counts(*, config: dict, input_path: Path, output_dir: Path, num_reads: int) -> None:
    counts = get_counts(input_path=input_path, num_reads=num_reads)

    list_of_selection_primer_ids = list(counts.keys())
    selection_ids = get_selection_ids(list_of_selection_primer_ids, config)
    counts = {selection_id: val for selection_id, val in zip(selection_ids, counts.values())}
    save_counts(counts, output_dir, config)

