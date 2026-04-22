from collections import defaultdict
import multiprocessing
from pathlib import Path

import pandas as pd
from isal import igzip
from tqdm import tqdm

def extract_ids(line: str):
    """Extract selection and barcode IDs from a cutadapt info line.

    Args:
        line: A line from the cutadapt info file.

    Returns:
        A dict with selection ID tuples and barcode tuples.
    """
    _, *adapters = line.strip().split('?')
    selection_ids = [i.split('.')[-1] for i in filter(lambda x: 'S' in x, adapters)]
    selection_ids = tuple(map(int, selection_ids))
    barcodes = tuple(int(i.split('.')[-1]) + 1 for i in filter(lambda x: 'B' in x, adapters))
    return {'selection_ids': selection_ids, 'barcodes': barcodes}


def save_counts(counts: dict, output_dir: Path, ids_to_name: dict = None,
                as_files: bool = True, sort_by_counts: bool = True) -> None:
    """Persist count tables to disk.

    Args:
        counts: Nested dict of selection IDs to barcode counts.
        output_dir: Directory to write output files.
        ids_to_name: Optional mapping from selection ID tuples to names.
        as_files: Whether to store counts as flat files or nested dirs.
        sort_by_counts: Whether to sort descending by count.
    """

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


def _count_batch(lines: list[str]) -> dict[tuple, dict[tuple, int]]:
    counts: dict[tuple, dict[tuple, int]] = {}
    for line in lines:
        ids = extract_ids(line)
        sel = ids['selection_ids']
        bc = ids['barcodes']
        if sel not in counts:
            counts[sel] = {}
        counts[sel][bc] = counts[sel].get(bc, 0) + 1
    return counts


def _merge(counts: defaultdict, partial: dict) -> None:
    for sel, barcodes in partial.items():
        inner = counts[sel]
        for bc, n in barcodes.items():
            inner[bc] = inner.get(bc, 0) + n


def get_counts(*, input_path: Path, num_reads: int,
               num_workers: int | None = None,
               batch_size: int = 50_000) -> dict:
    """Count barcode occurrences from a bgzipped read file in parallel.

    Args:
        input_path: Path to the bgzipped reads with adapter info.
        num_reads: Expected number of reads for progress tracking.
        num_workers: Number of worker processes (defaults to CPU count).
        batch_size: Lines per worker batch.

    Returns:
        A nested dict of selection IDs to barcode counts.
    """
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()

    counts: defaultdict = defaultdict(dict)

    def _batches(f):
        batch = []
        for line in f:
            batch.append(line)
            if len(batch) == batch_size:
                yield batch
                batch = []
        if batch:
            yield batch

    with igzip.open(input_path, 'rt') as f, \
            multiprocessing.Pool(num_workers) as pool:
        for partial in tqdm(
            pool.imap(_count_batch, _batches(f)),
            total=-(-num_reads // batch_size),
            ncols=100,
        ):
            _merge(counts, partial)

    return counts

