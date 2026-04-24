from collections import defaultdict
import gzip
import multiprocessing
from pathlib import Path

import pandas as pd
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
    barcodes = tuple(int(i.split('.')[-1]) for i in filter(lambda x: 'B' in x, adapters))
    return {'selection_ids': selection_ids, 'barcodes': barcodes}


def _count_bytes_chunk(chunk: bytes) -> dict:
    """Count barcodes in a raw bytes chunk (newline-separated header lines).

    Mirrors extract_ids logic but operates on bytes to avoid decode overhead.
    """
    counts = {}
    for line in chunk.split(b'\n'):
        if not line:
            continue
        _, *adapters = line.split(b'?')
        selection_ids = tuple(int(a.split(b'.')[-1]) for a in adapters if b'S' in a)
        barcodes = tuple(int(a.split(b'.')[-1]) for a in adapters if b'B' in a)
        if selection_ids not in counts:
            counts[selection_ids] = {barcodes: 1}
        else:
            bc_dict = counts[selection_ids]
            bc_dict[barcodes] = bc_dict.get(barcodes, 0) + 1
    return counts


def _iter_byte_chunks(input_path: Path, chunk_size_bytes: int):
    """Yield raw byte chunks from a gzip file, always ending at a newline."""
    from isal import igzip
    with igzip.open(input_path, 'rb') as f:
        while True:
            chunk = f.read(chunk_size_bytes)
            if not chunk:
                break
            if chunk[-1:] != b'\n':
                chunk += f.readline()
            yield chunk


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
    codon_cols = [f'code_{i}' for i in range(num_codes)]
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


def get_counts(*, input_path: Path, num_reads: int, num_workers: int = 1,
               chunk_size_bytes: int = 5_000_000) -> dict:
    """Count barcode occurrences from a gzipped read file.

    Args:
        input_path: Path to the gzipped reads with adapter info.
        num_reads: Expected number of reads for progress tracking.
        num_workers: Number of worker processes (1 = serial, original behaviour).
        chunk_size_bytes: Uncompressed bytes per work unit in parallel mode.

    Returns:
        A nested dict of selection IDs to barcode counts.
    """
    if num_workers == 1:
        with gzip.open(input_path, 'rt') as f:
            counts = defaultdict(lambda: defaultdict(int))
            for line in tqdm(f, total=num_reads, ncols=100):
                ids = extract_ids(line)
                counts[ids['selection_ids']][ids['barcodes']] += 1
        return counts

    counts: dict = {}
    with multiprocessing.Pool(num_workers) as pool:
        for partial in tqdm(
            pool.imap_unordered(_count_bytes_chunk, _iter_byte_chunks(input_path, chunk_size_bytes)),
            ncols=100,
        ):
            for sel, bc_counts in partial.items():
                if sel not in counts:
                    counts[sel] = bc_counts
                else:
                    existing = counts[sel]
                    for bc, n in bc_counts.items():
                        existing[bc] = existing.get(bc, 0) + n

    return counts
