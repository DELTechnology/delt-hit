import os
from pathlib import Path
import shutil
import subprocess
import typing as tp

import numpy as np

from utils import read_json, read_txt


def compute_counts(
        components: dict,
        config: dict,
) -> None:
    """
    Code    Count
    -------------
    GCCTCG  10
    TCCGAC  3
    CAAGTG  1
    """
    
    struct_file = config['struct_file']
    fastq_file = config['fastq_file']
    output_file = config['output_file']

    indices = map_sequences(components, fastq_file, output_file)
    print(indices)

    counts = {}

    return counts


def run_cutadapt(
        sequences: tp.List[str],
        mode: str,
        input_file: str,
        info_file: str,
        output_file: str = '/dev/null',
        max_error_rate: float = 0.1,
        min_overlap: float = 0.5,
) -> None:
    """
    Mode:
    -a: 3' adapter
    -g: 5' adapter
    """
    min_overlap = round(len(sequences[0]) * min_overlap)
    cmd = f'cutadapt -o {output_file} {input_file} --info-file={info_file} ' \
          f'-e {max_error_rate} -O {min_overlap} --trimmed-only' # --quiet
    for sequence in sequences:
        cmd += f' -{mode} {sequence}'
    subprocess.run(cmd.split(), stdout=subprocess.DEVNULL)


def update_indices(
        indices: np.ndarray,
        info_file: str,
        comp_idx: int,
) -> None:
    
    reads = read_txt(info_file)
    
    for read in reads:
        read = read.split('\t')
        match = int(read[1]) >= 0
        
        if match:
            read_idx = int(read[0].split('-')[-1])
            seq_idx = int(read[7]) - 1
            indices[comp_idx, read_idx] = seq_idx


def update_output():
    pass


def map_sequences(
        components: dict,
        fastq_file: str,
        output_file: str,
        max_error_rate: float = 0.1,
        min_overlap: int = 0.8,
) -> np.ndarray:
    
    parent = Path(output_file).parent
    info_file = parent / 'info.tsv'
    tmp = parent / 'tmp.fastq'
    shutil.copyfile(fastq_file, tmp)

    num_reads = len(read_txt(fastq_file)) // 4
    indices = np.full((len(components), num_reads), -1)
    pre_comps = []

    for comp_idx, (component, seqs) in enumerate(components.items()):

        # Map the rightmost region.
        if comp_idx == len(components) - 1:
            run_cutadapt(
                sequences=seqs,
                mode='g',
                input_file=tmp,
                info_file=info_file,
                max_error_rate=max_error_rate,
                min_overlap=min_overlap,
            )
            update_indices(indices, info_file, comp_idx)

        # Map the leftmost constant region of the remaining sequence.
        if component[0] != 'C':
            pre_comps += [component]
            continue
        
        run_cutadapt(
            sequences=seqs,
            mode='a',
            input_file=tmp,
            info_file=info_file,
            output_file=output_file,
            max_error_rate=max_error_rate,
            min_overlap=min_overlap,
        )
        update_indices(indices, info_file, comp_idx)
        
        # Map the sequence preceding the constant region.
        for i, pre_comp in enumerate(pre_comps):
            run_cutadapt(
                sequences=components[pre_comp],
                mode='g',
                input_file=output_file,
                info_file=info_file,
                max_error_rate=max_error_rate,
                min_overlap=min_overlap,
            )
            update_indices(indices, info_file, comp_idx - (len(pre_comps) - i))
        pre_comps = []

        # Keep the sequence following the constant region.
        run_cutadapt(
            sequences=seqs,
            mode='g',
            input_file=tmp,
            info_file=info_file,
            output_file=output_file,
            max_error_rate=max_error_rate,
            min_overlap=min_overlap,
        )
        shutil.copyfile(output_file, tmp)

    os.remove(tmp)
    return indices


if __name__ == '__main__':
    
    path = Path('/Users/Gary/Downloads/zivi/git/delt-core/data/simulation')
    components = read_json(path / 'structure/components.json')
    config = read_json(path / 'config/config_cutadapt.json')

    counts = compute_counts(components, config)

    # print(sum(counts.values()))

