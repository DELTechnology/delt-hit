from collections import OrderedDict
import os
from pathlib import Path
import shutil
import subprocess
import typing as tp

import numpy as np

from .utils import read_txt


def compute_counts(
        structure: tp.OrderedDict,
        input_file: str,
) -> tp.Dict:
    
    """
    Code    Count
    -------------
    GCCTCG  10
    TCCGAC  3
    CAAGTG  1
    """

    keys = list(structure.keys())
    num_reads = len(read_txt(input_file)) // 4
    indices = np.full((len(structure), num_reads), -1)
    
    map_sequences(structure, keys, indices, input_file)
    counts = create_table(structure, indices)

    return counts


def find_max_component(
        components: tp.Dict,
) -> str:
    return max(components, key=lambda k: len(components[k][0]))


def split_components(
        components: tp.OrderedDict,
        key: str,
) -> tp.Tuple[tp.OrderedDict, tp.OrderedDict]:
    split = list(components.keys()).index(key)
    comp1 = OrderedDict(list(components.items())[:split])
    comp2 = OrderedDict(list(components.items())[(split + 1):])
    return comp1, comp2


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
          f'-e {max_error_rate} -O {min_overlap} --trimmed-only --quiet'
    for sequence in sequences:
        cmd += f' -{mode} {sequence}'
    subprocess.run(cmd.split(), stdout=subprocess.DEVNULL)


def update_indices(
        indices: np.ndarray,
        info_file: str,
        output_file: str,
        comp_idx: int,
) -> None:
    
    info = read_txt(info_file)
    output = read_txt(output_file)
    
    for i, j in enumerate(range(1, len(output), 4)):
        read = info[i].split('\t')
        match = int(read[1]) >= 0
        
        if match:
            output[j] = read[6] + '\n'
            output[j+2] = read[10] + '\n'
            read_idx = int(read[0].split('-')[-1])
            seq_idx = int(read[7]) - 1
            indices[comp_idx, read_idx] = seq_idx

    update_output(output_file, output)


def update_output(
        output_file: str,
        output: str,
) -> None:
    copy = str(output_file).replace('a.fastq', 'b.fastq')
    shutil.copyfile(output_file, copy)
    with open(copy, 'w') as file:
        file.writelines(output)


def map_sequences(
        components: tp.Dict,
        keys: tp.List,
        indices: np.ndarray,
        input_file: str,
        max_error_rate: float = 0.1,
        min_overlap: int = 0.8,
) -> None:

    # Find the component with the longest sequence.
    max_component = find_max_component(components)

    parent = input_file.parent
    info_file = parent / f'info_{max_component}.tsv'
    output_file_a = parent / f'output_{max_component}a.fastq'
    output_file_b = Path(str(output_file_a).replace('a.fastq', 'b.fastq'))

    # Map the respective sequence.
    run_cutadapt(
        sequences=components[max_component],
        mode='a',
        input_file=input_file,
        info_file=info_file,
        output_file=output_file_a,
        max_error_rate=max_error_rate,
        min_overlap=min_overlap,
    )
    
    # Update the indices.
    comp_idx = keys.index(max_component)
    update_indices(indices, info_file, output_file_a, comp_idx)
    
    # Do the recursion.
    comp1, comp2 = split_components(components, max_component)
    if comp1:
        map_sequences(comp1, keys, indices, output_file_a)
    if comp2:
        map_sequences(comp2, keys, indices, output_file_b)

    os.remove(info_file)
    os.remove(output_file_a)
    os.remove(output_file_b)


def create_table(
        structure: tp.Dict,
        indices: np.ndarray,
):
    counts = {}
    for index in indices.T:
        if -1 in index:
            continue
        code = ''
        for i, sequences in zip(index, structure.values()):
            code += sequences[i]
        counts[code] = counts.get(code, 0) + 1
    return counts


