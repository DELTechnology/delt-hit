import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from utils import read_json, read_txt


def compute_counts(
        components: dict,
        const: dict,
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

    structure = read_json(struct_file)

    indices = map_sequences(components, consts, structure, fastq_file, output_file)
    print(indices)

    counts = {}
    return counts


def map_sequences(
        components: dict,
        consts: dict,
        structure: str,
        fastq_file: str,
        output_file: str,
        max_error_rate: float = .0,
        epsilon_len: int = 0,
        epsilon_dist: int = 0,
) -> np.ndarray:
    
    parent = Path(output_file).parent
    info_file = parent / 'info.tsv'
    tmp = parent / 'tmp.fastq'
    shutil.copyfile(fastq_file, tmp)

    num_reads = len(read_txt(fastq_file)) // 4
    indices = np.full((len(components), num_reads), -1)

    # for const_idx, (const, seqs) in enumerate(consts.items()):

    prev_comps = []

    for comp_idx, (component, seqs) in enumerate(components.items()):
        if component[0] != 'C':
            prev_comps += [component]
            continue
        
        const = seqs[0]

        # 3'
        min_overlap = len(seqs[0]) - epsilon_len
        cmd = f'cutadapt -o {output_file} {tmp} --info-file={info_file} ' \
              f'-e {max_error_rate} -O {min_overlap} --trimmed-only -a {const}'
        subprocess.run(cmd.split(), stdout=subprocess.DEVNULL)
        
        for prev_comp in prev_comps:
            seqs = components[prev_comp]
            min_overlap = len(seqs[0]) - epsilon_len
            cmd = f'cutadapt -o /dev/null {output_file} --info-file={info_file} ' \
                f'-e {max_error_rate} -O {min_overlap} --trimmed-only'
            for seq in seqs:
                cmd += f' -g {seq}'
            subprocess.run(cmd.split(), stdout=subprocess.DEVNULL)

            reads = read_txt(info_file)

            for read in reads:
                read = read.split('\t')
                read_idx = int(read[0].split('-')[-1])
                match = int(read[1]) >= 0

                if match:
                    seq_idx = int(read[7]) - 1
                    indices[comp_idx - 1, read_idx] = seq_idx
        

        # 5'
        min_overlap = len(seqs[0]) - epsilon_len
        cmd = f'cutadapt -o {output_file} {tmp} --info-file={info_file} ' \
              f'-e {max_error_rate} -O {min_overlap} --trimmed-only -g {const}' # --quiet

        # for seq in seqs:
                # cmd += f' -g X{seq}' # non-internal 5’ adapter
                # cmd += f' -g ^{seq}' # anchor 5’ adapter
                # cmd += f' -g {seq}' # regular 5’ adapter

        subprocess.run(cmd.split(), stdout=subprocess.DEVNULL)

        reads = read_txt(info_file)

        for read in reads:
            read = read.split('\t')
            read_idx = int(read[0].split('-')[-1])
            match = int(read[1]) >= 0

            if match:
                seq_idx = int(read[7]) - 1
                indices[comp_idx, read_idx] = seq_idx
                continue

                start = int(read[2])
                end = int(read[3])
                error_len = abs(end - start - len(seq))
                error_dist = abs(start - structure[component])
                
                if error_len <= epsilon_len and error_dist <= epsilon_dist:
                    indices[comp_idx, read_idx] = seq_idx
                    # print(seq, read_idx, seq_idx, read)
        
        shutil.copyfile(output_file, tmp)
    
    os.remove(tmp)
    return indices


if __name__ == '__main__':
    
    path = Path('/Users/Gary/Downloads/zivi/git/delt-core/data/simulation')
    components = read_json(path / 'structure/components.json')
    consts = read_json(path / 'structure/consts.json')
    config = read_json(path / 'config/config_cutadapt.json')

    counts = compute_counts(components, consts, config)

    # print(sum(counts.values()))

