import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from utils import read_json, read_txt, remove_reads


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

    structure = read_json(struct_file)

    indices = map_sequences(components, structure, fastq_file, output_file)
    print(indices)

    counts = {}
    return counts


def map_sequences(
        components: dict,
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

    for comp_idx, (component, seqs) in enumerate(components.items()):
        
        for seq_idx, seq in enumerate(seqs):

            if seq_idx > 10:
                break

            cmd = f'cutadapt -g {seq} -o {output_file} {tmp} --info-file={info_file} -e {max_error_rate}'
            subprocess.run(cmd.split(), stdout=subprocess.DEVNULL)

            reads = read_txt(info_file)

            for read in reads:
                read = read.split()
                read_idx = int(read[0].split('-')[-1])
                match = int(read[1]) >= 0
                
                if match:
                    start = int(read[2])
                    end = int(read[3])
                    error_len = abs(end - start - len(seq))
                    error_dist = abs(start - structure[component])
                    
                    if error_len <= epsilon_len and error_dist <= epsilon_dist:
                        indices[comp_idx, read_idx] = seq_idx
                        # print(seq, read_idx, seq_idx, read)
        
        # tmp = remove_reads(tmp, indices[comp_idx])
        # shutil.copyfile(output_file, tmp)
        # if comp_idx == 1:
        #     break
    
    os.remove(tmp)
    return indices


if __name__ == '__main__':
    
    path = Path('/Users/Gary/Downloads/zivi/git/delt-core/data/simulation')
    components = read_json(path / 'structure/components.json')
    config = read_json(path / 'config/config_cutadapt.json')

    counts = compute_counts(components, config)

    print(sum(counts.values()))

