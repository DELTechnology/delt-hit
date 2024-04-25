import json
from pathlib import Path
import subprocess


def read_json(
        path: str,
):
    with open(path, 'r') as file:
        return json.load(file)


def read_tsv(
        path: str,
):
    with open(path, 'r') as file:
        return file.readlines()


def compute_counts(
        components: dict,
        config: dict,
):
    """
    Code    Count
    -------------
    GCCTCG  10
    TCCGAC  3
    CAAGTG  1
    """

    
    fastq_file = config['fastq_file']
    info_file = config['info_file']

    b1 = components['B1']

    counts = {}

    for seq in b1:
        counts[seq] = 0

        if seq == 'ACACTC':
            break

        cmd = f'cutadapt -a {seq} -o /dev/null {fastq_file} --info-file={info_file}'
        subprocess.run(cmd.split(), stdout=subprocess.DEVNULL)
        reads = read_tsv(info_file)

        for read in reads:
            if int(read.split()[1]) >= 0 and len(read.split()[5]) >= 6:
                counts[seq] += 1
                print(seq, read.split())

    return counts



if __name__ == '__main__':
    
    path = Path('/Users/Gary/Downloads/zivi/git/delt-core/data/simulation')
    components = read_json(path / 'structure/components.json')
    config = read_json(path / 'config/config_cutadapt.json')

    counts = compute_counts(components, config)

    print(counts)
    print(sum(counts.values()))


    # with open(os.path.join(path, 'components.json'), 'w') as file:
    #     json.dump(components, file)


