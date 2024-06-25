from pathlib import Path
import subprocess

import yaml

from ... import compute as c
from ... import demultiplex as d


def init(
        root = None,
        selection_file: Path = 'selection.xlsx',
        fastq_file: Path = 'input.fastq.gz',
        library: Path = 'library.xlsx',
) -> None:
    if not root:
        root = Path.cwd()
    config = {
        'Root': str(root),
        'Selection': {
            'SelectionFile': selection_file,
            'FASTQFile': fastq_file,
            'Library': library,
        },
        'Structure': {},
    }
    max_error_rate = 0.0
    indels = 0
    structure = ['S1', 'C1', 'B1', 'C2', 'B2', 'C3', 'S2']
    for region in structure:
        config['Structure'][region] = {}
        config['Structure'][region]['MaxErrorRate'] = max_error_rate
        config['Structure'][region]['Indels'] = indels
    output_file = Path(root) / 'config.yml'
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)


def convert(
        struct_file: Path,
) -> None:
    d.convert_struct_file(struct_file)


def create_lists(
        config_file: Path,
        output_dir: Path = None,
) -> dict:
    config_file = Path(config_file).resolve()
    if not output_dir:
        output_dir = config_file.parent / 'codon_lists'
    Path(output_dir).mkdir(exist_ok=True)
    config = d.read_yaml(config_file)
    selections = d.get_selections(config)
    structure = config['Structure']
    keys = list(structure.keys())
    root = Path(config['Root'])
    lib_file = root / 'libraries' / config['Selection']['Library']
    bbs, _, _, consts = c.load_data(lib_file)
    # Building blocks.
    keys_b = [key for key in keys if key.startswith('B')]
    assert len(bbs) == len(keys_b)
    for bb in bbs:
        codes = bb['Codon']
        key = keys_b.pop(0)
        output_file = output_dir / f'{key}.txt'
        structure[key]['Path'] = output_file
        with open(output_file, 'w') as f:
            for code in codes:
                f.write(code)
                f.write('\n')
    # Constant regions.
    keys_c = [key for key in keys if key.startswith('C')]
    sequence = consts['Sequence'].squeeze()
    consts = sequence.split('{codon}')
    assert len(consts) == len(keys_c)
    for const in consts:
        key = keys_c.pop(0)
        output_file = output_dir / f'{key}.txt'
        structure[key]['Path'] = output_file
        with open(output_file, 'w') as f:
            f.write(const)
            f.write('\n')
    # Primers.
    keys_s = [key for key in keys if key.startswith('S')]
    primer_lists = [selections['FwdPrimer'], selections['RevPrimer']]
    assert len(primer_lists) == len(keys_s)
    for primer_list in primer_lists:
        key = keys_s.pop(0)
        output_file = output_dir / f'{key}.txt'
        structure[key]['Path'] = output_file
        with open(output_file, 'w') as f:
            for primer in primer_list.unique():
                f.write(primer)
                f.write('\n')
    return structure


def create_cutadapt_input(
        config_file: Path,
) -> None:
    structure = create_lists(config_file)
    config = d.read_yaml(config_file)
    root = Path(config['Root'])
    fastq_file = root / 'fastq_files' / config['Selection']['FASTQFile']
    if not d.is_gz_file(fastq_file):
        subprocess.run(['gzip', fastq_file])
        fastq_file = fastq_file.parent / (fastq_file.name + '.gz')
    d.generate_input_files(config_file, structure, root, fastq_file)


def compute_counts(
        config_file: Path,
        input_file: Path,
        output_dir: Path,
) -> None:
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    d.compute_counts(config_file, input_file, output_dir)


def run(
        config_file: Path,
) -> None:
    create_cutadapt_input(config_file)
    config = d.read_yaml(config_file)
    root = Path(config['Root'])
    input_file = root / 'cutadapt_input_files' / 'demultiplex.sh'
    subprocess.run(['bash', input_file])

