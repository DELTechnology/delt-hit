from pathlib import Path
import shutil
import subprocess

from ... import compute as c
from ... import demultiplex as d


def convert(
        struct_file: Path,
) -> None:
    d.convert_struct_file(struct_file)


def create_lists(
        config_file: Path,
        output_dir: Path = None,
) -> None:
    if not output_dir:
        output_dir = Path.cwd() / 'sequences'
    Path(output_dir).mkdir(exist_ok=True)
    config_file = Path(config_file).resolve()
    config = d.read_yaml(config_file)
    structure = config['Structure']
    keys = list(structure.keys())
    selection = d.get_selection(config)
    lib_file = selection['Library'][0]
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
    sequence = consts['Sequence'][0]
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
    selections = [selection['FwdPrimer1'][0], selection['RevPrimer1'][0]]
    assert len(selections) == len(keys_s)
    for s in selections:
        key = keys_s.pop(0)
        output_file = output_dir / f'{key}.txt'
        structure[key]['Path'] = output_file
        shutil.copy(s, output_file)
    return structure


def create_cutadapt_input(
        config_file: Path,
) -> None:
    structure = create_lists(config_file)
    config = d.read_yaml(config_file)
    selection = d.get_selection(config)
    fastq_file = selection['FASTQFile'][0]
    fastq_file = Path(fastq_file).resolve()
    if not d.is_gz_file(fastq_file):
        subprocess.run(['gzip', fastq_file])
        fastq_file = fastq_file.parent / (fastq_file.name + '.gz')
    d.generate_input_files(structure, fastq_file)


def compute_counts(
        input_file: Path,
        output_dir: Path,
) -> None:
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    d.compute_counts(input_file, output_dir)


def run(
        config_file: Path,
) -> None:
    create_cutadapt_input(config_file)
    config = d.read_yaml(config_file)
    selection = d.get_selection(config)
    fastq_file = selection['FASTQFile'][0]
    input_file = Path(fastq_file).parent / 'cutadapt_input_files' / 'demultiplex.sh'
    subprocess.run(['bash', input_file])

