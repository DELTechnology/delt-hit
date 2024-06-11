from pathlib import Path
import subprocess

import pandas as pd

from ... import compute as c
from ... import demultiplex as d


def convert(
        struct_file: Path,
) -> None:
    d.convert_struct_file(struct_file)


def init(
        excel_file: Path,
        struct_file: Path,
        fastq_file: Path,
) -> None:
    output_dir = Path(struct_file).parent / 'sequences'
    create_lists(excel_file, output_dir)
    create_cutadapt_input(struct_file, fastq_file)


def create_lists(
        input_file: Path,
        output_dir: Path = None,
) -> None:
    # TODO: Add primer sequences.
    if not output_dir:
        output_dir = Path(input_file).parent / 'sequences'
    Path(output_dir).mkdir(exist_ok=True)
    bbs, _, _, consts = c.load_data(input_file)
    for i, bb in enumerate(bbs, 1):
        codes = bb['Codon']
        with open(output_dir / f'b{i}.txt', 'w') as f:
            for code in codes:
                f.write(code)
                f.write('\n')
    sequence = consts['Sequence'][0]
    consts = sequence.split('{codon}')
    for i, const in enumerate(consts, 1):
        with open(output_dir / f'c{i}.txt', 'w') as f:
            f.write(const)
            f.write('\n')


def create_cutadapt_input(
        struct_file: Path,
        fastq_file: Path = None,
) -> None:
    struct_file = Path(struct_file).resolve()
    if not fastq_file:
        fastq_file = struct_file.parent / 'input.fastq.gz'
    fastq_file = Path(fastq_file).resolve()
    d.generate_input_files(struct_file, fastq_file)


def compute_counts(
        input_file: Path,
        output_dir: Path,
) -> None:
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    d.compute_counts(input_file, output_dir)


def run(
        struct_file: Path,
        fastq_file: Path = None,
) -> None:
    struct_file = Path(struct_file).resolve()
    create_cutadapt_input(struct_file, fastq_file)
    input_file = struct_file.parent / 'cutadapt_input_files' / 'demultiplex.sh'
    subprocess.run(['bash', input_file])

