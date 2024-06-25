import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from delt_core.cli.demultiplex.cmds import init


def read_txt(
        path: Path,
) -> None:
    with open(path, 'r') as file:
        return file.readlines()


def write_txt(
        data: list,
        path: Path,
) -> None:
    if Path(path).suffix == '.gz':
        with gzip.open(path, 'wt') as f:
            f.writelines(data)
    else:
        with open(path, 'w') as f:
            f.writelines(data)


def read_yaml(
        path: Path,
) -> None:
    with open(path, 'r') as file:
        return yaml.safe_load(file)


def write_yaml(
        data: dict,
        path: Path,
) -> None:
    with open(path, 'w') as file:
        yaml.dump(data, file, default_flow_style=False, sort_keys=False)


def create_config_file(
        root: Path,
        config_file: Path,
        selection_file: Path,
        library: Path,
        fastq_file: Path,
        output_file: Path,
        num_reads: int,
) -> None:
    if not root:
        root = Path.cwd()
    init(root, config_file, selection_file, fastq_file, library)
    config = read_yaml(Path(root) / config_file)
    config['Simulation'] = {}
    config['Simulation']['OutputFile'] = output_file
    config['Simulation']['NumReads'] = num_reads
    config['Simulation']['Errors'] = []
    write_yaml(config, Path(root) / config_file)


def generate_codons(
        num_codons: int = 100,
        length: int = 10,
) -> str:
    bases = ['A', 'T', 'C', 'G']
    rng = np.random.default_rng()
    codons = []
    for _ in range(num_codons):
        codons += [''.join([rng.choice(bases) for _ in range(length)])]
    return codons


def create_library_template(
        root: Path,
        library: Path,
) -> Path:
    if not root:
        root = Path.cwd()
    library = str(Path(root) / library)
    num_entries = 100
    step1 = {
        'ID': np.arange(1, num_entries + 1),
        'SMILES': '',
        'ScaffoldID': '',
        'Codon': generate_codons(num_entries),
        'ReactionType': '',
    }
    step2 = {
        'ID': np.arange(1, num_entries + 1),
        'SMILES': '',
        'Codon': generate_codons(num_entries),
        'ReactionType': '',
    }
    scaffolds = {
        'ScaffoldID': [],
        'SMILES': [],
    }
    smarts = {
        'ReactionType': [],
        'SMARTS': [],
    }
    const = {
        'Sequence': ['{codon}'.join(generate_codons(3))],
        'Reverse': [0],
        'Complement': [0],
    }
    sheets = [step1, step2, scaffolds, smarts, const]
    names = ['step1', 'step2', 'scaffolds', 'smarts', 'const']
    with pd.ExcelWriter(library) as writer:
        for sheet, name in zip(sheets, names):
            pd.DataFrame(sheet).to_excel(writer, sheet_name=name, index=False)
    return library


def create_selection_template(
        root: Path,
        selection_file: Path,
        library: Path,
        fastq_file: Path,
) -> Path:
    if not root:
        root = Path.cwd()
    selection_file = str(Path(root) / selection_file)
    num_fwd_primers = 10
    num_rev_primers = 10
    num_entries = num_fwd_primers * num_rev_primers
    selection = {
        'SelectionID': np.arange(1, num_entries + 1),
        'Library': num_entries * [Path(library).name],
        'FwdPrimer': num_rev_primers * [*generate_codons(num_fwd_primers)],
        'RevPrimer': [codon for codon in generate_codons(num_rev_primers) for _ in range(num_fwd_primers)],
        'FASTQFile': num_entries * [Path(fastq_file).name],
    }
    with pd.ExcelWriter(selection_file) as writer:
        pd.DataFrame(selection).to_excel(writer, sheet_name='selections', index=False)
    return selection_file

