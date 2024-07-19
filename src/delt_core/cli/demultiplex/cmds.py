import json
from pathlib import Path
import subprocess

from ... import compute as c
from ... import demultiplex as d
from delt_core.demultiplex.utils import hash_dict, is_gz_file, init_config, Config


def init(
        root: Path,
        experiment_name: str,
        selection_file: Path,
        fastq_file: Path,
        library: Path,
        simulation: dict = None,
) -> None:
    if not root:
        root = Path.cwd()
    if Path(library).exists():
        bbs, _, _, _ = c.load_data(library)
        structure = ['S1']
        for i in range(1, len(bbs) + 1):
            structure += [f'C{i}']
            structure += [f'B{i}']
        structure += [f'C{len(bbs) + 1}', 'S2']
    else:
        structure = ['S1', 'C1', 'B1', 'C2', 'B2', 'C3', 'S2']
    init_config(
        structure=structure,
        root=root,
        experiment_name=experiment_name,
        selection_file=selection_file,
        fastq_file=fastq_file,
        library=library,
        simulation=simulation,
    )


def convert(
        struct_file: Path,
) -> None:
    with open(struct_file, 'r') as f:
        lines = f.readlines()[2:]
    lines = [line.strip().split() for line in lines]
    lines = sorted(filter(None, lines), key=lambda x: int(x[0]))
    indices = {}
    structure = []
    for line in lines:
        _type = line[2]
        indices[_type] = indices.get(_type, 0) + 1
        structure += [f'{_type}{indices[_type]}']
    init_config(
        structure=structure,
    )


def create_lists(
        config_file: Path,
        selection_id: int = None,
        output_dir: Path = None,
) -> dict:
    config_file = Path(config_file).resolve()
    config = Config.from_yaml(config_file).model_dump()
    root = config['Root']
    structure = config['Structure']

    selections = d.get_selections(config, selection_id)

    # WARNING: We cannot do this, this alters the content of the S1/2 lists and thus the indices of the primers
    #   and leads to mappings to the wrong selection ids
    hash_value = hash_dict(structure)
    # for selection_id in selections['SelectionID']:
    #     path = root / 'evaluations' / f'selection-{selection_id}' / f'{hash_value}.txt'
    #     if path.exists():
    #         selections = selections[selections['SelectionID'] != selection_id]
    # if selections.empty:
    #     exit()

    keys = list(structure.keys())
    lib_file = root / config['Selection']['Library']
    bbs, _, _, consts = c.load_data(lib_file)
    if not output_dir:
        output_dir = config_file.parent / 'codon_lists'
    Path(output_dir).mkdir(exist_ok=True)

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
    consts = list(filter(None, sequence.split('{codon}')))
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
        *,
        config_file: Path,
        selection_id: int = None,
        write_json_file: bool = True,
        write_info_file: bool = True,
        fast_dev_run: bool = False,
) -> None:

    structure = create_lists(config_file, selection_id)
    config = Config.from_yaml(config_file).model_dump()
    root_dir = config['Root']
    path_input_fastq = root_dir / config['Selection']['FASTQFile']
    if not is_gz_file(path_input_fastq):
        subprocess.run(['gzip', path_input_fastq])
        path_input_fastq = path_input_fastq.parent / (path_input_fastq.name + '.gz')
    d.generate_input_files(config_file=config_file, structure=structure, root_dir=root_dir,
                           path_input_fastq=path_input_fastq,
                           write_json_file=write_json_file, write_info_file=write_info_file, fast_dev_run=fast_dev_run)


def compute_counts(
        config_file: Path,
        input_file: Path,
        output_dir: Path,
) -> None:
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = input_file.parent
    config = Config.from_yaml(config_file).model_dump()
    num_reads = json.load(open(
        sorted(input_dir.glob('*.cutadapt.json'))[-1]
    ))['read_counts']['output']
    d.compute_counts(config=config, input_file=input_file, num_reads=num_reads, output_dir=output_dir)


def run(
        *,
        config_file: Path,
        selection_id: int = None,
        write_json_file: bool = True,
        write_info_file: bool = False,
        fast_dev_run: bool = False,
) -> None:
    create_cutadapt_input(config_file=config_file, selection_id=selection_id,
                          write_json_file=write_json_file, write_info_file=write_info_file, fast_dev_run=fast_dev_run)
    config = Config.from_yaml(config_file).model_dump()
    root = config['Root']
    experiment_name = config['Experiment']['Name']
    input_file = root / 'experiments' / experiment_name / 'cutadapt_input_files' / 'demultiplex.sh'
    subprocess.run(['bash', input_file])

