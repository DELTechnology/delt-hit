import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import stat
import textwrap

import pandas as pd
from pydantic import BaseModel, computed_field
import yaml


class Region(BaseModel):
    name: str
    codons: list[str]
    max_error_rate: float
    indels: int
    position_in_construct: int = None

    @computed_field
    @property
    def region_id(self) -> str:
        return f'{self.position_in_construct}-{self.name}'


def is_gz_file(
        path: Path,
) -> bool:
    with open(path, 'rb') as file:
        return file.read(2) == b'\x1f\x8b'


def read_yaml(
        path: Path,
) -> dict:
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def convert_dict(
        data: dict,
) -> tuple:
    if isinstance(data, dict):
        return tuple((key, convert_dict(value)) for key, value in data.items())
    elif isinstance(data, int):
        # TODO: fix return type to (float(data), )?
        return float(data)
    else:
        return data


def hash_dict(
        data: dict,
) -> tuple:
    data = convert_dict(data)
    data_str = json.dumps(data)
    hash_object = hashlib.sha256()
    hash_object.update(data_str.encode())
    return hash_object.hexdigest()


def get_selections(
        config: dict,
) -> pd.DataFrame:
    root = Path(config['Root'])
    config_selection = config['Selection']
    selection_file = root / config_selection['SelectionFile']
    fastq_file = str(Path(config_selection['FASTQFile']).name)
    library = str(Path(config_selection['Library']).name)
    selections = pd.read_excel(selection_file)
    return selections[(selections['FASTQFile'] == fastq_file) & (selections['Library'] == library)]


def convert_struct_file(
        struct_file: Path,
) -> None:
    with open(struct_file, 'r') as f:
        struct = f.readlines()[2:]
    struct = sorted(struct, key=lambda line: int(line.split('\t')[0]))
    config = {
        'Root': '~/',
        'Selection': {
            'SelectionFile': 'selection.xlsx',
            'FASTQFile': 'input.fastq.gz',
            'Library': 'library.xlsx'
        },
        'Structure': {},
    }
    max_error_rate = 0.0
    indels = 0
    indices = {}
    for line in struct:
        line = line.strip().split('\t')
        _type = line[2]
        indices[_type] = indices.get(_type, 0) + 1
        region = f'{_type}{indices[_type]}'
        config['Structure'][region] = {}
        config['Structure'][region]['MaxErrorRate'] = max_error_rate
        config['Structure'][region]['Indels'] = indels
    output_file = Path(struct_file).parent / 'config.yml'
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)


def get_regions(
        structure: dict,
) -> list[Region]:
    regions = []
    for key, value in structure.items():
        with open(value['Path'], 'r') as f:
            codons = f.read().split('\n')
            codons = filter(len, codons)
        region = Region(
            name=key,
            codons=codons,
            max_error_rate=value['MaxErrorRate'],
            indels=value['Indels']
        )
        regions.append(region)
    return regions


def write_fastq_files(
        regions: list[Region],
        path: str,
) -> None:
    fastq = ''
    for i, region in enumerate(regions):
        region.position_in_construct = i
        fastq = [f'>{region.region_id}.{index}\n{codon}'
                for index, codon in enumerate(region.codons)]
        fastq = '\n'.join(fastq)
        with open(path / f'{region.region_id}.fastq', 'w') as f:
            f.write(fastq)


def generate_input_files(
        config_file: Path,
        structure: dict,
        root_dir: Path,
        path_input_fastq: Path,
) -> None:
    
    path_input_dir = root_dir / 'cutadapt_input_files'
    path_input_dir.mkdir(parents=True, exist_ok=True)
    path_output_dir = root_dir / 'cutadapt_output_files'
    path_demultiplex_exec = path_input_dir / 'demultiplex.sh'
    path_final_reads = path_output_dir / 'reads_with_adapters.gz'
    path_output_fastq = path_output_dir / 'out.fastq.gz'
    path_counts = root_dir / 'evaluations'

    regions = get_regions(structure)
    write_fastq_files(regions, path_input_dir)

    with open(path_demultiplex_exec, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('# make sure you installed pigz with `brew install pigz` to enable parallel processing\n\n')
        f.write(f'mkdir "{path_output_dir}"\n')
        f.write(f'mkdir "{path_counts}"\n')
        # NOTE: we symlink the fastq file we want to demultiplex
        f.write(f'ln -sf "{path_input_fastq}" "{path_output_fastq}"\n')

    rename_command = '{id} {comment}?{adapter_name}'
    n_cores = multiprocessing.cpu_count()

    for region in regions:
        error_rate = region.max_error_rate
        indels = f' --no-indels' if not int(region.indels) else ''
        path_adapters = path_input_dir / f'{region.region_id}.fastq'
        # NOTE: from now on we use the output of the previous step as input
        path_input_fastq = path_output_dir / 'input.fastq.gz'

        report_file_name = path_output_dir / f'{region.region_id}.cutadapt.json'
        stdout_file_name = path_output_dir / f'{region.region_id}.cutadapt.log'
        info_file_name = path_output_dir / f'{region.region_id}.cutadapt.info.gz'

        with open(path_demultiplex_exec, 'a') as f:
            cmd = f"""
                mv "{path_output_fastq}" "{path_input_fastq}"
                
                cutadapt "{path_input_fastq}" \\
                -o "{path_output_fastq}" \\
                -e {error_rate}{indels} \\
                -g "^file:{path_adapters}" \\
                --rename '{rename_command}' \\
                --discard-untrimmed \\
                --json="{report_file_name}" \\
                --info-file="{info_file_name}" \\
                --cores={n_cores} 2>&1 | tee "{stdout_file_name}"
                
                """
            cmd = textwrap.dedent(cmd)
            f.write(cmd)

    with open(path_demultiplex_exec, 'a') as f:
        f.write(f'zgrep @ "{path_output_fastq}" | gzip -c > "{path_final_reads}"\n')
        f.write(f'delt-cli demultiplex compute-counts {config_file} "{path_final_reads}" "{path_counts}"\n')
        f.write(f'rm "{path_output_fastq}" "{path_input_fastq}"\n')

    os.chmod(path_demultiplex_exec, os.stat(path_demultiplex_exec).st_mode | stat.S_IEXEC)

