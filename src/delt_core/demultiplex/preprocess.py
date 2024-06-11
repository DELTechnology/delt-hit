import multiprocessing
import os
from pathlib import Path
import stat
import textwrap

import pandas as pd
from pydantic import BaseModel, computed_field


class Region(BaseModel):
    name: str
    codons: list[str]
    max_error_rate: float
    indels: bool
    position_in_construct: int = None

    @computed_field
    @property
    def region_id(self) -> str:
        return f'{self.position_in_construct}-{self.name}'

    
def read_struct(
        path: str,
) -> dict:
    data = pd.read_excel(path).to_dict()
    struct = {}
    for key, value in data['Region'].items():
        # TODO: create data model for struct
        struct[value] = {
            'path': data['Path'][key],
            'max_error_rate': data['MaxErrorRate'][key],
            'indels': data['Indels'][key]
            # TODO: indel flag
        }
    return struct


def get_regions(
        struct_file: Path,
) -> list[Region]:
    struct = read_struct(struct_file)
    regions = []
    for key, value in struct.items():
        with open(struct_file.parent / value['path'], 'r') as f:
            codons = f.read().split('\n')
            codons = filter(len, codons)
        region = Region(
            name=key,
            codons=codons,
            max_error_rate=value['max_error_rate'],
            indels=value['indels']
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
        path_struct_file: Path,
        path_input_fastq: Path,
) -> None:
    dir = path_struct_file.parent
    path_input_dir = dir / 'cutadapt_input_files'
    path_input_dir.mkdir(parents=True, exist_ok=True)
    path_output_dir = dir / 'cutadapt_output_files'
    path_demultiplex_exec = path_input_dir / 'demultiplex.sh'
    path_final_reads = path_output_dir / 'reads_with_adapters.gz'
    path_output_fastq = path_output_dir / 'out.fastq.gz'
    path_counts = dir / 'counts'

    regions = get_regions(path_struct_file)
    write_fastq_files(regions, path_input_dir)

    with open(path_demultiplex_exec, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('# make sure you installed pigz with `brew install pigz` to enable parallel processing\n\n')
        f.write(f'mkdir {path_output_dir}\n')
        f.write(f'mkdir {path_counts}\n')
        # NOTE: we symlink the fastq file we want to demultiplex
        f.write(f'ln -sf {path_input_fastq} {path_output_fastq}\n')

    rename_command = '{id} {comment}?{adapter_name}'
    n_cores = multiprocessing.cpu_count()

    for region in regions:
        error_rate = region.max_error_rate
        path_adapters = path_input_dir / f'{region.region_id}.fastq'
        # NOTE: from now on we use the output of the previous step as input
        path_input_fastq = path_output_dir / 'input.fastq.gz'

        report_file_name = path_output_dir / f'{region.region_id}.cutadapt.json'
        stdout_file_name = path_output_dir / f'{region.region_id}.cutadapt.log'
        info_file_name = path_output_dir / f'{region.region_id}.cutadapt.info.gz'

        with open(path_demultiplex_exec, 'a') as f:
            cmd = f"""
                mv {path_output_fastq} {path_input_fastq}
                
                cutadapt {path_input_fastq} \\
                -o {path_output_fastq} \\
                -e {error_rate} \\
                -g ^file:{path_adapters} \\
                --rename '{rename_command}' \\
                --discard-untrimmed \\
                --json='{report_file_name}' \\
                --info-file='{info_file_name}' \\
                --cores={n_cores} 2>&1 | tee '{stdout_file_name}'
                
                """
            cmd = textwrap.dedent(cmd)
            f.write(cmd)

    with open(path_demultiplex_exec, 'a') as f:
        f.write(f'zgrep @ {path_output_fastq} | gzip -c > {path_final_reads}\n')
        f.write(f'delt-cli demultiplex compute-counts {path_final_reads} {path_counts}\n')
        f.write(f'rm {path_output_fastq} {path_input_fastq}')

    os.chmod(path_demultiplex_exec, os.stat(path_demultiplex_exec).st_mode | stat.S_IEXEC)

