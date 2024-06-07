import multiprocessing
import os
import stat
import textwrap
import typing as tp

from pydantic import BaseModel, computed_field


class Region(BaseModel):
    start: int
    end: int
    type: str
    codons: list[str]
    position_in_construct: int = None

    @computed_field
    @property
    def region_id(self) -> str:
        return f'{self.position_in_construct}.{self.type}'


def get_regions(
        struct_file = str,
):
    with open(struct_file, 'r') as f:
        struct = f.read().split('\n')
    
    regions = []
    for i in struct:
        if len(i) == 0:
            continue
        idx_start, idx_end, _type, relative_path_to_codons = i.split('\t')

        with open(struct_file.parent / relative_path_to_codons, 'r') as f:
            codons = f.read().split('\n')
            codons = filter(len, codons)

        region = Region(start=int(idx_start) - 1, end=int(idx_end) - 1, type=_type, codons=codons)
        regions.append(region)

    return sorted(regions, key=lambda x: x.start)


def write_fastq_files(
        regions = tp.List,
        path = str,
):
    fastq = ''
    for i, region in enumerate(regions):
        region.position_in_construct = i
        fastq = [f'>{region.region_id}.{index}\n{codon}'
                for index, codon in enumerate(region.codons)]
        fastq = '\n'.join(fastq)

        with open(path / f'{region.region_id}.fastq', 'w') as f:
            f.write(fastq)


def generate_input_files(
        path_input_fastq: str,
        path_struct_file: str,
        path_demultiplex_exec: str,
        path_final_reads: str,
):
    path_input_dir = path_demultiplex_exec.parent
    path_output_dir = path_final_reads.parent
    path_input_dir.mkdir(parents=True, exist_ok=True)
    path_output_dir.mkdir(parents=True, exist_ok=True)
    path_output_fastq = path_output_dir / 'out.fastq.gz'

    regions = get_regions(path_struct_file)
    write_fastq_files(regions, path_input_dir)

    with open(path_demultiplex_exec, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('# make sure you installed pigz with `brew install pigz` to enable parallel processing\n\n')
        f.write(f'mkdir "{path_output_dir}"\n')
        # we symlink the fastq file we want to demultiplex
        f.write(f'ln -sf "{path_input_fastq}" "{path_output_fastq}"\n')

    error_rate = 0
    rename_command = '{id} {comment}?{adapter_name}'
    n_cores = multiprocessing.cpu_count()

    for region in regions:
        # path_output_info = path_output_dir / 'info.tsv'
        path_adapters = path_input_dir / f'{region.region_id}.fastq'
        # note: from now on we use the output of the previous step as input
        path_input_fastq = path_output_dir / 'input.fastq.gz'

        report_file_name = path_output_dir / f'{region.region_id}.cutadapt.json'
        stdout_file_name = path_output_dir / f'{region.region_id}.cutadapt.log'
        info_file_name = path_output_dir / f'{region.region_id}.cutadapt.info.gz'

        with open(path_demultiplex_exec, 'a') as f:
            cmd = f"""
                mv {path_output_fastq} {path_input_fastq}
                
                cutadapt "{path_input_fastq}" \\
                -o "{path_output_fastq}" \\
                -e {error_rate} \\
                -g "^file:{path_adapters}" \\
                --rename '{rename_command}' \\
                --discard-untrimmed \\
                --json='{report_file_name}' \\
                --info-file='{info_file_name}' \\
                --cores={n_cores} 2>&1 | tee '{stdout_file_name}'
                
                """
            cmd = textwrap.dedent(cmd)
            f.write(cmd)

    with open(path_demultiplex_exec, 'a') as f:
        f.write(f'zgrep "@" "{path_output_fastq}" | gzip -c > "{path_final_reads}"\n')
        f.write(f'rm {path_output_fastq} {path_input_fastq}')

    os.chmod(path_demultiplex_exec, os.stat(path_demultiplex_exec).st_mode | stat.S_IEXEC)

