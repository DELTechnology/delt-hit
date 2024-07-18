import multiprocessing
import os
from pathlib import Path
import stat
import textwrap

import pandas as pd

from .utils import read_yaml
from .validation import validate, Region, SelectionFile


def get_selections(
        config: dict,
        selection_id: int = None,
) -> pd.DataFrame:
    root = config['Root']
    config_selection = config['Selection']
    selection_file = root / config_selection['SelectionFile']
    fastq_file = str(Path(config_selection['FASTQFile']).name)
    library = str(Path(config_selection['Library']).name)
    selections = validate(pd.read_excel(selection_file), SelectionFile)
    if selection_id:
        selections = selections[selections['SelectionID'] == selection_id]
    return selections[(selections['FASTQFile'] == fastq_file) & (selections['Library'] == library)]


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
        write_json_file: bool,
        write_info_file: bool,
        fast_dev_run: bool = False,
) -> None:
    config = read_yaml(config_file)
    experiment_name = config['Experiment']['Name']
    cutadapt_input_files_dir = root_dir / 'experiments' / experiment_name / 'cutadapt_input_files'
    cutadapt_output_files_dir = root_dir / 'experiments' / experiment_name / 'cutadapt_output_files'

    cutadapt_input_files_dir.mkdir(parents=True, exist_ok=True)
    path_demultiplex_exec = cutadapt_input_files_dir / 'demultiplex.sh'

    path_final_reads = cutadapt_output_files_dir / 'reads_with_adapters.gz'
    path_output_fastq = cutadapt_output_files_dir / 'out.fastq.gz'
    path_counts = root_dir / 'evaluations'

    regions = get_regions(structure)
    write_fastq_files(regions, cutadapt_input_files_dir)

    with open(path_demultiplex_exec, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('# make sure you installed pigz with `brew install pigz` to enable parallel processing\n\n')
        f.write(f'mkdir "{cutadapt_output_files_dir}"\n')
        f.write(f'mkdir "{path_counts}"\n')
        # NOTE: we symlink the fastq file we want to demultiplex
        f.write(f'ln -sf "{path_input_fastq}" "{path_output_fastq}"\n')

        if fast_dev_run:
            n_reads_for_fast_dev_run = 1000
            n_lines = 4 * n_reads_for_fast_dev_run
            f.write(f'# fast-dev-run enabled\n')

            cmd = f"""
            tmp_file=$(mktemp)
            gzcat "{path_output_fastq}" | head -n {n_lines} | gzip -c > "$tmp_file"
            mv $tmp_file "{path_output_fastq}"
            """

            cmd = textwrap.dedent(cmd)
            f.write(cmd)

    rename_command = '{id} {comment}?{adapter_name}'
    n_cores = multiprocessing.cpu_count()

    for region in regions:
        error_rate = region.max_error_rate
        indels = f' --no-indels' if not int(region.indels) else ''
        path_adapters = cutadapt_input_files_dir / f'{region.region_id}.fastq'
        # NOTE: from now on we use the output of the previous step as input
        path_input_fastq = cutadapt_output_files_dir / 'input.fastq.gz'

        report_file_name = cutadapt_output_files_dir / f'{region.region_id}.cutadapt.json'
        stdout_file_name = cutadapt_output_files_dir / f'{region.region_id}.cutadapt.log'
        info_file_name = cutadapt_output_files_dir / f'{region.region_id}.cutadapt.info.gz'

        with open(path_demultiplex_exec, 'a') as f:
            cmd = f"""
                mv "{path_output_fastq}" "{path_input_fastq}"
                
                cutadapt "{path_input_fastq}" \\
                -o "{path_output_fastq}" \\
                -e {error_rate}{indels} \\
                -g "^file:{path_adapters}" \\
                --rename '{rename_command}' \\
                --discard-untrimmed \\
                """

            cmd = textwrap.dedent(cmd)

            if write_json_file:
                cmd += f'--json="{report_file_name}" \\\n'
            if write_info_file:
                cmd += f'--info-file="{info_file_name}" \\\n'

            cmd += f'--cores={n_cores} 2>&1 | tee "{stdout_file_name}"\n'

            f.write(cmd)

    with open(path_demultiplex_exec, 'a') as f:
        f.write(f'zgrep @ "{path_output_fastq}" | gzip -c > "{path_final_reads}"\n')
        f.write(f'delt-cli demultiplex compute-counts "{config_file}" "{path_final_reads}" "{path_counts}"\n')
        f.write(f'rm "{path_output_fastq}" "{path_input_fastq}"\n')

    os.chmod(path_demultiplex_exec, os.stat(path_demultiplex_exec).st_mode | stat.S_IEXEC)
