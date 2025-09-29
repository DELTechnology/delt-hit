import multiprocessing
import os
import stat
import textwrap
from pathlib import Path

import pandas as pd

from delt_core.utils import read_yaml
from delt_core.demultiplex.validation import SelectionFile, Region, validate

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


def get_codons(name: str, whitelists: dict) -> list[str]:
    return [item['codon'] for item in whitelists[name]]


def get_regions(structure: list[dict], whitelists: dict) -> list[Region]:
    def unique_codons(codons: list[str]) -> list[str]:
        seen = set()
        unique = []
        for codon in codons:
            if codon not in seen:
                unique.append(codon)
                seen.add(codon)
        return unique

    return [Region(
        name=item['name'],
        index=i,
        codons=unique_codons(get_codons(item['name'], whitelists)),
        max_error_rate=item['max_error_rate'],
        indels=item['indels']
    )
        for i, item in enumerate(structure)]


def write_fastq_files(regions: list[Region], save_path: Path) -> None:
    for i, region in enumerate(regions):
        fastq = [f'>{region.id}.{index}\n{codon}'
                 for index, codon in enumerate(region.codons)]
        fastq = '\n'.join(fastq)
        with open(save_path / f'{region.id}.fastq', 'w') as f:
            f.write(fastq)


def generate_input_files(
        config_path: Path,
        write_json_file: bool = True,
        write_info_file: bool = False,
        fast_dev_run: bool = False,
) -> None:
    config = read_yaml(config_path)

    save_dir = Path(config['experiment']['save_dir']).expanduser().resolve()
    path_input_fastq = Path(config['experiment']['fastq_path']).expanduser().resolve()

    num_cores = config['experiment']['num_cores']
    num_cores = multiprocessing.cpu_count() if pd.isna(num_cores) else num_cores

    structure = config['structure']
    whitelists = config['whitelists']

    experiment_name = config['experiment']['name']
    cutadapt_input_files_dir = save_dir / experiment_name / 'demultiplex' / 'cutadapt_input_files'
    cutadapt_output_files_dir = save_dir / experiment_name / 'demultiplex' / 'cutadapt_output_files'

    cutadapt_input_files_dir.mkdir(parents=True, exist_ok=True)
    path_demultiplex_exec = cutadapt_input_files_dir / 'demultiplex.sh'

    path_final_reads = cutadapt_output_files_dir / 'reads_with_adapters.gz'
    path_output_fastq = cutadapt_output_files_dir / 'out.fastq.gz'
    path_counts = save_dir / 'selections'

    regions = get_regions(structure=structure, whitelists=whitelists)
    write_fastq_files(regions, save_path=cutadapt_input_files_dir)

    with open(path_demultiplex_exec, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('# make sure you installed pigz with `brew install pigz` to enable parallel processing\n\n')
        f.write(f'mkdir "{cutadapt_output_files_dir}"\n')
        f.write(f'mkdir "{path_counts}"\n')
        # NOTE: we symlink the fastq file we want to demultiplex
        f.write(f'ln -sf "{path_input_fastq}" "{path_output_fastq}"\n')

        if fast_dev_run:
            n_reads_for_fast_dev_run = 10000
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

    for region in regions:
        error_rate = region.max_error_rate
        indels = f' --no-indels' if not int(region.indels) else ''
        path_adapters = cutadapt_input_files_dir / f'{region.id}.fastq'
        # NOTE: from now on we use the output of the previous step as input
        path_input_fastq = cutadapt_output_files_dir / 'input.fastq.gz'

        report_file_name = cutadapt_output_files_dir / f'{region.id}.cutadapt.json'
        stdout_file_name = cutadapt_output_files_dir / f'{region.id}.cutadapt.log'
        info_file_name = cutadapt_output_files_dir / f'{region.id}.cutadapt.info.gz'

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

            cmd += f'--cores={num_cores} 2>&1 | tee "{stdout_file_name}"\n'

            f.write(cmd)

    with open(path_demultiplex_exec, 'a') as f:
        f.write(f'\nzgrep @ "{path_output_fastq}" | gzip -c > "{path_final_reads}" || exit\n')
        f.write(f'delt-hit demultiplex process --config_path="{config_path}" || exit\n')
        f.write(f'rm "{path_output_fastq}" "{path_input_fastq}"\n')

    os.chmod(path_demultiplex_exec, os.stat(path_demultiplex_exec).st_mode | stat.S_IEXEC)

    return path_demultiplex_exec