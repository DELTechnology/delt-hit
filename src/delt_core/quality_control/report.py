import json
from pathlib import Path
import textwrap

import numpy as np


def print_report(
        experiment_path: Path,
) -> None:
    report_files = (experiment_path / 'cutadapt_output_files').glob('*.cutadapt.json')
    reports = {i.name.replace('.cutadapt.json', ''): json.load(i.open('r')) for i in report_files}

    pipeline_report = []
    for region_id, report in reports.items():
        c_input = report['read_counts']['input']
        c_output = report['read_counts']['output']
        if c_input == 0:
            c_output_proportion = np.nan
        else:
            c_output_proportion = c_output / c_input

        region_report = {
            'region_id': region_id,
            'total_reads_processed': c_input,
            'reads_with_adapters': c_output,
            'proportion_of_reads_with_adapters': c_output_proportion,
            'reads_discarded': c_input - c_output,
            'proportion_of_reads_discarded': 1 - c_output_proportion
        }

        pipeline_report.append(region_report)

    pipeline_report = sorted(pipeline_report, key=lambda x: x['region_id'])

    # ANSI escape codes for colors
    RED = '\033[91m'
    GREEN = '\033[92m'
    RESET = '\033[0m'

    for rr in pipeline_report:
        s = f"""Region {rr['region_id']}: {rr['total_reads_processed']:,} {RED}→{RESET} {rr['reads_discarded']:,} ({RED if rr['proportion_of_reads_discarded'] > 0.1 else ''}{rr['proportion_of_reads_discarded']:.2%}{RESET if rr['proportion_of_reads_discarded'] > 0.1 else ''}) reads discarded
            \t{GREEN}↓{RESET} ({rr['proportion_of_reads_with_adapters']:.2%} reads with adapters)"""
        s = textwrap.dedent(s)
        print(s)

    r0 = pipeline_report[0]
    rr = pipeline_report[-1]

    overall_proportion_reads_with_adapter = rr['reads_with_adapters'] / r0['total_reads_processed']
    overall_reads_discarded = r0['total_reads_processed'] - rr['reads_with_adapters']
    overall_proportion_reads_discarded = 1 - overall_proportion_reads_with_adapter

    s = \
        f"""{'Overall':<10}: {GREEN}{rr['reads_with_adapters']:,}{RESET} reads with adapters ({overall_proportion_reads_with_adapter:.2%})
    {'':<11}{RED}{overall_reads_discarded:,}{RESET} ({overall_proportion_reads_discarded:.2%}) reads discarded"""
    s = textwrap.dedent(s)
    with open(experiment_path / 'report.txt', 'w') as f:
        f.write(s)
    print(s)

