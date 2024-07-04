import json
from pathlib import Path

from matplotlib import pyplot as plt
import pandas as pd


def get_stats(
        experiment_path: Path,
) -> list:
    items = []
    report = json.load(experiment_path.open('r'))
    stats = report['adapters_read1']
    reads_in = report['read_counts']['input']
    reads_out = report['read_counts']['output']
    for stat in stats:
        if stat['five_prime_end']['trimmed_lengths']:
            item = {'reads_in': reads_in, 'reads_out': reads_out}
            item['region_id'] = '_'.join(stat['name'].split('.')[:-1])
            item['index'] = int(stat['name'].split('.')[-1])
            ec = pd.concat([pd.DataFrame(i) for i in stat['five_prime_end']['trimmed_lengths']])
            ec = ec.reset_index().rename(columns={'index': 'number_of_errors', 'counts': 'error_counts'})
            ec = ec.assign(index=item['index'])
            item['error_counts'] = ec
            items += [item]
    return items


def plot_hits(
        experiment_dir: Path,
        output_dir: Path,
) -> None:
    stats = []
    report_paths = (experiment_dir / 'cutadapt_output_files').glob('*.cutadapt.json')
    for report_path in report_paths:
        stats += [*get_stats(report_path)]
    df = pd.DataFrame(stats)

    for grp_name, grp_dat in df.groupby('region_id'):
        pdat = pd.concat(grp_dat['error_counts'].tolist())
        pdat = pdat[['index', 'number_of_errors', 'error_counts']] \
            .groupby(['index', 'number_of_errors']) \
            .agg('sum') \
            .reset_index() \
            .pivot(index='index', columns='number_of_errors', values='error_counts')

        _in, _out, = grp_dat.reads_in.iloc[0], grp_dat.reads_out.iloc[0]
        # expected = grp_dat.expected.iloc[0]
        error_sums = {col: pdat[col].sum() for col in pdat}
        error_stats = ' '.join(f'e{key}: {val:.3E} ({val/_out:.2%})' for key, val in error_sums.items())

        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        fig.suptitle(f'{grp_name}, {_out:.3E} / {_in:.3E} ({_out/_in:.2%})\n{error_stats}')

        pdat.plot(kind='bar', stacked=True, ax=axs[0])
        pdat.plot(kind='bar', stacked=True, ax=axs[1])

        axs[0].set_ylabel('counts')
        _max = max(axs[1].get_ylim())
        _max += _max / 2
        ylim = (1, _max)
        axs[1].set_yscale('log')
        axs[1].set_ylim(ylim)

        for ax in axs:
            ax.set_xticklabels([])
        fig.tight_layout()
        fig.savefig(output_dir / f'hits_{grp_name}.pdf')

