from pathlib import Path
import json

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from Levenshtein import distance
from itertools import product, combinations
from matplotlib import pyplot as plt
import numpy as np


path_exp_dir = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e1')
path_eval_dir = path_exp_dir / 'eval'

f = path_exp_dir / 'cutadapt_input_files' / '2.B.fastq'
b1 = []
with open(f, 'r') as f:
    while True:
        _id = f.readline()

        if not _id:
            break

        b1.append(f.readline().strip())

assert len(b1) == len(set(b1))


d = [{'i': i[0], 'j': i[1], 'd': distance(b1[i[0]], b1[i[1]])} for i in combinations(range(len(b1)), 2)]
d = pd.DataFrame(d)

fig, ax = plt.subplots()
labels, counts = np.unique(d['d'], return_counts=True)
bars = ax.bar(labels, counts, align='center')
ax.set_title('2.B codon set edit distance')

# Add counts below each bar
for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
                ha='center', va='bottom', color='black')
fig.savefig(path_eval_dir / '2.B_edit_distances.pdf')
fig.show()

report_path = path_eval_dir / '2.B.cutadapt.json'

report = json.load(report_path.open('r'))
stats = report['adapters_read1']
reads_in = report['read_counts']['input']
reads_out = report['read_counts']['output']

cont = []
for stat in stats:
    if stat['five_prime_end']['trimmed_lengths']:
        item = {'reads_in': reads_in, 'reads_out': reads_out}

        item['region_id'] = '_'.join(stat['name'].split('.')[:-1])
        item['index'] = int(stat['name'].split('.')[-1])

        ec = pd.concat([pd.DataFrame(i) for i in stat['five_prime_end']['trimmed_lengths']])
        ec = ec.reset_index().rename(columns={'index': 'number_of_errors', 'counts': 'error_counts'})
        ec = ec.assign(index=item['index'])
        item['error_counts'] = ec

        cont.append(item)

df = pd.DataFrame(cont)
ec = pd.concat(df.error_counts.tolist())

most_counts_woE = ec.loc[ec.number_of_errors == 0, :].sort_values('error_counts', ascending=False)
most_counts_woE = most_counts_woE.assign(codon = [b1[i] for i in most_counts_woE['index']])

more_than_expect = ec[ec.expect < ec.error_counts]
print(f'Number of more than expected counts by number of errors {more_than_expect.number_of_errors.value_counts()}')

more_than_expect_wE = more_than_expect[more_than_expect.number_of_errors > 0]
more_than_expect_wE = more_than_expect_wE.sort_values('error_counts', ascending=False)
more_than_expect_wE = more_than_expect_wE.assign(codon = [b1[i] for i in more_than_expect_wE['index']])

def get_edit_distances_for_codon(index):
    return d[(d.i == index) | (d.j ==index)]

e = get_edit_distances_for_codon(17)


for grp_name, grp_dat in df.groupby('region_id'):
    # grp_dat = pd.concat([grp_dat, grp_dat.assign(number_of_errors=1)])

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

    # axs[1].plot([0, grp_dat.index.max()], [expected, expected], 'k--')


    # g = sns.barplot(data=grp_dat,
    #                 x='index', y='error_counts',
    #                 hue='number_of_errors',
    #                 estimator='sum',
    #                 ax=axs[0])
    # g = sns.barplot(data=grp_dat,
    #                     x='index', y='error_counts',
    #                     hue='number_of_errors',
    #                 estimator='sum',
    #                     ax=axs[1])

    for ax in axs:
        ax.set_xticklabels([])

    fig.tight_layout()
    fig.savefig(path_eval_dir / f'hits_{grp_name}.pdf')
    fig.show()


report = json.load((path_eval_dir / '2.B.cutadapt.json').open('r'))
stats = report['adapters_read1']
stats = sorted(stats, key=lambda x: x['total_matches'])
stats[-1]
