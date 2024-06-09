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

path_exp_dir = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B')
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
# fig.show()

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

# %% const region associated with 2.B.317
import gzip
k_n_most_mismatched = 4
index = more_than_expect_wE['index'].iloc[k_n_most_mismatched]
adapter_name = f'2.B.{index}'
adapter_seq = more_than_expect_wE['codon'].iloc[k_n_most_mismatched]

p_input = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B/eval/2.B.cutadapt.info.gz')
sequence_prefixes = []
prefix = slice(0, 25)
with gzip.open(p_input, 'rb') as f:
    _id = f.readline().decode('utf-8')
    while _id:
        is_trimmed = _id.split('\t')[1]
        if is_trimmed != '-1':
            _name, _, _, _, _, _seq, *_ = _id.split('\t')
            if adapter_name == _name.split('?')[-1]:
                sequence_prefixes.append(_seq)
        _id = f.readline().decode('utf-8')

labels, counts = np.unique(sequence_prefixes, return_counts=True)
idcs = np.argsort(counts)[::-1]
labels, counts = labels[idcs], counts[idcs]
labels = labels.astype(str)

top_k = 5
fig, ax = plt.subplots(figsize=(5, 1.5))
g = ax.barh(labels[:top_k], counts[:top_k], color='orange')
ax.set_yticklabels([])

x_pad = ax.get_xlim()[1] * 0.5 / len(adapter_seq)
for i, seq in enumerate(labels[:top_k]):
    t = ax.get_yticklabels()[i]
    x, y = t.get_position()
    for j, char in enumerate(seq):
        x += x_pad
        # color = 'green' if char == adapter_seq[j] else 'black'
        color = 'black'
        ax.text(x, y, char, verticalalignment='center', color=color)

t = ax.get_yticklabels()[i]
x, y = t.get_position()
y += 1
for j, char in enumerate(adapter_seq):
    x += x_pad
    ax.text(x, y, char, verticalalignment='center', color='black')

fig.tight_layout()
fig.show()
fig.savefig(path_eval_dir / f'3.C.{index}_trimmed_sequences_distribution.pdf')