import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

p_adapter = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-3.C/cutadapt_input_files/3.C.fastq')
adapter_seq = open(p_adapter).readlines()[-1].strip()

p_input = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-3.C/eval/input.fastq.gz')
sequence_prefixes = []
prefix = slice(0, 25)
with gzip.open(p_input, 'rb') as f:
    _id = f.readline()
    while _id:
        _seq = f.readline()
        _ = f.readline()
        _phred = f.readline()
        sequence_prefixes.append(_seq[prefix])
        _id = f.readline()

labels, counts = np.unique(sequence_prefixes, return_counts=True)
idcs = np.argsort(counts)[::-1]
labels, counts = labels[idcs], counts[idcs]
labels = labels.astype(str)

def format_sequences(L, labels):
    formatted_labels = [L]
    for seq in labels:
        formatted_seq = ''.join([f'\033[32m{char}\033[0m' if char == L[i] else char for i, char in enumerate(seq)])
        formatted_labels.append(formatted_seq)
    return formatted_labels


# %%
top_k = 5
fig, ax = plt.subplots(figsize=(5, 1.5))
g = ax.barh(labels[:top_k], counts[:top_k], color='orange')
# ax.set_yticklabels(ax.get_yticklabels(), horizontalalignment='left')
ax.set_yticklabels([])
# ax.yaxis.set_tick_params(pad=-10)

x_pad = ax.get_xlim()[1] * 0.8 / len(adapter_seq)
for i, seq in enumerate(labels[:top_k]):
    t = ax.get_yticklabels()[i]
    x, y = t.get_position()
    for j, char in enumerate(seq):
        x += x_pad
        color = 'green' if char == adapter_seq[j] else 'black'
        ax.text(x, y, char, verticalalignment='center', color=color)

t = ax.get_yticklabels()[i]
x, y = t.get_position()
y += 1
for j, char in enumerate(adapter_seq):
    x += x_pad
    ax.text(x, y, char, verticalalignment='center', color='black')

# %%

formatted_labels = format_sequences(adapter_seq, labels[:top_k].astype(str))
text_to_display = '\n'.join(formatted_labels)
print(text_to_display)

p_eval_dir = p_input.parent
with open(p_eval_dir / '3.C.sequence_alignment.txt', 'w') as f:
    f.write(text_to_display + '\n')

fig.tight_layout()
fig.show()
fig.savefig(p_eval_dir / '3.C.prefix_distribution.pdf')
