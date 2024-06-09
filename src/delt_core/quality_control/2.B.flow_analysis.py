import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

p_adapter = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B/cutadapt_input_files/3.C.fastq')
adapter_seq = open(p_adapter).readlines()[-1].strip()

p_report = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B/eval/2.B.cutadapt.json')
report = json.load(open(p_report, 'r'))

# %% most common adapter

adapters = report['adapters_read1']
adapters = sorted(adapters, key=lambda x: x['total_matches'])
msg = f"Adapter {adapters[-1]['name']} accumulates {adapters[-1]['total_matches'] / sum(i['total_matches'] for i in adapters):.2%} of the matches\n\n"
print(msg)

# %%
p_outputs = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B/eval/').glob(
    'out.*.fastq.gz')
trimmed = {}
adapter_id = adapters[-1]['name']
for p_out in p_outputs:
    c = 0
    with gzip.open(p_out, 'rb') as f:
        _id = f.readline()
        while _id:
            if adapter_id in _id.decode('utf-8'):
                c += 1

            _seq = f.readline()
            _ = f.readline()
            _phred = f.readline()

            _id = f.readline()

    trimmed[p_out.stem.replace('out.', '').replace('.fastq', '')] = c

trimmed

# %%
p_outputs = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B/eval/').glob(
    '*untrimmed.output.gz')
untrimmed = {}
adapter_id = adapters[-1]['name']
for p_out in p_outputs:
    c = 0
    with gzip.open(p_out, 'rb') as f:
        _id = f.readline()
        while _id:
            if adapter_id in _id.decode('utf-8'):
                c += 1

            _seq = f.readline()
            _ = f.readline()
            _phred = f.readline()

            _id = f.readline()

    untrimmed[p_out.stem.replace('.untrimmed.output', '')] = c

untrimmed

# %% const region associated with 2.B.317
p_input = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B/eval/3.C.untrimmed.output.gz')
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


top_k = 5
fig, ax = plt.subplots(figsize=(5, 1.5))
g = ax.barh(labels[:top_k], counts[:top_k], color='orange')
ax.set_yticklabels([])

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
with open(p_eval_dir / '3.C_untrimmed_sequence_alignment.txt', 'w') as f:
    f.write(text_to_display + '\n')

fig.tight_layout()
fig.show()
fig.savefig(p_eval_dir / '3.C_untrimmed_prefix_distribution.pdf')


# %% print flow
s = 'Flow of reads with 2.B.317 adapters in subsequent trimming rounds\n'
_region = "2.B"
_in = trimmed[_region]
for i in untrimmed.keys():
    _t = trimmed[i]
    _ut = untrimmed[i]
    s += f'Region {i}: {_in} → {_ut} ({_ut/_in:.2%}) discarded\n'
    s += f'              ↓ {_t} ({_t/_in:.2%}) trimmed\n'

    _region = i
    _in = trimmed[i]

print(s)
with open(p_eval_dir / '2.B_flow_analysis.txt', 'w') as f:
    f.write(msg)
    f.write(s)


# %%

# untrimmed_node_ids = {i: int(i.split('.')[0] + '0') for i in untrimmed.keys()}
# trimmed_node_ids = {i: int(i.split('.')[0] + '1') for i in trimmed.keys()}
#
# trimmed_v2 = {int(i.split('.')[0] + '1'): v for i, v in trimmed.items()}
# untrimmed_v2 = {int(i.split('.')[0] + '0'): v for i, v in untrimmed.items()}
# values = {**trimmed_v2, **untrimmed_v2}
#
# source = sorted(trimmed_node_ids.values())
# target = sorted([s + 9 for s in source] + [s + 10 for s in source])[:-2]
# source = sorted(source * 2)[:-2]
# value = [values[t] for t in target]
# target = [i if i % 2 else 0 for i in target]
#
# import plotly.graph_objects as go
#
# link = dict(source=source, target=target, value=value)
# data = go.Sankey(link=link)
# 
# fig = go.Figure(data)
#
# fig.write_html(p_eval_dir / "sankey-diagram-plotly1.html")