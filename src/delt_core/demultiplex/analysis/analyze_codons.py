f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/code1_NF2.txt'
with open(f, 'r') as f:
    b1 = f.read().split('\n')
assert len(b1) == len(set(b1))
b1 = set(filter(len, b1))

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/code2_NF2.txt'
with open(f, 'r') as f:
    b2 = f.read().split('\n')
assert len(b2) == len(set(b2))
b2 = set(filter(len, b2))

b1.intersection(b2)

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/PrimerCodes1.txt'
with open(f, 'r') as f:
    p1 = f.read().split('\n')
assert len(p1) == len(set(p1))
p1 = set(filter(len, p1))

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/NF2_primers1b.txt'
with open(f, 'r') as f:
    p2 = f.read().split('\n')
assert len(p2) == len(set(p2))
p2 = set(filter(len, p2))

# not true, i.e. we have duplicates
len(b1.union(b2).union(p1).union(p2)) == len(b1) + len(b2) + len(p1) + len(p2)

# here we have an intesection
b1.intersection(p2)

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/const1.txt'
with open(f, 'r') as f:
    c1 = f.read().split('\n')
assert len(c1) == 1
# c1 = c1[0]

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/const2.txt'
with open(f, 'r') as f:
    c2 = f.read().split('\n')
assert len(c2) == 1
# c2 = c2[0]

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/const3.txt'
with open(f, 'r') as f:
    c3 = f.read().split('\n')
assert len(c3) == 1
# c3 = c3[0]

# %% dict of codons
codons = dict(b1=b1, b2=b2, p1=p1, p2=p2, c1=c1, c2=c2, c3=c3)

# %% find edit distance between codons
from Levenshtein import distance
from itertools import product, combinations
from matplotlib import pyplot as plt
import numpy as np

keys = set(filter(lambda x: 'c' not in x, set(codons.keys())))
for key in keys:
    other = keys - {key}
    codons_key = codons[key]
    codons_other = []
    for i in other:
        codons_other.extend(codons[i])

    d = [distance(*i) for i in combinations(codons_key, 2)]
    d_other = [distance(*i) for i in product(codons_key, codons_other)]

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    labels, counts = np.unique(d, return_counts=True)
    bars = axs[0].bar(labels, counts, align='center')
    axs[0].set_title('Intra codon set edit distance')

    # Add counts below each bar
    for bar, count in zip(bars, counts):
        axs[0].text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
                    ha='center', va='bottom', color='black')

    labels, counts = np.unique(d_other, return_counts=True)
    bars = axs[1].bar(labels, counts, align='center')
    axs[1].set_title('Inter codon set edit distance')
    # Add counts below each bar
    for bar, count in zip(bars, counts):
        axs[1].text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
                    ha='center', va='bottom', color='black')

    fig.suptitle(f'{key} codon set')
    fig.subplots_adjust(top=0.85)
    fig.tight_layout()
    fig.show()
    fig.savefig(f'edit-distance-{key}.png')

# %% find number of codons that end with the same base as the const region starts
from itertools import pairwise
from functools import reduce
regions_ord = ['p1', 'c1', 'b1', 'c2', 'b2', 'c3', 'p2']

# 0, 1, 2
reduce(lambda x,y: x['a'] + [y], range(3), {'a': []})

reduce(lambda x, y: x['a'] + [y], range(3), {'a': []})

def compute_stats(d, a, b):
    d['number_of_comparisons'] += 1
    d['number_of_matches'] += a == b
    return d

overlap = {}
for a,b in pairwise(regions_ord):
    last_base = [i[-1] for i in codons[a]]
    first_base = [i[0] for i in codons[b]]

    item = {'number_of_comparisons': 0, 'number_of_matches': 0}
    res = reduce(lambda x,y: compute_stats(x, a=y[0], b=y[1]), product(last_base, first_base), item)
    overlap[(a, b)] = res

# NOTE: as expected around 1/4 of the codons in the const region end with the same base as the next region starts
for k, v in overlap.items():
    overlap[k]['propportion_of_matches'] = v['number_of_matches'] / v['number_of_comparisons']

# %%
# The overall probability that a codon is not completely removed is given by
# p = prob_of_insert * prop_of_matching_pre_postfix
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0976-y
# p ~= 5*10^-6 * 0.25 ~= 10^-6