from functools import reduce
from itertools import combinations, pairwise, product
from pathlib import Path

from Levenshtein import distance
from matplotlib import pyplot as plt
import numpy as np


def read_txt(
        path: Path,
) -> list:
    with open(path, 'r') as f:
        return f.readlines()


def plot_edit_distances(
        codons: dict,
        output_dir: Path,
) -> None:
    
    keys = set(filter(lambda x: 'C' not in x, set(codons.keys())))
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
        output_file = output_dir / f'edit-distance-{key}.png'
        fig.savefig(output_file)


def compute_stats(item, a, b):
    item['number_of_comparisons'] += 1
    item['number_of_matches'] += a == b
    item['proportion_of_matches'] = item['number_of_matches'] / item['number_of_comparisons']
    return item


def compute_overlap(
        codons: dict,
) -> dict:
    overlap = {}
    for a, b in pairwise(codons.keys()):
        last_base = [i.strip()[-1] for i in codons[a]]
        first_base = [i.strip()[0] for i in codons[b]]
        item = {'number_of_comparisons': 0, 'number_of_matches': 0, 'proportion_of_matches': 0}
        res = reduce(lambda x, y: compute_stats(x, a=y[0], b=y[1]), product(last_base, first_base), item)
        overlap[(a, b)] = res
    return overlap


def print_overlap(
        codons: dict,
) -> None:
    # NOTE: as expected around 1/4 of the codons in the const region end with the same base as the next region starts
    # The overall probability that a codon is not completely removed is given by
    # p = prob_of_insert * prop_of_matching_pre_postfix
    # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0976-y
    # p ~= 5*10^-6 * 0.25 ~= 10^-6
    overlap = compute_overlap(codons)
    for key, stats in overlap.items():
        print(f'\n{key}:')
        for stat, value in stats.items():
            print(f'{stat}:\t{value}')
    print('\n')


def analyze_codons(
        structure: dict,
        output_dir: Path,
) -> None:
    codons = {}
    for key, value in structure.items():
        codons[key] = read_txt(value['Path'])
    
    plot_edit_distances(codons, output_dir)
    print_overlap(codons)

