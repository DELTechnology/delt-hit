from pathlib import Path

import pandas as pd

from .. import demultiplex as d


def counts_are_identical(results: Path, legacy_results):
    legacy_results = legacy_results.assign(Code1=legacy_results.Code1 - 1, Code2=legacy_results.Code2 - 1)

    cols = ['Count', 'Code1', 'Code2']
    results = results[cols].sort_values(cols)
    legacy_results = legacy_results[cols].sort_values(cols)

    # TODO: remove after fixing the issue with the NF library
    bb1_legacy = pd.read_csv('/Users/adrianomartinelli/polybox/decl-data/raw-files/code1_NF2.txt')
    bb1_new = pd.read_csv(
        '/Users/adrianomartinelli/polybox/decl-data/db/experiments/OST1-e0-2024-07-12-20-29-10/codon_lists/B1.txt')
    pd.concat((bb1_legacy, bb1_new, bb1_legacy == bb1_new), axis=1)
    tmp = pd.concat((bb1_legacy, bb1_new, bb1_legacy == bb1_new), axis=1)

    results = results[results.Code1 <= 199]
    legacy_results = legacy_results[legacy_results.Code1 <= 199]

    return (results == legacy_results).all().all()


def get_selection_primer_ids_from_legacy_identifier(selection_name: str):
    return [int(i) - 1 for i in filter(len, selection_name.replace('selection_', '').split('_'))]


def compare_counts_with_legacy(config: dict, legacy_results_dir: Path):
    root = Path(config['Root']) / 'evaluations'
    hash = d.utils.hash_dict(config['Structure'])

    if not len(list(legacy_results_dir.glob('selection*.txt'))) == len(list(root.glob(f'selection-*/{hash}.txt'))):
        raise ValueError('Number of selections does not match between legacy and new results.')

    for legacy_result_path in legacy_results_dir.glob('selection*.txt'):
        selection_primer_ids = get_selection_primer_ids_from_legacy_identifier(legacy_result_path.stem)
        selection_id = d.postprocess.get_selection_ids([selection_primer_ids], config=config)[0]

        results_legacy = pd.read_csv(legacy_result_path, sep='\t')
        result_path = root / f'selection-{selection_id}' / f'{hash}.txt'
        assert result_path.exists()
        results = pd.read_csv(result_path, sep='\t')

        if not counts_are_identical(results, results_legacy):
            print(f'Evaluation results differ for original {legacy_result_path.name}')
