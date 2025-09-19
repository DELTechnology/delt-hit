from pathlib import Path

import pandas as pd

import delt_core.utils
from .. import demultiplex as d


def counts_are_identical(results: Path, legacy_results):
    legacy_results = legacy_results.assign(Code1=legacy_results.Code1 - 1, Code2=legacy_results.Code2 - 1)
    cols = ['Count', 'Code1', 'Code2']
    results = results[cols].sort_values(cols)
    legacy_results = legacy_results[cols].sort_values(cols)
    return (results == legacy_results).all().all()


def get_selection_primer_ids_from_legacy_identifier(selection_name: str):
    return [int(i) - 1 for i in filter(len, selection_name.replace('selection_', '').split('_'))]


def compare_counts_with_legacy(config: dict, legacy_results_dir: Path):
    root = Path(config['Root']) / 'evaluations'
    _hash = delt_core.utils.hash_dict(config['Structure'])

    for legacy_result_path in legacy_results_dir.glob('selection*.txt'):
        selection_primer_ids = get_selection_primer_ids_from_legacy_identifier(legacy_result_path.stem)
        selection_id = d.postprocess.get_selection_ids([selection_primer_ids], config=config)[0]

        results_legacy = pd.read_csv(legacy_result_path, sep='\t')
        result_path = root / f'selection-{selection_id}' / f'{_hash}.txt'
        if not result_path.exists():
            print(f'⛔️ no results found ({selection_id} ↔️ {legacy_result_path.name})')
            continue

        results = pd.read_csv(result_path, sep='\t')

        if not counts_are_identical(results, results_legacy):
            print(f'⛔️ results differ ({selection_id} ↔️ {legacy_result_path.name})')
        else:
            print(f'✅ results identical ({selection_id} ↔️ {legacy_result_path.name})')

