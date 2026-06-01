from pathlib import Path

import pandas as pd
import pytest

from delt_hit.demultiplex.parser import config_from_excel, experiment_from_excel, structure_from_excel, whitelists_from_excel


def _write_excel_config(path, *, include_structure_flag_columns: bool) -> None:
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame(
            {
                'name': ['sel_a', 'sel_b'],
                'S0': ['AAAA', 'CCCC'],
                'S1': ['GGGG', 'TTTT'],
            }
        ).to_excel(writer, sheet_name='selection', index=False)

        pd.DataFrame(
            {
                'name': ['C0'],
                'codon': ['ACAC'],
            }
        ).to_excel(writer, sheet_name='constant', index=False)

        structure_df = pd.DataFrame(
            {
                'name': ['B0', 'B1', 'S0', 'S1', 'C0'],
                'type': ['building_block', 'building_block', 'selection', 'selection', 'constant'],
            }
        )
        if include_structure_flag_columns:
            structure_df['reverse'] = [True, None, 'TRUE', False, None]
            structure_df['complement'] = [True, None, 'TRUE', False, None]
        structure_df.to_excel(writer, sheet_name='structure', index=False)

        pd.DataFrame(
            {
                'codon': ['GCCTCG', 'AATTCC', 'TTGGCC'],
            }
        ).to_excel(writer, sheet_name='B0', index=False)

        pd.DataFrame({'codon': ['AGTCGA']}).to_excel(writer, sheet_name='B1', index=False)


def test_structure_from_excel_normalizes_reverse_and_complement_flags(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_excel_config(excel_path, include_structure_flag_columns=True)

    structure = structure_from_excel(excel_path)

    assert structure == [
        {'name': 'B0', 'type': 'building_block', 'reverse': True, 'complement': True},
        {'name': 'B1', 'type': 'building_block', 'reverse': False, 'complement': False},
        {'name': 'S0', 'type': 'selection', 'reverse': True, 'complement': True},
        {'name': 'S1', 'type': 'selection', 'reverse': False, 'complement': False},
        {'name': 'C0', 'type': 'constant', 'reverse': False, 'complement': False},
    ]


def test_structure_from_excel_defaults_missing_reverse_and_complement_columns_to_false(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_excel_config(excel_path, include_structure_flag_columns=False)

    structure = structure_from_excel(excel_path)

    assert structure == [
        {'name': 'B0', 'type': 'building_block', 'reverse': False, 'complement': False},
        {'name': 'B1', 'type': 'building_block', 'reverse': False, 'complement': False},
        {'name': 'S0', 'type': 'selection', 'reverse': False, 'complement': False},
        {'name': 'S1', 'type': 'selection', 'reverse': False, 'complement': False},
        {'name': 'C0', 'type': 'constant', 'reverse': False, 'complement': False},
    ]


def test_experiment_from_excel_strips_surrounding_whitespace_from_string_entries(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    with pd.ExcelWriter(excel_path) as writer:
        pd.DataFrame(
            {
                'variable': ['name', 'save_dir', 'fastq_path'],
                'value': ['  demo  ', f'  {tmp_path}  ', '  reads/sample.fastq.gz  '],
            }
        ).to_excel(writer, sheet_name='experiment', index=False)

    experiment = experiment_from_excel(excel_path)

    assert experiment == {
        'name': 'demo',
        'save_dir': str(tmp_path),
        'fastq_path': 'reads/sample.fastq.gz',
    }


def test_whitelists_from_excel_does_not_preserve_reverse_or_complement_flags(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_excel_config(excel_path, include_structure_flag_columns=True)

    whitelists = whitelists_from_excel(excel_path)

    assert whitelists['B0'][:3] == [
        {'index': 0, 'codon': 'GCCTCG'},
        {'index': 1, 'codon': 'AATTCC'},
        {'index': 2, 'codon': 'TTGGCC'},
    ]
    assert whitelists['B1'] == [{'index': 0, 'codon': 'AGTCGA'}]


def test_example_dual_display_extended_workbook_includes_structure_level_flags():
    workbook_path = (
        Path(__file__).resolve().parents[2]
        / 'supporting_material'
        / 'experiments'
        / 'example-dual-display'
        / 'example-dual-display-extended.xlsx'
    )

    structure = structure_from_excel(workbook_path)
    structure_by_name = {entry['name']: entry for entry in structure}

    for region_name in ('B0', 'B1'):
        assert structure_by_name[region_name]['reverse'] is False
        assert structure_by_name[region_name]['complement'] is False

    for region_name in ('B2', 'B3'):
        assert structure_by_name[region_name]['reverse'] is True
        assert structure_by_name[region_name]['complement'] is True


def _write_structure_excel_config(path, *, strands: list[str | None], include_selection_strand: bool = False) -> None:
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame(
            {
                'variable': ['name', 'save_dir'],
                'value': ['dual', str(path.parent)],
            }
        ).to_excel(writer, sheet_name='experiment', index=False)

        pd.DataFrame(
            {
                'name': ['sel_a'],
                'S0': ['AAAA'],
                'S1': ['CCCC'],
            }
        ).to_excel(writer, sheet_name='selection', index=False)

        pd.DataFrame(
            {
                'name': ['C0'],
                'codon': ['ACAC'],
            }
        ).to_excel(writer, sheet_name='constant', index=False)

        selection_strand = 'A' if include_selection_strand else None
        pd.DataFrame(
            {
                'name': ['B0', 'B1', 'B2', 'B3', 'S0', 'S1', 'C0'],
                'type': ['building_block', 'building_block', 'building_block', 'building_block', 'selection', 'selection', 'constant'],
                'strand': [*strands, selection_strand, None, None],
            }
        ).to_excel(writer, sheet_name='structure', index=False)

        for idx, sheet in enumerate(['B0', 'B1', 'B2', 'B3']):
            pd.DataFrame(
                {
                    'codon': [f'CODE{idx}'],
                    'smiles': ['C'],
                    'reaction': ['rxn_join'],
                    'product': ['product_1'],
                    'educt': ['scaffold'],
                }
            ).to_excel(writer, sheet_name=sheet, index=False)

        pd.DataFrame(
            {
                'name': ['scaffold'],
                'smiles': ['CC'],
            }
        ).to_excel(writer, sheet_name='compounds', index=False)

        pd.DataFrame(
            {
                'name': ['rxn_join'],
                'smirks': ['[C:1].[C:2]>>[C:1][C:2]'],
            }
        ).to_excel(writer, sheet_name='reactions', index=False)

        pd.DataFrame(
            {
                'educt_1': ['scaffold'],
                'educt_2': [None],
                'reaction': ['rxn_join'],
                'product': ['product_1'],
            }
        ).to_excel(writer, sheet_name='reaction_graph', index=False)


def test_structure_from_excel_accepts_strand_for_building_blocks(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_structure_excel_config(excel_path, strands=['A', 'A', 'B', 'B'])

    structure = structure_from_excel(excel_path)

    building_blocks = [entry for entry in structure if entry['type'] == 'building_block']
    assert [entry['strand'] for entry in building_blocks] == ['A', 'A', 'B', 'B']


def test_structure_from_excel_rejects_strand_for_non_building_block_rows(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_structure_excel_config(excel_path, strands=['A', 'A', 'B', 'B'], include_selection_strand=True)

    with pytest.raises(AssertionError, match="only valid for rows with `type == 'building_block'`"):
        structure_from_excel(excel_path)


def test_config_from_excel_requires_complete_dual_display_strands(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_structure_excel_config(excel_path, strands=['A', None, 'B', 'B'])

    with pytest.raises(AssertionError, match='require `strand` for every building_block row'):
        config_from_excel(excel_path)


def test_config_from_excel_rejects_dual_display_educts_from_other_strand(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    with pd.ExcelWriter(excel_path) as writer:
        pd.DataFrame(
            {'variable': ['name', 'save_dir'], 'value': ['dual', str(tmp_path)]}
        ).to_excel(writer, sheet_name='experiment', index=False)
        pd.DataFrame(
            {'name': ['sel_a'], 'S0': ['AAAA'], 'S1': ['CCCC']}
        ).to_excel(writer, sheet_name='selection', index=False)
        pd.DataFrame(
            {'name': ['C0'], 'codon': ['ACAC']}
        ).to_excel(writer, sheet_name='constant', index=False)
        pd.DataFrame(
            {
                'name': ['B0', 'B1', 'B2', 'B3'],
                'type': ['building_block', 'building_block', 'building_block', 'building_block'],
                'strand': ['A', 'A', 'B', 'B'],
            }
        ).to_excel(writer, sheet_name='structure', index=False)
        pd.DataFrame(
            {'name': ['scaffold'], 'smiles': ['CC']}
        ).to_excel(writer, sheet_name='compounds', index=False)
        pd.DataFrame(
            {'name': ['rxn_1', 'rxn_2'], 'smirks': ['[C:1].[C:2]>>[C:1][C:2]', '[C:1].[C:2]>>[C:1][C:2]']}
        ).to_excel(writer, sheet_name='reactions', index=False)
        pd.DataFrame(
            {
                'educt_1': ['scaffold', 'product_1', 'scaffold', 'product_1'],
                'educt_2': ['B0', 'B1', 'B2', 'B3'],
                'reaction': ['rxn_1', 'rxn_2', 'rxn_1', 'rxn_2'],
                'product': ['product_1', 'product_2', 'product_3', 'product_2'],
            }
        ).to_excel(writer, sheet_name='reaction_graph', index=False)
        pd.DataFrame(
            {'codon': ['AAAA'], 'smiles': ['C'], 'educt': ['scaffold'], 'reaction': ['rxn_1'], 'product': ['product_1']}
        ).to_excel(writer, sheet_name='B0', index=False)
        pd.DataFrame(
            {'codon': ['BBBB'], 'smiles': ['C'], 'educt': ['product_1'], 'reaction': ['rxn_2'], 'product': ['product_2']}
        ).to_excel(writer, sheet_name='B1', index=False)
        pd.DataFrame(
            {'codon': ['CCCC'], 'smiles': ['C'], 'educt': ['scaffold'], 'reaction': ['rxn_1'], 'product': ['product_3']}
        ).to_excel(writer, sheet_name='B2', index=False)
        pd.DataFrame(
            {'codon': ['DDDD'], 'smiles': ['C'], 'educt': ['product_1'], 'reaction': ['rxn_2'], 'product': ['product_2']}
        ).to_excel(writer, sheet_name='B3', index=False)

    with pytest.raises(AssertionError, match='outside the B-strand branch'):
        config_from_excel(excel_path)
