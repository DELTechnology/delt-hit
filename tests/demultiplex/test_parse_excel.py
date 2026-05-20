from pathlib import Path

import pandas as pd
import pytest

from delt_hit.demultiplex.parser import config_from_excel, structure_from_excel, whitelists_from_excel


def _write_excel_config(path, *, include_b1_flag_columns: bool) -> None:
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

        pd.DataFrame(
            {
                'codon': ['GCCTCG', 'AATTCC', 'TTGGCC'],
                'reverse': [True, None, 'TRUE'],
                'complement': [True, None, 'TRUE'],
            }
        ).to_excel(writer, sheet_name='B0', index=False)

        b1_df = pd.DataFrame({'codon': ['AGTCGA']})
        if include_b1_flag_columns:
            b1_df['reverse'] = [None]
            b1_df['complement'] = [None]
        b1_df.to_excel(writer, sheet_name='B1', index=False)


def test_whitelists_from_excel_normalizes_reverse_and_complement_flags_for_building_blocks(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_excel_config(excel_path, include_b1_flag_columns=True)

    whitelists = whitelists_from_excel(excel_path)

    assert whitelists['B0'][:3] == [
        {'index': 0, 'codon': 'GCCTCG', 'reverse': True, 'complement': True},
        {'index': 1, 'codon': 'AATTCC', 'reverse': False, 'complement': False},
        {'index': 2, 'codon': 'TTGGCC', 'reverse': True, 'complement': True},
    ]
    assert whitelists['B1'] == [{'index': 0, 'codon': 'AGTCGA', 'reverse': False, 'complement': False}]
    assert whitelists['S0'] == [
        {'name': 'sel_a', 'codon': 'AAAA'},
        {'name': 'sel_b', 'codon': 'CCCC'},
    ]
    assert whitelists['C0'] == [{'codon': 'ACAC'}]


def test_whitelists_from_excel_defaults_missing_reverse_and_complement_columns_to_false(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_excel_config(excel_path, include_b1_flag_columns=False)

    whitelists = whitelists_from_excel(excel_path)

    assert whitelists['B1'] == [{'index': 0, 'codon': 'AGTCGA', 'reverse': False, 'complement': False}]


def test_example_dual_display_workbook_includes_reverse_and_complement_flags():
    workbook_path = (
        Path(__file__).resolve().parents[2]
        / 'supporting_material'
        / 'experiments'
        / 'example-dual-display'
        / 'example-dual-display-2.xlsx'
    )

    whitelists = whitelists_from_excel(workbook_path)

    for sheet_name in ('B0', 'B1'):
        assert all(
            row['reverse'] is False and row['complement'] is False
            for row in whitelists[sheet_name]
        )

    for sheet_name in ('B2', 'B3'):
        assert all(
            row['reverse'] is True and row['complement'] is True
            for row in whitelists[sheet_name]
        )


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
