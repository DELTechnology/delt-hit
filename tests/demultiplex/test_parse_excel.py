import pandas as pd

from delt_hit.demultiplex.parser import whitelists_from_excel


def _write_excel_config(path, *, include_b1_complement_column: bool) -> None:
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
                'complement': [True, None, 'TRUE'],
            }
        ).to_excel(writer, sheet_name='B0', index=False)

        b1_df = pd.DataFrame({'codon': ['AGTCGA']})
        if include_b1_complement_column:
            b1_df['complement'] = [None]
        b1_df.to_excel(writer, sheet_name='B1', index=False)


def test_whitelists_from_excel_normalizes_complement_flags_for_building_blocks(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_excel_config(excel_path, include_b1_complement_column=True)

    whitelists = whitelists_from_excel(excel_path)

    assert whitelists['B0'][:3] == [
        {'index': 0, 'codon': 'GCCTCG', 'complement': True},
        {'index': 1, 'codon': 'AATTCC', 'complement': False},
        {'index': 2, 'codon': 'TTGGCC', 'complement': True},
    ]
    assert whitelists['B1'] == [{'index': 0, 'codon': 'AGTCGA', 'complement': False}]
    assert whitelists['S0'] == [
        {'name': 'sel_a', 'codon': 'AAAA'},
        {'name': 'sel_b', 'codon': 'CCCC'},
    ]
    assert whitelists['C0'] == [{'codon': 'ACAC'}]


def test_whitelists_from_excel_defaults_missing_complement_column_to_false(tmp_path):
    excel_path = tmp_path / 'config.xlsx'
    _write_excel_config(excel_path, include_b1_complement_column=False)

    whitelists = whitelists_from_excel(excel_path)

    assert whitelists['B1'] == [{'index': 0, 'codon': 'AGTCGA', 'complement': False}]
