from pathlib import Path
import pandas as pd

base_dir = Path(__file__).resolve().parent
selection_name = '1_1'

latest_output = pd.read_csv(base_dir / 'lane-1' / 'selections' / selection_name / 'counts.txt', sep='\t')
latest_output = latest_output.rename(columns={'count': 'latest_c'})
latest_output.loc[:, ['code_0', 'code_1']] = latest_output.loc[:, ['code_0', 'code_1']] + 1

legacy_output = pd.read_csv(base_dir / 'published' / '1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt')
legacy_output = legacy_output[['CodeA', 'CodeB', selection_name]].rename(
    columns={'CodeA': 'code_0', 'CodeB': 'code_1', selection_name: 'legacy_c'}
)

index_cols = ['code_0', 'code_1']
latest_output = latest_output.set_index(index_cols)
legacy_output = legacy_output.set_index(index_cols)

counts = pd.concat([legacy_output, latest_output], axis=1)
counts = counts[counts.legacy_c != 0]
counts = counts.convert_dtypes()
counts['identical'] = counts.legacy_c == counts.latest_c
counts.to_csv(base_dir / f'compare_selection_{selection_name}.csv', index=True)
