from pathlib import Path
import pandas as pd

delt_204 = pd.read_csv('/users/amarti51/projects/delt-hit/supporting_material/experiments/favalli-lane-1/selections/AG24_204_counts.txt', sep='\t')
delt_204 = delt_204.rename(columns={'count': 'delt_hit'})
delt_204.loc[:, ['code_0', 'code_1']] = delt_204.loc[:, ['code_0', 'code_1']] + 1

new_57_2 = pd.read_csv('/users/amarti51/projects/delt-hit/supporting_material/experiments/favalli/altesC++/1907_NF2GB2_s1_R1_260424JSselection_57_2_.txt', sep='\t')
new_57_2 = new_57_2.rename(columns={'Code1': 'code_0', 'Code2': 'code_1', 'Count': 'new_c'})

old_c = pd.read_csv('/users/amarti51/projects/delt-hit/supporting_material/experiments/favalli/altesC++/1907_NF2GB2_s1_R1_260424JS_2026_4_24_16_20_51_eval.txt')

old_57_2 = old_c[['CodeA', 'CodeB', '57_2']].rename(columns={'CodeA': 'code_0', 'CodeB': 'code_1', '57_2': 'old_c'})

index_cols = ['code_0', 'code_1']
delt_204 = delt_204.set_index(index_cols)
new_57_2 = new_57_2.set_index(index_cols)
old_57_2 = old_57_2.set_index(index_cols)

counts = pd.concat([old_57_2, new_57_2, delt_204], axis=1)
counts = counts[counts.old_c != 0]
counts = counts.convert_dtypes()
counts['identical'] = (counts.old_c == counts.new_c) & (counts.new_c == counts.delt_hit)
counts.to_csv('/users/amarti51/projects/delt-hit/supporting_material/experiments/favalli/compare_selection_57_2.csv', index=True)
