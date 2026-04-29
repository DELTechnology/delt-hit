import pandas as pd
from pathlib import Path

sel_dir = Path('/users/amarti51/projects/delt-hit/supporting_material/experiments/favalli-lane-2/selections')
sel_dir = Path('/users/amarti51/projects/delt-hit/supporting_material/experiments/favalli-lane-1/selections')
sel_path = Path('/users/amarti51/projects/delt-hit/supporting_material/experiments/favalli-lane-2/selections/AG24_40_counts.txt')
code_0, code_1, count = 3-1, 22-1, 21
code_0, code_1, count = 3-1, 23-1, 13
# code_0, code_1, count = 96-1, 239-1, 9

combs = [(3-1, 23-1, 13), (3-1, 22-1, 21)]

for sel_path in sel_dir.glob('*.txt'):
    df = pd.read_csv(sel_path, sep='\t')

    filter_ = None
    for code_0, code_1, count in combs:
        # filter_ = (df.code_0 == code_0) & (df.code_1 == code_1)
        f_ = (df.code_0 == code_0) & (df.code_1 == code_1) & (df['count'] == count)
        filter_ = f_ if filter_ is None else f_ | filter_

    tmp = df[filter_]

    if not tmp.empty:
        print(sel_path)
        print(tmp)

# AG24_204_counts.txt corresponds to 1907_NF2GB2_s1_R1_260424JSselection_57_2_.txt