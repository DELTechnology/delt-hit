from pathlib import Path

import pandas as pd

o_dir1 = Path(
    '/Users/adrianomartinelli/polybox - Adriano Martinelli (adriano.martinelli@pharma.ethz.ch)@polybox.ethz.ch/decl-data/raw-files/Evaluation_2301_704504_NF2-AG_yOST/')
o_dir2 = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/')

assert len(list(o_dir1.glob('selection*.txt'))) == len(list(o_dir2.glob('*_*.tsv')))

for s1_path in o_dir1.glob('selection_*.txt'):
    s1 = pd.read_csv(s1_path, sep='\t')

    s1_name = s1_path.stem
    s2_name = '_'.join(map(lambda x: str(int(x) - 1), filter(len, s1_name.split('_')[1:])))

    s2_path = o_dir2 / (s2_name + '.tsv')
    assert s2_path.exists()

    s2 = pd.read_csv(s2_path, sep='\t')
    s2 = s2.rename(columns=dict(zip(['count', 'barcode1', 'barcode2'], ['Count', 'Code1', 'Code2'])))
    s2 = s2.assign(Code1=s2.Code1 + 1, Code2=s2.Code2 + 1)

    cols = ['Count', 'Code1', 'Code2']
    s1 = s1[cols].sort_values(cols)
    s2 = s2[cols].sort_values(cols)

    if not (s1 == s2).all().all():
        print(f'Evaluation results differ for original {s1_path}.name')
