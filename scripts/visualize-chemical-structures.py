from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from matplotlib import pyplot as plt

save_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-3/analyses/analysis-1')
hits_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-3/analyses/analysis-1/counts/hits.csv')
enrichment_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-3/analyses/analysis-1/edgeR/enrichment_hits.csv')
library = pd.read_parquet('/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-merged/library.parquet')
library.rename(columns = dict(code_0='code_1', code_1='code_2'), inplace=True)

def plot_chemical_structures(path: Path, library, num: int = 100):
    """
    Plots a grid of chemical structures for the top N hits from a CSV file, merging with a compound library.

    Args:
        path (Path): Path to the CSV file containing hit compounds with 'code_1' and 'code_2' columns.
        library (pd.DataFrame): DataFrame containing the compound library with 'code_1', 'code_2', and 'smiles' columns.
        num (int, optional): Number of top compounds to plot. Defaults to 100.

    Returns:
        PIL.Image.Image: The generated grid image of chemical structures.

    Side Effects:
        Saves the generated image as a PNG file in the same directory as the input CSV, with the same name but a .png extension.
    """
    hits = pd.read_csv(path)

    data = pd.merge(hits, library, on=['code_1', 'code_2'], how='inner')
    data = data[:num]
    mols = [Chem.MolFromSmiles(i) for i in data.smiles]
    legends = (data.code_1.astype(str) + '-' + data.code_2.astype(str)).tolist()

    img = Draw.MolsToGridImage(mols, molsPerRow=10, legends=legends)
    img.save(path.with_suffix('.png'))
    return img

img = plot_chemical_structures(path=hits_path, library=library, num=100)
plt.imshow(img).figure.show()

img = plot_chemical_structures(path=enrichment_path, library=library, num=100)
plt.imshow(img).figure.show()

# %% COMPARISON BETWEEN HITS AND ENRICHMENT

hits = pd.read_csv(hits_path)
enrichment = pd.read_csv(enrichment_path)
topK = 100

hit_compounds = set([tuple(i.values()) for i in hits[['code_1', 'code_2']].to_dict('records')])
enrichment_compounds = set([tuple(i.values()) for i in enrichment[['code_1', 'code_2']].to_dict('records')])
intx = hit_compounds.intersection(enrichment_compounds)

enrichment = enrichment[:len(hits)]
enrichment_compounds_topK = set([tuple(i.values()) for i in enrichment[['code_1', 'code_2']].to_dict('records')])

intx_top100 = hit_compounds.intersection(enrichment_compounds_topK)

pdat = pd.DataFrame()

df = hits.code_1.value_counts().to_frame().reset_index().rename(columns = dict(code_1='code'))
df['method'] = 'counts'
df['code_index'] = 'code_1'
pdat = pd.concat([pdat, df])

df = enrichment.code_1.value_counts().to_frame().reset_index().rename(columns = dict(code_1='code'))
df['method'] = 'edgeR'
df['code_index'] = 'code_1'
pdat = pd.concat([pdat, df])

df = hits.code_2.value_counts().to_frame().reset_index().rename(columns = dict(code_2='code'))
df['method'] = 'counts'
df['code_index'] = 'code_2'
pdat = pd.concat([pdat, df])

df = enrichment.code_2.value_counts().to_frame().reset_index().rename(columns = dict(code_2='code'))
df['method'] = 'edgeR'
df['code_index'] = 'code_2'
pdat = pd.concat([pdat, df])

pdat.reset_index(drop=True, inplace=True)
assert not pdat.isna().any().any()
pdat.to_csv(save_path / 'comparison_counts_edgeR.csv', index=False)
color_dict = dict(edgeR='tab:blue', counts='tab:orange')

import seaborn as sns
for code_index in pdat.code_index.unique():
    fig, ax = plt.subplots(figsize=(16, 4))
    sns.barplot(x='code', y='count', hue='method', data=pdat[pdat.code_index == code_index], ax=ax, palette=color_dict)
    ax.set_title(code_index)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.figure.tight_layout()
    ax.figure.savefig(save_path / f'{code_index}.png')
    ax.figure.show()
