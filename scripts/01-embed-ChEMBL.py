from itertools import batched
from pathlib import Path

import numpy as np
import pandas as pd
import umap.plot
from ai4bmr_datasets import ChEMBL
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

chembl_dir = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/embeddings/ChEMBL')
chembl_dir.mkdir(parents=True, exist_ok=True)

nf2_dir = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/embeddings/NF2')
nf2_dir.mkdir(parents=True, exist_ok=True)

# %% DATA
logger.info('Loading ChEMBL')

query = """SELECT
    md.chembl_id,
    cs.canonical_smiles,
    md.molecule_type,
    di.efo_term,
    di.mesh_heading
FROM
    molecule_dictionary md
LEFT JOIN
    compound_structures cs ON md.molregno = cs.molregno
LEFT JOIN
    drug_indication di ON md.molregno = di.molregno
WHERE
    md.molecule_type IS NOT NULL
"""
# LIMIT 500000;"""

ds = ChEMBL()
ds.setup(query=query)
smiles = ds[0]['canonical_smiles']

data = ds.db.drop_duplicates('chembl_id')
data = data.set_index('chembl_id')

# %% MORGAN
def get_morgan_fp(smiles, radius=2, n_bits=2048):

    if smiles is None:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mfpgen = AllChem.GetMorganGenerator(radius=radius,fpSize=n_bits)
    fp = mfpgen.GetFingerprint(mol)
    fp = np.asarray(fp)
    return fp

smiles = data.canonical_smiles
save_path = chembl_dir / 'morgan_fp.parquet'
if save_path.exists():
    logger.info(f'Loading ChEMBL morgan_fp from cache: {save_path}')
    morgan_fp = pd.read_parquet(save_path, engine='fastparquet')
else:
    logger.info(f'Computing ChEMBL morgan_fp...')
    morgan_fp = smiles.map(get_morgan_fp)
    logger.info(f'Caching ChEMBL morgan_fp...')
    filter_ = morgan_fp.notna()
    morgan_fp = morgan_fp[filter_]
    morgan_fp = pd.DataFrame(np.stack(morgan_fp.tolist()), index=morgan_fp.index)
    morgan_fp.columns = morgan_fp.columns.astype(str)
    morgan_fp.to_parquet(save_path, engine='fastparquet')

# %% BERT
from transformers import BertTokenizerFast, BertModel
import torch

def get_bert_fp(smiles: list[str]):
    checkpoint = 'unikei/bert-base-smiles'
    tokenizer = BertTokenizerFast.from_pretrained(checkpoint)
    model = BertModel.from_pretrained(checkpoint)
    model.to('cuda')

    bert_fp = []
    batch_size = 512
    for batch in tqdm(batched(smiles, batch_size), total= len(smiles) // batch_size + 1):
        tokens = tokenizer(batch, return_tensors='pt', padding=True, truncation=True, max_length=512)
        tokens = {k: v.to('cuda') for k, v in tokens.items()}
        with torch.no_grad():
            predictions = model(**tokens)
        bert_fp.append(predictions.pooler_output.cpu().numpy())

    return bert_fp

smiles = data.loc[morgan_fp.index].canonical_smiles

smiles = data.canonical_smiles
filter_ = smiles.notna()
smiles = smiles[filter_]

save_path = chembl_dir / 'bert_fp.parquet'
if save_path.exists():
    logger.info(f'Loading ChEMBL bert_fp from cache: {save_path}')
    bert_fp = pd.read_parquet(save_path, engine='fastparquet')
else:
    logger.info(f'Computing ChEMBL bert_fp...')
    bert_fp = get_bert_fp(smiles.tolist())
    bert_fp = pd.DataFrame(np.vstack(bert_fp), index=morgan_fp.index)
    bert_fp.columns = bert_fp.columns.astype(str)
    bert_fp.to_parquet(save_path, engine='fastparquet')

# %% NF2
import pandas as pd
path = "/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/libraries/smiles/NF_smiles.txt.gz"
df = pd.read_csv(path, compression='gzip', sep='\t')

df[['Scaffold_L1', 'Product_L1']]

smiles = df.Product_L1
save_path = nf2_dir / 'morgan_fp.parquet'
if save_path.exists():
    logger.info(f'Loading NF2 morgan_fp from cache: {save_path}')
    morgan_fp = pd.read_parquet(save_path, engine='fastparquet')
else:
    logger.info(f'Computing NF2 morgan_fp')
    morgan_fp = smiles.map(get_morgan_fp)
    filter_ = morgan_fp.notna()
    morgan_fp = morgan_fp[filter_]
    morgan_fp = pd.DataFrame(np.stack(morgan_fp.tolist()), index=morgan_fp.index)
    morgan_fp.columns = morgan_fp.columns.astype(str)
    morgan_fp.to_parquet(save_path, engine='fastparquet')

save_path = nf2_dir / 'bert_fp.parquet'
if save_path.exists():
    logger.info(f'Loading NF2 bert_fp from cache: {save_path}')
    bert_fp = pd.read_parquet(save_path, engine='fastparquet')
else:
    logger.info(f'Computing NF2 bert_fp')
    bert_fp = get_bert_fp(smiles.tolist())
    bert_fp = pd.DataFrame(np.vstack(bert_fp))
    bert_fp.columns = bert_fp.columns.astype(str)
    bert_fp.to_parquet(save_path, engine='fastparquet')

# %%
num_obs = 100_000
save_dir = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/embeddings')
# for fp in ['morgan_fp', 'bert_fp']:
for fp in ['morgan_fp']:
    logger.info(f'Computing UMAP for {fp}')

    chembl_fp = pd.read_parquet(chembl_dir / f'{fp}.parquet', engine='fastparquet')
    chembl_fp['lib'] = 'chembl'
    chembl_fp = chembl_fp.sample(num_obs)

    nf2_fp = pd.read_parquet(nf2_dir / f'{fp}.parquet', engine='fastparquet')
    nf2_fp['lib'] = 'nf2'
    nf2_fp = nf2_fp.sample(num_obs)

    pdat = pd.concat([chembl_fp, nf2_fp], axis=0)
    labels =  pdat.pop('lib')

    mapper = umap.UMAP(metric='jaccard').fit(pdat)
    ax = umap.plot.points(mapper, labels=labels, background='black')
    ax.figure.show()
    ax.figure.savefig(save_dir / f'{fp}-num_obs={num_obs}.png', dpi=300)


