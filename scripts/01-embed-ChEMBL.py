import pickle
from itertools import batched
from pathlib import Path

import numpy as np
import umap.plot
from ai4bmr_datasets import ChEMBL
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy import sparse
from tqdm import tqdm

# base_dir = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/datasets/ChEMBL')
base_dir = None
# data_dir = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB')
data_dir = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/DECLT-DB/')

chembl_dir = data_dir / 'embeddings' / 'ChEMBL'
chembl_dir.mkdir(parents=True, exist_ok=True)

nf2_dir = data_dir / 'embeddings' / 'NF2'
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
LIMIT 5000;"""
# """

ds = ChEMBL(base_dir=base_dir)
ds.prepare_data()
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
    fp = sparse.csr_array(fp, dtype=np.uint8)
    return fp

smiles = data.canonical_smiles
save_path = chembl_dir / 'morgan_fp.npz'
index_path = chembl_dir / 'index.pkl'
if save_path.exists():
    logger.info(f'Loading ChEMBL morgan_fp from cache: {save_path}')
    morgan_fp = sparse.load_npz(save_path)
    with open(index_path, 'rb') as f:
        index = pickle.load(f)
else:
    logger.info(f'Computing ChEMBL morgan_fp...')
    morgan_fps = [get_morgan_fp(s) for s in tqdm(smiles, total=len(smiles))]

    filter_ = [fp is not None for fp in morgan_fps]
    index = smiles.index[filter_].tolist()
    morgan_fps = [fp for fp in morgan_fps if fp is not None]
    morgan_fps = sparse.vstack(morgan_fps)
    sparse.save_npz(save_path, morgan_fps)
    with open(index_path, 'wb') as f:
        pickle.dump(index, f)

# %% BERT
from transformers import BertTokenizerFast, BertModel
import torch

def get_bert_fp(smiles: list[str], device='cuda'):
    checkpoint = 'unikei/bert-base-smiles'
    tokenizer = BertTokenizerFast.from_pretrained(checkpoint)
    model = BertModel.from_pretrained(checkpoint)
    model.to(device)

    bert_fp = []
    batch_size = 128
    for batch in tqdm(batched(smiles, batch_size), total= len(smiles) // batch_size + 1):
        tokens = tokenizer(batch, return_tensors='pt', padding=True, truncation=True, max_length=512)
        tokens = {k: v.to(device) for k, v in tokens.items()}
        with torch.no_grad():
            predictions = model(**tokens)
        bert_fp.append(predictions.pooler_output.cpu().numpy())

    return bert_fp

smiles = data.loc[index].canonical_smiles

save_path = chembl_dir / 'bert_fp.npy'
device = 'mps'
if save_path.exists():
    logger.info(f'Loading ChEMBL bert_fp from cache: {save_path}')
    bert_fp = np.load(save_path, allow_pickle=True)
else:
    logger.info(f'Computing ChEMBL bert_fp...')
    bert_fp = get_bert_fp(smiles.tolist(), device=device)
    bert_fp = np.vstack(bert_fp)
    np.save(save_path, bert_fp)

# %% NF2
import pandas as pd
path = data_dir / "libraries/smiles/NF_smiles.txt.gz"
df = pd.read_csv(path, compression='gzip', sep='\t')

df[['Scaffold_L1', 'Product_L1']]

smiles = df.Product_L1[:100]
save_path = nf2_dir / 'morgan_fp.npz'
index_path = nf2_dir / 'index.pkl'
if save_path.exists():
    logger.info(f'Loading NF2 morgan_fp from cache: {save_path}')
    morgan_fp = sparse.load_npz(save_path)
    with open(index_path, 'rb') as f:
        index = pickle.load(f)
else:
    logger.info(f'Computing NF2 morgan_fp')
    morgan_fps = [get_morgan_fp(s) for s in tqdm(smiles, total=len(smiles))]

    filter_ = [fp is not None for fp in morgan_fps]
    index = smiles.index[filter_].tolist()
    morgan_fps = [fp for fp in morgan_fps if fp is not None]
    morgan_fps = sparse.vstack(morgan_fps)

    sparse.save_npz(save_path, morgan_fps)
    with open(index_path, 'wb') as f:
        pickle.dump(index, f)

save_path = nf2_dir / 'bert_fp.npy'
if save_path.exists():
    logger.info(f'Loading NF2 bert_fp from cache: {save_path}')
    bert_fp = np.load(save_path, allow_pickle=True)
else:
    logger.info(f'Computing NF2 bert_fp')
    bert_fp = get_bert_fp(smiles.tolist(), device=device)
    bert_fp = np.vstack(bert_fp)
    np.save(save_path, bert_fp)

# %%
num_obs = 100_000
save_dir = data_dir / 'embeddings'
rng = np.random.default_rng(42)
is_sparse = False
# for fp in ['morgan_fp']:
for fp in ['morgan_fp.npz', 'bert_fp.npy']:
    logger.info(f'Computing UMAP for {fp}')

    is_sparse = 'npz' in fp

    chembl_fps = sparse.load_npz(chembl_dir / fp) if is_sparse else np.load(chembl_dir / fp, allow_pickle=True)
    n = chembl_fps.shape[0]
    idc = rng.choice(n, size=num_obs, replace=False)
    chembl_fps = chembl_fps[idc]
    labels = ['ChEMBL'] * num_obs

    nf2_fps = sparse.load_npz(nf2_dir / fp) if is_sparse else np.load(nf2_dir / fp, allow_pickle=True)
    n = nf2_fps.shape[0]
    idc = rng.choice(n, size=num_obs, replace=False)
    nf2_fps = nf2_fps[idc]
    labels += ['NF2'] * num_obs

    pdat = sparse.vstack([chembl_fps, nf2_fps]) if is_sparse else np.vstack([chembl_fps, nf2_fps])
    pdat = pdat.toarray() if is_sparse else pdat
    labels = np.array(labels)

    metric = 'jaccard' if fp == 'morgan_fp' else 'cosine'
    mapper = umap.UMAP(metric=metric).fit(pdat)
    ax = umap.plot.points(mapper, labels=labels, background='black')
    ax.set_title(f'{fp.split('.')[0]} - num_obs={num_obs}', fontsize=16)
    ax.figure.show()
    ax.figure.savefig(save_dir / f'{fp}-num_obs={num_obs}.png', dpi=300)


