from pathlib import Path
import pandas as pd
import networkx as nx
from loguru import logger
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from delt_core.demultiplex.parser import catalog_from_excel
from delt_core.utils import read_yaml

import pickle
from itertools import batched
from pathlib import Path

import numpy as np
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy import sparse
from tqdm import tqdm

class Represent:

    def run(self, *, config_path: Path, method: str = 'morgan'):
        cfg = read_yaml(config_path)

        exp_dir = Path(cfg['experiment']['save_dir']).expanduser().resolve() / cfg['experiment']['name']
        save_dir = exp_dir / 'representations'
        save_dir.mkdir(parents=True, exist_ok=True)

        lib_path = exp_dir / 'library.parquet'
        df = pd.read_parquet(lib_path)
        smiles = df.smiles

        match method:
            case 'morgan':
                run_morgan(smiles, save_path =save_dir / 'morgan_fps.npz')

    def deep(self, *, model_name: str, path: Path, save_path: Path, device='cuda'):

        df = pd.read_parquet(path)
        smiles = df.smiles.tolist()

        if model_name == 'bert':
            fps = get_bert_fp(smiles, device=device)
            fps = np.vstack(fps)
        else:
            raise ValueError(f'Unknown model name: {model_name}')

        save_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(save_path, fps)

        logger.info(f"Representations saved to {save_path}")

def get_bert_fp(smiles: list[str], device='cuda'):
    from transformers import BertTokenizerFast, BertModel
    import torch

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


def run_morgan(smiles: list[str], save_path: Path):
    fps = []
    for smiles in tqdm(smiles):
        fp = get_morgan_fp(smiles)
        fps.append(fp)

    fps = sparse.vstack(fps, format="csr")
    save_path.parent.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(save_path, fps)

    logger.info(f"Fingerprints saved to {save_path}")

def get_morgan_fp(smiles, radius=2, n_bits=2048) -> sparse.csr_array:

    if smiles is None:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mfpgen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = mfpgen.GetFingerprint(mol)
    fp = sparse.csr_array(fp, dtype=np.uint8)
    return fp