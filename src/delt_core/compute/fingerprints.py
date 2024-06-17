from collections import defaultdict
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import QED
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.AllChem import GetRDKitFPGenerator
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from sklearn.manifold import TSNE

from utils import read_txt


def compute_fingerprint(
        mol: Chem.rdchem.Mol,
        num_bits: int = 2048,
) -> ExplicitBitVect:
    fpgen = GetRDKitFPGenerator(fpSize=num_bits)
    return fpgen.GetFingerprint(mol)


def compute_similarity(
        bv1: ExplicitBitVect,
        bv2: ExplicitBitVect,
) -> float:
    return DataStructs.TanimotoSimilarity(bv1, bv2)


def plot_fingerprints(
        fingerprints: list,
        counts: list,
        output_file: str = 'fingerprints.png',
) -> None:
    cmap = plt.cm.YlGnBu
    plt.figure(figsize=(10, 5))
    x, y = zip(*fingerprints)
    sc = plt.scatter(x, y, c=counts, edgecolor='k', s=20.0, alpha=1.0, cmap=cmap)
    plt.colorbar(sc, label='Count')
    plt.savefig(output_file, dpi=300)


def compute_tpsa(
        mol: Chem.rdchem.Mol,
) -> float:
    return rdMolDescriptors.CalcTPSA(mol)


def compute_properties(
        mol_structures: list,
) -> dict:
    qed_properties = [QED.properties(mol) for mol in mol_structures]
    mw = [p.MW for p in qed_properties]
    alogp = [p.ALOGP for p in qed_properties]
    psa = [p.PSA for p in qed_properties]
    hbd = [p.HBD for p in qed_properties]
    hba = [p.HBA for p in qed_properties]
    rotb = [p.ROTB for p in qed_properties]
    return {
        'MW [g/mol]': mw,
        'AlogP': alogp,
        'TPSA [A^2]': psa,
        'HBDs': hbd,
        'HBAs': hba,
        'NRotBs': rotb,
    }


def plot_properties(
        properties: list,
        num_bins: list = None,
        output_file: str = 'properties.png',
) -> None:
    plt.figure(figsize=(10, 5))
    for i, ((key, values), bins) in enumerate(zip(properties.items(), num_bins), 1):
        plt.subplot(2, 3, i)
        if not bins:
            bins = np.arange(np.min(values), np.max(values) + 1) - 0.5
        plt.hist(values, bins=bins)
        plt.xlabel(key)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.savefig(output_file, dpi=300)


def main(
        input_file: Path,
        counts_file: Path,
) -> None:

    num_bb1 = 3
    num_bb2 = 548

    data = read_txt(input_file)[1:]
    data = data[:num_bb1 * num_bb2]
    table = read_txt(counts_file)[1:]
    smiles = [s.split('\t')[-2] for s in data]
    mols = [Chem.MolFromSmiles(m) for m in smiles]

    """"
    properties = compute_properties(mols)
    num_bins = [30, 30, 30, None, None, None]
    plot_properties(properties, num_bins)
    """

    counts_dict = defaultdict(int)
    for line in table:
        count, bb1, bb2 = line.strip().split('\t')
        counts_dict[f'{bb1}-{bb2}'] = int(count)
    
    fps = [compute_fingerprint(mol) for mol in mols]
    labels = [f'{bb1}-{bb2}' for bb1 in range(num_bb1) for bb2 in range(num_bb2)]
    counts = [counts_dict[label] for label in labels]

    tsne = TSNE(n_components=2, random_state=42)
    fps_2d = tsne.fit_transform(np.array(fps))
    plot_fingerprints(fps_2d, counts)



if __name__ == '__main__':

    input_file = 'smiles_nf.txt'
    counts_file = 'counts/0_0.tsv'
    main(input_file, counts_file)

