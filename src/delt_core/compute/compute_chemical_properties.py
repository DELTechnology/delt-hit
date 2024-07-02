from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import QED

from .utils import read_txt


def compute_qed(
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
        output_file: Path,
        num_bins: list = None,
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


def compute_properties(
        input_file: Path,
) -> None:

    data = read_txt(input_file)[1:]
    smiles = [s.split('\t')[-2] for s in data]
    mols = [Chem.MolFromSmiles(m) for m in smiles]

    properties = compute_qed(mols)
    num_bins = [30, 30, 30, None, None, None]
    output_file = input_file.parent / 'properties.png'
    plot_properties(properties, output_file, num_bins)

