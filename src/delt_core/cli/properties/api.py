from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

from delt_core.utils import read_yaml

config_path = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml')
cfg = read_yaml(config_path)

from pathlib import Path
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors as RD, QED


# ---------- 1) Create & save a dummy dataset of 50 diverse SMILES ----------
def write_dummy_smiles_csv() -> pd.DataFrame:
    smiles = [
        # aromatics / simple rings
        "c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1", "c1ccncc1",
        "c1cc2cccc2c1", "c1ccccc1O", "c1ccccc1N", "c1ccccc1C(=O)O", "c1cccc(c1)Cl",
        # small heterocycles
        "C1CCOC1", "C1COCCN1", "C1=NC=CN1", "C1=CC(=O)NC(=O)N1", "C1CCN(CC1)C",
        # common drugs-ish
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # aspirin
        "CN1C(=O)N(C)c2ncn(C)c2C1=O",  # caffeine
        "CC(=O)NC1=CC=C(O)C=C1",  # acetaminophen
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # ibuprofen
        "COC1=CC=CC=C1C(=O)OCC",  # naproxen-ish (simplified)
        # aliphatic, alcohols, acids
        "CCCC", "CCCCCC", "CCO", "CCCO", "CCCCO",
        "CC(=O)O", "CCC(=O)O", "CC(C)O", "CC(C)(C)O", "CCOC(=O)C",
        # amines / amides / nitriles
        "CCN", "CCCN", "CCCCN", "CCNC", "CC(=O)N",
        "CC#N", "CCC#N", "N#CCOC", "CC(C#N)O", "CNC(=O)C",
        # halogens & sulfur
        "CCCl", "CCBr", "CCI", "CCS", "CCSC",
        # heteroaromatics more
        "c1ncccc1", "c1ccsc1", "c1cncnc1", "c1nccs1", "c1occc1",
        # polyfunctional
        "O=C(O)CC(O)C(O)CO",  # glyceric-acid-like
        "CC(C)C(C(=O)O)N",  # valine-like
        "C(C(=O)O)N",  # glycine
        "CC(C(=O)O)O",  # lactic acid
        "COC(=O)C(O)C(O)CO"  # sugar-like ester
    ]
    # Give them simple names
    df = pd.DataFrame({
        'coda_1': range(len(smiles)),
        'coda_2': range(len(smiles)),
        'smiles': smiles})
    return df


def get_smiles_from_mol(mol):
    return Chem.MolToSmiles(mol) if mol is not None else None


# ---------- 2) Properties class ----------
class Properties:

    def run(self, *, config_path: Path):
        cfg = read_yaml(config_path)

        exp_dir = Path(cfg['experiment']['save_dir']).expanduser().resolve() / cfg['experiment']['name']
        lib_path = exp_dir / 'library.parquet'

        self.save_dir = lib_path.parent / 'properties'
        self.save_dir.mkdir(parents=True, exist_ok=True)

        self.df = pd.read_parquet(lib_path)
        self.mols = [Chem.MolFromSmiles(s) for s in self.df['smiles'].tolist()]
        self.compute_properties()

        prop_names = [col for col in self.df.columns if col.startswith('prop_')]
        for name in prop_names:
            self.plot(name=name)

    def compute_properties(self) -> None:
        self.df["prop_mw"] = [Descriptors.MolWt(m) for m in self.mols]
        self.df["prop_logP"] = [Crippen.MolLogP(m) for m in self.mols]
        self.df["prop_HBD"] = [Lipinski.NumHDonors(m) for m in self.mols]
        self.df["prop_HBA"] = [Lipinski.NumHAcceptors(m) for m in self.mols]
        self.df["prop_rotB"] = [Lipinski.NumRotatableBonds(m) for m in self.mols]
        self.df["prop_TPSA"] = [RD.CalcTPSA(m) for m in self.mols]
        self.df["prop_RBonds"] = [RD.CalcNumRotatableBonds(m) for m in self.mols]
        self.df["prop_ARings"] = [RD.CalcNumAromaticRings(m) for m in self.mols]
        self.df["prop_rings"] = [RD.CalcNumRings(m) for m in self.mols]
        self.df["prop_heavyAtoms"] = [Descriptors.HeavyAtomCount(m) for m in self.mols]
        self.df["prop_formalCharge"] = [Chem.GetFormalCharge(m) for m in self.mols]
        self.df["prop_heteroAtoms"] = [Descriptors.NumHeteroatoms(m) for m in self.mols]
        self.df["prop_fractionCsp3"] = [RD.CalcFractionCSP3(m) for m in self.mols]
        self.df["prop_QED"] = [QED.qed(m) for m in self.mols]

    def plot(self, name: str):
        if name not in self.df.columns:
            raise ValueError(f"Column {name} not found in dataframe")

        plt.figure(figsize=(8, 6))
        ax = sns.histplot(self.df[name].dropna(), kde=False, discrete=self.df[name].dtype == int)
        ax.set_title(f"Distribution of {name}")
        ax.set_xlabel(name)
        ax.set_ylabel("Frequency")
        ax.grid(True)
        ax.figure.tight_layout()

        ax.figure.savefig(self.save_dir / f"{name}.png")

    def save(self):
        self.df.to_parquet(self.save_dir / 'properties.parquet', index=False)


# ---------------- Example usage ----------------

# # 1) write dummy dataset
config_path = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml')
out_path = file_path = Path(
    "/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/library.parquet")
smiles_csv = write_dummy_smiles_csv(out_path)

# 2) compute properties
props = self = Properties()
props.run(config_path=config_path)
props.compute_properties()
props.plot(name="prop_mw")
