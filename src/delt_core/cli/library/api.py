from itertools import batched
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors as RD, QED
from rdkit.Chem import rdChemReactions
from scipy import sparse
from tqdm import tqdm

from delt_core.utils import read_yaml


class Library:

    def get_experiment_dir(self, *, config_path: Path) -> Path:
        cfg = read_yaml(config_path)
        exp_dir = Path(cfg['experiment']['save_dir']).expanduser().resolve() / cfg['experiment']['name']
        return exp_dir

    def get_library_path(self, *, config_path: Path) -> Path:
        exp_dir = self.get_experiment_dir(config_path=config_path)
        lib_path = exp_dir / 'library.parquet'
        return lib_path

    def enumerate(self, *, config_path: Path):
        lib_path = self.get_library_path(config_path=config_path)
        if lib_path.exists():
            logger.info(f'Library {lib_path} exists')
            return

        config_path = Path(
            '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config_with_steps.yaml')
        cfg = read_yaml(config_path)

        steps = cfg['library']['steps']
        catalog = cfg['catalog']
        reactions = catalog['reactions']
        compounds = catalog['compounds']
        products = {k: dict(smiles=None) for k in cfg['library']['products']['names']}
        terminal = cfg['library']['products']['terminal']

        G = get_reaction_graph(steps=steps, reactions=reactions, compounds=compounds, products=products)
        ax = visualize_reaction_graph(G)
        ax.figure.savefig(lib_path.parent / 'reaction_graph.png', dpi=300)

        G = complete_reaction_graph(G)
        product = G.nodes[terminal]['smiles']
        ax = visualize_smiles([product])
        ax.figure.savefig(lib_path.parent / 'products.png', dpi=300)

        df = get_dummy_library()
        df.to_parquet(lib_path, index=False)

    def properties(self, *, config_path: Path):
        lib_path = self.get_library_path(config_path=config_path)

        save_dir = lib_path.parent / 'properties'
        save_dir.mkdir(parents=True, exist_ok=True)

        df = pd.read_parquet(lib_path)
        df = self.compute_properties(data=df)
        df.to_parquet(save_dir / 'properties.parquet', index=False)

        prop_names = [col for col in df.columns if col.startswith('prop_')]
        for name in prop_names:
            ax = self.plot_property(data=df, name=name)
            ax.figure.savefig(save_dir / f"{name}.png")
            plt.close(ax.figure)

    def compute_properties(self, data: pd.DataFrame) -> pd.DataFrame:
        mols = [Chem.MolFromSmiles(s) for s in data['smiles'].tolist()]
        data["prop_mw"] = [Descriptors.MolWt(m) for m in mols]
        data["prop_logP"] = [Crippen.MolLogP(m) for m in mols]
        data["prop_HBD"] = [Lipinski.NumHDonors(m) for m in mols]
        data["prop_HBA"] = [Lipinski.NumHAcceptors(m) for m in mols]
        data["prop_rotB"] = [Lipinski.NumRotatableBonds(m) for m in mols]
        data["prop_TPSA"] = [RD.CalcTPSA(m) for m in mols]
        data["prop_RBonds"] = [RD.CalcNumRotatableBonds(m) for m in mols]
        data["prop_ARings"] = [RD.CalcNumAromaticRings(m) for m in mols]
        data["prop_rings"] = [RD.CalcNumRings(m) for m in mols]
        data["prop_heavyAtoms"] = [Descriptors.HeavyAtomCount(m) for m in mols]
        data["prop_formalCharge"] = [Chem.GetFormalCharge(m) for m in mols]
        data["prop_heteroAtoms"] = [Descriptors.NumHeteroatoms(m) for m in mols]
        data["prop_fractionCsp3"] = [RD.CalcFractionCSP3(m) for m in mols]
        data["prop_QED"] = [QED.qed(m) for m in mols]
        return data

    def plot_property(self, data: pd.DataFrame, name: str) -> plt.Axes:
        if name not in data.columns:
            raise ValueError(f"Column {name} not found in dataframe")

        ax = sns.histplot(data[name].dropna(), kde=False, discrete=data[name].dtype == int)
        ax.set_title(f"Distribution of {name}")
        ax.set_xlabel(name)
        ax.set_ylabel("Frequency")
        ax.grid(True)
        ax.figure.tight_layout()
        return ax

    def represent(self, *, config_path: Path, method: str = 'morgan'):
        exp_dir = self.get_experiment_dir(config_path=config_path)

        save_dir = exp_dir / 'representations'
        save_dir.mkdir(parents=True, exist_ok=True)

        lib_path = exp_dir / 'library.parquet'
        df = pd.read_parquet(lib_path)
        smiles = df.smiles

        match method:
            case 'morgan':
                run_morgan(smiles, save_path=save_dir / 'morgan.npz')
            case 'bert':
                run_morgan(smiles, save_path=save_dir / 'bert.npz')


def run_bert(*, model_name: str, path: Path, save_path: Path, device='cuda'):
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
    for batch in tqdm(batched(smiles, batch_size), total=len(smiles) // batch_size + 1):
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
    mol = Chem.MolFromSmiles(smiles)
    mfpgen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = mfpgen.GetFingerprint(mol)
    fp = sparse.csr_array(fp, dtype=np.uint8)
    return fp


def get_dummy_library() -> pd.DataFrame:
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


def get_reaction_graph(steps: list, reactions: dict, compounds: dict, products: dict) -> nx.DiGraph:
    G = nx.DiGraph()

    G.add_edges_from(steps)

    attrs = {k: {**v, 'type': 'reaction'} for k, v in reactions.items()}
    nx.set_node_attributes(G, attrs)

    attrs = {k: {**v, 'type': 'compound'} for k, v in compounds.items()}
    nx.set_node_attributes(G, attrs)

    attrs = {k: {**v, 'type': 'product'} for k, v in products.items()}
    nx.set_node_attributes(G, attrs)

    return G


def visualize_reaction_graph(G: nx.DiGraph) -> plt.Axes:
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    compounds = [n for n, d in G.nodes(data=True) if d.get("type") == "compound"]
    reactions = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    products = [n for n, d in G.nodes(data=True) if d.get("type") == "product"]

    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")

    # Draw compounds (blue circles)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=compounds,
        node_color="lightblue",
        node_shape="o",  # circle
        node_size=500,
        ax=ax
    )

    # Draw compounds (blue circles)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=products,
        node_color="salmon",
        node_shape="o",  # circle
        node_size=500,
        ax=ax
    )

    # Draw reactions (green squares)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=reactions,
        node_color="lightgreen",
        node_shape="s",  # square
        node_size=600,
        ax=ax
    )

    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)
    nx.draw_networkx_edges(G, pos, ax=ax, arrows=True)

    fig.show()
    return ax


def find_next_reaction(G: nx.DiGraph):
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    for node in reaction_nodes:
        preds = sorted(G.predecessors(node))
        succ, = sorted(G.successors(node))  # note: currently only one product per reaction

        if G.nodes[succ]['smiles'] is None and all([G.nodes[i]['smiles'] is not None for i in preds]):
            return {'reactants': preds, 'reaction': node, 'product': succ}

    return None


def perform_reaction(smirks: str, reactants: list[str]) -> list[str]:
    mols = [Chem.MolFromSmiles(i) for i in reactants]
    rxn = rdChemReactions.ReactionFromSmarts(smirks, useSmiles=True)

    # note: order matters
    product_sets = rxn.RunReactants([*mols])

    products = set()
    for tup in product_sets:
        for pmol in tup:
            Chem.SanitizeMol(pmol)
            products.add(Chem.MolToSmiles(pmol, canonical=True, kekuleSmiles=False, isomericSmiles=False))

    return sorted(products)


def complete_reaction_graph(G: nx.DiGraph) -> nx.DiGraph:
    while True:

        try:
            next_reaction = find_next_reaction(G)
            if next_reaction is None:
                break

            smirks = G.nodes[next_reaction['reaction']]['smirks']
            reactants = [G.nodes[i]['smiles'] for i in next_reaction['reactants']]
            products = perform_reaction(smirks, reactants)

            # Take the first product if multiple are formed
            product = {next_reaction['product']: dict(smiles=products[0])}
            nx.set_node_attributes(G, product)
            print(f"Reaction {next_reaction['reaction']}: {reactants} -> {products[0]}")
        except Exception as e:
            print(f"Error processing reaction {next_reaction}: {e}")
            break

    return G

def visualize_smiles(smiles: list[str], nrow: int = 25):
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    # could provide legends=product_names
    nrow = min(nrow, len(mols))
    img = Draw.MolsToGridImage(mols, molsPerRow=nrow, subImgSize=(200, 200))
    plt.figure(figsize=(10, 6))
    ax = plt.imshow(img)
    ax.axes.set_axis_off()
    ax.axes.set_title('Product Structures')
    ax.figure.tight_layout()
    return ax
