from itertools import batched
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors as RD, QED
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from scipy import sparse
from tqdm import tqdm

from delt_hit.utils import read_yaml

class Library:

    def get_experiment_dir(self, *, config_path: Path) -> Path:
        cfg = read_yaml(config_path)
        exp_dir = Path(cfg['experiment']['save_dir']).expanduser().resolve() / cfg['experiment']['name']
        return exp_dir

    def get_library_path(self, *, config_path: Path) -> Path:
        exp_dir = self.get_experiment_dir(config_path=config_path)
        lib_path = exp_dir / 'library' / 'library.parquet'
        return lib_path

    def enumerate(self, *, config_path: Path, debug: str = 'False', overwrite: bool = False,
                  graph_only: bool = False, errors: str = 'raise', building_block_ids: list[str] | None = None):

        lib_path = self.get_library_path(config_path=config_path)
        if lib_path.exists() and not overwrite:
            logger.info(f'Library {lib_path} exists')
            return

        cfg = read_yaml(config_path)

        building_block_edges = cfg['library']['bb_edges']
        other_edges = cfg['library']['other_edges']
        steps = building_block_edges + other_edges
        building_blocks = {k: dict(smiles=None) for k in cfg['library']['building_blocks']}
        products = {k: dict(smiles=None) for k in cfg['library']['products']}

        catalog = cfg['catalog']
        reactions = catalog['reactions']
        compounds = catalog['compounds']

        bb_G = get_reaction_graph(steps=building_block_edges,
                                  reactions=reactions,
                                  building_blocks=building_blocks,
                                  compounds=compounds,
                                  products=products)

        ax = visualize_reaction_graph(bb_G)
        # ax.figure.show()
        lib_path.parent.mkdir(parents=True, exist_ok=True)
        ax.figure.savefig(lib_path.parent / 'building_block_reactions_graph.png', dpi=300)

        add_G = get_reaction_graph(steps=other_edges,
                                   reactions=reactions,
                                   building_blocks=building_blocks,
                                   compounds=compounds,
                                   products=products)
        ax = visualize_reaction_graph(add_G)
        # ax.figure.show()
        ax.figure.savefig(lib_path.parent / 'additional_reactions_graph.png', dpi=300)

        G = get_reaction_graph(steps=steps,
                               reactions=reactions,
                               building_blocks=building_blocks,
                               compounds=compounds,
                               products=products)
        ax = visualize_reaction_graph(G)
        ax.figure.savefig(lib_path.parent / 'reaction_graph.png', dpi=300)

        logger.info(f'Saved reaction graph visualizations to {lib_path.parent}')

        if graph_only:
            return

        building_block_names = sorted(building_blocks)
        if building_block_ids:
            building_block_names = list(filter(lambda x: x in building_block_ids, building_block_names))
        lists = [cfg['whitelists'][bbn] for bbn in building_block_names]
        combs = list(product(*lists))

        library = []
        logger.info(f'Starting enumeration of library...')
        for i, comb in tqdm(enumerate(combs)):

            bb_edges = [(bb, c['reaction']) for bb, c in zip(building_block_names, comb)]
            bb_edges += [(c['reaction'], c['product']) for c in comb]
            bb_edges += [(c['educt'], c['reaction']) for c in comb]
            bb_nodes = set([n for e in bb_edges for n in e])

            # NOTE: we add all subgraphs from the additional reaction if a component in the building block
            #   reaction graph is a sink (out_degree == 0) in the subgraph. This means the additional reactions
            #   contain instructions on how to synthesize that component.
            additional_edges = set()
            subgraphs = list(nx.weakly_connected_components(add_G))
            for n in bb_nodes:
                for sgn in subgraphs:
                    sg = add_G.subgraph(sgn).copy()
                    if n in sg and sg.out_degree(n) == 0:
                        additional_edges.update(sg.edges)

            additional_edges = [tuple(e) for e in additional_edges]

            edges = bb_edges + additional_edges
            edges = [tuple(e) for e in edges]
            nodes = [n for e in edges for n in e]

            rnx = set([n for n in nodes if n in reactions])
            prods = set([n for n in nodes if n in products])
            comps = set([n for n in nodes if n in compounds])

            rnx = {r: cfg['catalog']['reactions'][r] for r in rnx}
            prods = {p: dict(smiles=None) for p in prods}
            comps = {c: cfg['catalog']['compounds'][c] for c in comps}
            bbs = {bbn: bb
                   for bbn, bb in zip(building_block_names, comb)
                   if not pd.isna(bb['smiles'])}

            nodes = {**comps, **prods, **bbs, **rnx}

            g = G.edge_subgraph(edges=edges).copy()
            sinks = [n for n, d in g.out_degree() if d == 0]
            is_valid = len(sinks) == 1

            if (debug == 'all') or (debug == 'valid') and is_valid:
                ax = visualize_reaction_graph(g)
                ax.figure.savefig(
                    lib_path.parent / f'reaction_graph_combination={i}_{"_".join(str(c["index"]) for c in comb)}.png',
                    dpi=300)
                plt.close('all')
                # ax.figure.show()

            if is_valid:
                assert len(sinks) == 1, f"Expected exactly one sink node, found {len(sinks)}"
                terminal = sinks[0]
            else:
                logger.warning(
                    f'More than one terminal node for combination: {i}',
                    f"Run with debug='all' to visualize reaction graphs with multiple terminal nodes.")
                # NOTE: this means the combination is invalid,
                #   the participating reactions and compounds are not linearly connected. This happens for example if a
                #   certain product_1 is only used in a subset of subsequent reactions.
                continue

            nx.set_node_attributes(g, nodes)

            try:
                g = complete_reaction_graph(g, errors=errors)
            except Exception as e:
                if debug == 'invalid':
                    ax = visualize_reaction_graph(g)
                    ax.figure.savefig(
                        lib_path.parent / f'reaction_graph_combination={i}_{"_".join(str(c["index"]) for c in comb)}.png',
                        dpi=300)
                    plt.close('all')

                if errors == 'raise':
                    raise e
                elif errors == 'ignore':
                    continue

            smiles = g.nodes[terminal]['smiles']
            record = {f'code_{i}': c['index'] for i, c in enumerate(comb)}
            record['smiles'] = smiles
            library.append(record)

        df = pd.DataFrame(library)
        df = df[df.smiles.notna()]
        # df = get_dummy_library()
        df.to_parquet(lib_path, index=False)

    def properties(self, *, config_path: Path, library_path: Path | None = None):
        lib_path = library_path or self.get_library_path(config_path=config_path)

        save_dir = lib_path.parent / 'properties'
        save_dir.mkdir(parents=True, exist_ok=True)

        df = pd.read_parquet(lib_path)
        df = self.compute_properties(data=df)
        df.to_parquet(save_dir / 'properties.parquet', index=False)

        prop_names = [col for col in df.columns if col.startswith('prop_')]
        plt.close('all')
        for name in prop_names:
            ax = self.plot_property(data=df, name=name)
            ax.figure.savefig(save_dir / f"{name}.png")
            plt.close(ax.figure)

    def compute_properties(self, data: pd.DataFrame) -> pd.DataFrame:

        records = []
        for smiles in tqdm(data['smiles']):
            record = {}
            m = Chem.MolFromSmiles(smiles)
            record["prop_mw"] = Descriptors.MolWt(m)
            record["prop_logP"] = Crippen.MolLogP(m)
            record["prop_HBD"] = Lipinski.NumHDonors(m)
            record["prop_HBA"] = Lipinski.NumHAcceptors(m)
            record["prop_rotB"] = Lipinski.NumRotatableBonds(m)
            record["prop_TPSA"] = RD.CalcTPSA(m)
            record["prop_RBonds"] = RD.CalcNumRotatableBonds(m)
            record["prop_ARings"] = RD.CalcNumAromaticRings(m)
            record["prop_rings"] = RD.CalcNumRings(m)
            record["prop_heavyAtoms"] = Descriptors.HeavyAtomCount(m)
            record["prop_formalCharge"] = Chem.GetFormalCharge(m)
            record["prop_heteroAtoms"] = Descriptors.NumHeteroatoms(m)
            record["prop_fractionCsp3"] = RD.CalcFractionCSP3(m)
            record["prop_QED"] = QED.qed(m)
            records.append(record)

        props = pd.DataFrame(records)
        props = pd.concat([data, props], axis=1)
        return props

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

    def represent(self, *, config_path: Path, method: str = 'morgan', library_path: Path | None = None):
        exp_dir = self.get_experiment_dir(config_path=config_path)

        save_dir = exp_dir / 'representations'
        save_dir.mkdir(parents=True, exist_ok=True)

        lib_path = library_path or self.get_library_path(config_path=config_path)
        df = pd.read_parquet(lib_path)
        smiles = df.smiles

        match method:
            case 'morgan':
                run_morgan(smiles, save_path=save_dir / 'morgan.npz')
            case 'bert':
                run_morgan(smiles, save_path=save_dir / 'bert.npz')


# self = Library()


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


def get_reaction_graph(steps: list,
                       reactions: dict,
                       compounds: dict,
                       products: dict,
                       building_blocks: dict = None) -> nx.DiGraph:
    G = nx.DiGraph()

    G.add_edges_from(steps)

    attrs = {k: {**v, 'type': 'reaction'} for k, v in reactions.items()}
    nx.set_node_attributes(G, attrs)

    attrs = {k: {**v, 'type': 'compound'} for k, v in compounds.items()}
    nx.set_node_attributes(G, attrs)

    attrs = {k: {**v, 'type': 'product'} for k, v in products.items()}
    nx.set_node_attributes(G, attrs)

    if building_blocks is not None:
        attrs = {k: {**v, 'type': 'building_block'} for k, v in building_blocks.items()}
        nx.set_node_attributes(G, attrs)

    return G


def visualize_reaction_graph(G: nx.DiGraph) -> plt.Axes:
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    compounds = [n for n, d in G.nodes(data=True) if d.get("type") == "compound"]
    reactions = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    products = [n for n, d in G.nodes(data=True) if d.get("type") == "product"]
    building_blocks = [n for n, d in G.nodes(data=True) if d.get("type") == "building_block"]

    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")

    # Draw compounds
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=compounds,
        node_color="lightblue",
        node_shape="o",  # circle
        node_size=500,
        ax=ax
    )

    # Draw building blocks
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=building_blocks,
        node_color="mediumorchid",
        node_shape="o",  # circle
        node_size=500,
        ax=ax
    )

    # Draw compounds
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=products,
        node_color="salmon",
        node_shape="o",  # circle
        node_size=500,
        ax=ax
    )

    # Draw reactions
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

    # fig.show()
    return ax


def find_next_reaction(G: nx.DiGraph):
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    for node in reaction_nodes:
        preds = sorted(G.predecessors(node))
        succ, = sorted(G.successors(node))  # note: currently only one product per reaction

        if G.nodes[succ]['smiles'] is None and all([G.nodes[i]['smiles'] is not None for i in preds]):
            return {'reactants': preds, 'reaction': node, 'product': succ}

    return None


def perform_reaction(smirks: str, reactants: list[str], use_smiles: bool = False) -> list[str]:
    mols = [Chem.MolFromSmiles(i) for i in reactants]
    rxn = rdChemReactions.ReactionFromSmarts(smirks, useSmiles=use_smiles)

    # note: order matters
    product_sets = rxn.RunReactants([*mols])

    products = set()
    for tup in product_sets:
        for pmol in tup:
            Chem.SanitizeMol(pmol)
            products.add(Chem.MolToSmiles(pmol, canonical=True, kekuleSmiles=False, isomericSmiles=False))

    return sorted(products)


def complete_reaction_graph(G: nx.DiGraph, errors: str = 'raise') -> nx.DiGraph:
    while True:

        try:
            next_reaction = find_next_reaction(G)
            if next_reaction is None:
                break

            reactants = [G.nodes[i]['smiles'] for i in next_reaction['reactants']]
            smirks = G.nodes[next_reaction['reaction']]['smirks']

            if pd.isna(smirks):
                assert len(reactants) == 1, "PASS reaction should have exactly one reactant"
                products = reactants
            else:
                products = perform_reaction(smirks, reactants)

            if len(products) == 0:
                products = perform_reaction(smirks, reactants[::-1])

            assert len(products) == 1, f"Expected exactly one product, found {len(products)} for reaction {next_reaction}"
            product = {next_reaction['product']: dict(smiles=products[0])}
            nx.set_node_attributes(G, product)

        except Exception as e:

            logger.error(f"Error processing reaction {next_reaction}: {e}\n")
            logger.error(f"Current products: {products}\n")
            logger.error(f'Reaction graph at  error: {G.nodes(data=True)}\n')
            data = G.nodes(data=True)

            data = [i for i in G.nodes(data=True) if i[0] == 'B0']
            print(data)
            data = [i for i in G.nodes(data=True) if i[0] == 'product_1B']
            print(data)
            data = [i for i in G.nodes(data=True) if i[0] == 'B1']
            print(data)

            # idx = [i[0] for i in data]
            # data = [i[1] for i in data]
            # df = pd.DataFrame(data, index=idx)
            # print(df.to_markdown())
            # nodes = [n[1]['smiles'] for n in G.nodes(data=True) if n[0] == 'B1']
            # logger.error(f'SMILES: {nodes}')

            if errors == 'raise':
                exit()
            elif errors == 'ignore':
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
