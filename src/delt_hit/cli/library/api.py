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
        """Resolve the experiment output directory.

        Args:
            config_path: Path to the YAML config file.

        Returns:
            Experiment directory path.
        """
        cfg = read_yaml(config_path)
        exp_dir = Path(cfg['experiment']['save_dir']).expanduser().resolve() / cfg['experiment']['name']
        return exp_dir

    def get_library_path(self, *, config_path: Path) -> Path:
        """Resolve the default library parquet path.

        Args:
            config_path: Path to the YAML config file.

        Returns:
            Path to the library parquet file.
        """
        exp_dir = self.get_experiment_dir(config_path=config_path)
        lib_path = exp_dir / 'library' / 'library.parquet'
        return lib_path

    def get_named_library_path(self, *, config_path: Path, library_name: str) -> Path:
        """Resolve a named library parquet path inside the experiment library dir."""
        exp_dir = self.get_experiment_dir(config_path=config_path)
        return exp_dir / 'library' / f'{library_name}.parquet'

    def get_library_dir(self, *, config_path: Path) -> Path:
        """Resolve the experiment library output directory."""
        exp_dir = self.get_experiment_dir(config_path=config_path)
        return exp_dir / 'library'

    def enumerate(self, *, config_path: Path, debug: str = 'False', overwrite: bool = False,
                  graph_only: bool = False, errors: str = 'raise', building_block_ids: list[str] | None = None,
                  counts_path: Path | None = None, top_n: int | None = None, library_name: str | None = None):
        """Enumerate the combinatorial library from a config.

        Args:
            config_path: Path to the YAML config file.
            debug: Debug mode ('False', 'all', 'valid', 'invalid').
            overwrite: Whether to overwrite an existing library file.
            graph_only: Whether to stop after writing reaction graphs.
            errors: Error handling mode ('raise' or 'ignore').
            building_block_ids: Optional list of building block IDs to keep.
            counts_path: Optional path to a file with observed combinations.
            top_n: Optional cap on the number of input combinations to enumerate.
            library_name: Optional output parquet base name for filtered mode.
        """
        if counts_path is not None:
            assert library_name, "Filtered enumeration requires `library_name`"
            lib_path = self.get_named_library_path(config_path=config_path, library_name=library_name)
        else:
            lib_path = self.get_library_path(config_path=config_path)
        if lib_path.exists() and not overwrite:
            logger.info(f'Library {lib_path} exists')
            return

        if top_n is not None:
            assert top_n > 0, "`top_n` must be a positive integer"

        cfg = read_yaml(config_path)
        graph_bundle = prepare_graph_bundle(cfg=cfg)

        lib_path.parent.mkdir(parents=True, exist_ok=True)
        save_graph_visualizations(graph_bundle=graph_bundle, save_dir=lib_path.parent)

        logger.info(f'Saved reaction graph visualizations to {lib_path.parent}')

        if graph_only:
            return

        building_block_names = sorted(graph_bundle['building_blocks'])
        if building_block_ids:
            building_block_names = list(filter(lambda x: x in building_block_ids, building_block_names))

        df = enumerate_library_dataframe(
            cfg=cfg,
            graph_bundle=graph_bundle,
            building_block_names=building_block_names,
            counts_path=counts_path,
            top_n=top_n,
            debug=debug,
            errors=errors,
            save_dir=lib_path.parent,
        )
        df.to_parquet(lib_path, index=False)

    def visualize(self, *, config_path: Path, counts_path: Path | None = None, top_n: int = 12,
                  library_name: str | None = None, building_block_ids: list[str] | None = None,
                  output_name: str = "visualization", library_path: Path | None = None,
                  overwrite: bool = False, nrow: int = 10):
        """Generate reaction and molecule visualization panels.

        Args:
            config_path: Path to the YAML config file.
            counts_path: Optional demultiplex-style counts file for top-hit visualization.
            top_n: Number of rows/molecules to visualize from the filtered input or library.
            library_name: Optional named library to load or create in filtered mode.
            building_block_ids: Optional subset of building block IDs to consider.
            output_name: Base name for the generated PNG files.
            library_path: Optional explicit library parquet path to visualize.
            overwrite: Whether to overwrite an existing filtered library in visualization mode.
            nrow: Number of molecules per row in structure grids.
        """
        assert top_n > 0, "`top_n` must be a positive integer"
        cfg = read_yaml(config_path)
        lib_dir = self.get_library_dir(config_path=config_path)
        lib_dir.mkdir(parents=True, exist_ok=True)

        graph_bundle = prepare_graph_bundle(cfg=cfg)
        save_graph_visualizations(graph_bundle=graph_bundle, save_dir=lib_dir)

        reactions_dir = lib_dir / 'reactions'
        reactions_dir.mkdir(parents=True, exist_ok=True)
        visualize_reaction_schemes(cfg['catalog']['reactions'], save_dir=reactions_dir)

        building_block_names = sorted(graph_bundle['building_blocks'])
        if building_block_ids:
            building_block_names = list(filter(lambda x: x in building_block_ids, building_block_names))

        selected_library_path = library_path
        if selected_library_path is None:
            if counts_path is not None:
                assert library_name, "Filtered visualization requires `library_name`"
                selected_library_path = self.get_named_library_path(config_path=config_path, library_name=library_name)
                if overwrite or not selected_library_path.exists():
                    df = enumerate_library_dataframe(
                        cfg=cfg,
                        graph_bundle=graph_bundle,
                        building_block_names=building_block_names,
                        counts_path=counts_path,
                        top_n=top_n,
                        errors='raise',
                        save_dir=lib_dir,
                    )
                    df.to_parquet(selected_library_path, index=False)
            else:
                selected_library_path = self.get_library_path(config_path=config_path)

        assert selected_library_path.exists(), f"Library file not found at {selected_library_path}"
        library_df = pd.read_parquet(selected_library_path).head(top_n)

        product_ax = visualize_smiles(
            smiles=library_df['smiles'].dropna().tolist(),
            legends=[f"row_{idx}" for idx in library_df.index],
            title='Product Structures',
            nrow=nrow,
        )
        product_ax.figure.savefig(lib_dir / f"{output_name}_products.png", dpi=300)
        close_figure(product_ax.figure)

        for bb_idx, bb_name in enumerate(building_block_names):
            if f'code_{bb_idx}' not in library_df.columns:
                continue

            whitelist = cfg['whitelists'][bb_name]
            smiles = []
            legends = []
            for code in library_df[f'code_{bb_idx}']:
                if pd.isna(code):
                    continue
                entry = whitelist[int(code)]
                smiles.append(entry['smiles'])
                legends.append(f"{bb_name}:{int(code)}")

            if not smiles:
                continue

            ax = visualize_smiles(
                smiles=smiles,
                legends=legends,
                title=f'{bb_name} Building Blocks',
                nrow=nrow,
            )
            ax.figure.savefig(lib_dir / f"{output_name}_{bb_name}.png", dpi=300)
            close_figure(ax.figure)

    def properties(self, *, config_path: Path, library_name: str | None = None,
                   library_path: Path | None = None):
        """Compute molecular properties for a library and plot histograms.

        Args:
            config_path: Path to the YAML config file.
            library_name: Optional named library parquet to load from the experiment library dir.
            library_path: Optional library parquet override.
        """
        if library_path is not None:
            lib_path = library_path
            output_name = library_path.stem
        elif library_name is not None:
            lib_path = self.get_named_library_path(config_path=config_path, library_name=library_name)
            output_name = library_name
        else:
            lib_path = self.get_library_path(config_path=config_path)
            output_name = 'properties'

        assert lib_path.exists(), f"Library file not found at {lib_path}"

        save_dir = lib_path.parent / 'properties'
        save_dir.mkdir(parents=True, exist_ok=True)

        df = pd.read_parquet(lib_path)
        df = self.compute_properties(data=df)
        df.to_parquet(save_dir / f'{output_name}.parquet', index=False)

        prop_names = [col for col in df.columns if col.startswith('prop_')]
        plt.close('all')
        for name in prop_names:
            ax = self.plot_property(data=df, name=name)
            ax.figure.savefig(save_dir / f"{name}.png")
            plt.close(ax.figure)

    def compute_properties(self, data: pd.DataFrame) -> pd.DataFrame:
        """Compute RDKit property columns for each SMILES entry.

        Args:
            data: DataFrame with a ``smiles`` column.

        Returns:
            DataFrame with appended ``prop_*`` columns.
        """

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
        """Plot a histogram for a single property column.

        Args:
            data: DataFrame with property columns.
            name: Property column name to plot.

        Returns:
            Matplotlib Axes with the plot.

        Raises:
            ValueError: If the property column is missing.
        """
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
        """Generate molecular representations for the library.

        Args:
            config_path: Path to the YAML config file.
            method: Representation type ('morgan' or 'bert').
            library_path: Optional library parquet override.
        """
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
    """Compute BERT representations for a SMILES library.

    Args:
        model_name: Name of the model to use.
        path: Path to a parquet file with a ``smiles`` column.
        save_path: Path to write the numpy array.
        device: Torch device name.

    Raises:
        ValueError: If the model name is unknown.
    """
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
    """Generate pooled BERT embeddings for SMILES strings.

    Args:
        smiles: List of SMILES strings.
        device: Torch device name.

    Returns:
        A list of numpy arrays with pooled embeddings.
    """
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
    """Compute Morgan fingerprints and save them.

    Args:
        smiles: List of SMILES strings.
        save_path: Path to write the sparse matrix.
    """
    fps = []
    for smiles in tqdm(smiles):
        fp = get_morgan_fp(smiles)
        fps.append(fp)

    fps = sparse.vstack(fps, format="csr")
    save_path.parent.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(save_path, fps)

    logger.info(f"Fingerprints saved to {save_path}")


def get_morgan_fp(smiles, radius=2, n_bits=2048) -> sparse.csr_array:
    """Create a Morgan fingerprint for a SMILES string.

    Args:
        smiles: SMILES string to featurize.
        radius: Morgan radius.
        n_bits: Number of fingerprint bits.

    Returns:
        A sparse CSR fingerprint vector.
    """
    mol = Chem.MolFromSmiles(smiles)
    mfpgen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = mfpgen.GetFingerprint(mol)
    fp = sparse.csr_array(fp, dtype=np.uint8)
    return fp


def get_dummy_library() -> pd.DataFrame:
    """Return a small demo library of SMILES strings.

    Returns:
        DataFrame with demo ``smiles`` values.
    """
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


def get_combinations(*,
                     cfg: dict,
                     building_block_names: list[str],
                     counts_path: Path | None = None,
                     top_n: int | None = None) -> list[tuple[dict, ...]]:
    """Return enumeration combinations from either whitelist product or filtered input."""
    if counts_path is None:
        lists = [cfg['whitelists'][bbn] for bbn in building_block_names]
        return list(product(*lists))

    return load_filtered_combinations(
        cfg=cfg,
        building_block_names=building_block_names,
        counts_path=counts_path,
        top_n=top_n,
    )


def load_filtered_combinations(*,
                               cfg: dict,
                               building_block_names: list[str],
                               counts_path: Path,
                               top_n: int | None = None) -> list[tuple[dict, ...]]:
    """Load explicit building-block combinations from a demultiplex-style counts file."""
    counts_path = Path(counts_path).expanduser().resolve()
    assert counts_path.exists(), f"Combination file not found at {counts_path}"

    df = pd.read_csv(counts_path, sep=None, engine='python')

    code_columns = [f'code_{i}' for i in range(len(building_block_names))]
    assert all(name in df.columns for name in code_columns), (
        "Filtered enumeration expects the same `code_*` columns produced by demultiplex counts files. "
        f"Missing one or more required columns from {code_columns}."
    )

    if top_n is not None:
        df = df.head(top_n)

    combinations = []
    for row_idx, row in df.loc[:, code_columns].iterrows():
        combination = []
        for bb_idx, bb_name in enumerate(building_block_names):
            whitelist = cfg['whitelists'][bb_name]
            raw_value = row[f'code_{bb_idx}']

            assert not pd.isna(raw_value), f"Missing value for {bb_name} at row {row_idx}"

            try:
                value = int(raw_value)
            except (TypeError, ValueError) as exc:
                raise AssertionError(f"Invalid value for {bb_name} at row {row_idx}: {raw_value}") from exc

            assert 0 <= value < len(whitelist), (
                f"Value {value} for {bb_name} at row {row_idx} is outside whitelist range [0, {len(whitelist) - 1}]"
            )
            combination.append(whitelist[value])
        combinations.append(tuple(combination))

    return combinations


def get_reaction_graph(steps: list,
                       reactions: dict,
                       compounds: dict,
                       products: dict,
                       building_blocks: dict = None) -> nx.DiGraph:
    """Build a directed reaction graph from inputs.

    Args:
        steps: Edge list of the reaction graph.
        reactions: Reaction metadata keyed by name.
        compounds: Compound metadata keyed by name.
        products: Product metadata keyed by name.
        building_blocks: Optional building block metadata.

    Returns:
        A populated ``networkx`` directed graph.
    """
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


def prepare_graph_bundle(cfg: dict) -> dict:
    """Prepare graph inputs shared by enumeration and visualization."""
    building_block_edges = cfg['library']['bb_edges']
    other_edges = cfg['library']['other_edges']
    steps = building_block_edges + other_edges
    building_blocks = {k: dict(smiles=None) for k in cfg['library']['building_blocks']}
    products = {k: dict(smiles=None) for k in cfg['library']['products']}
    reactions = cfg['catalog']['reactions']
    compounds = cfg['catalog']['compounds']

    bb_G = get_reaction_graph(
        steps=building_block_edges,
        reactions=reactions,
        building_blocks=building_blocks,
        compounds=compounds,
        products=products,
    )
    add_G = get_reaction_graph(
        steps=other_edges,
        reactions=reactions,
        building_blocks=building_blocks,
        compounds=compounds,
        products=products,
    )
    G = get_reaction_graph(
        steps=steps,
        reactions=reactions,
        building_blocks=building_blocks,
        compounds=compounds,
        products=products,
    )

    return {
        'bb_G': bb_G,
        'add_G': add_G,
        'G': G,
        'building_blocks': building_blocks,
        'products': products,
        'reactions': reactions,
        'compounds': compounds,
    }


def save_graph_visualizations(*, graph_bundle: dict, save_dir: Path) -> None:
    """Write the standard reaction graph PNGs to disk."""
    for graph, filename in [
        (graph_bundle['bb_G'], 'building_block_reactions_graph.png'),
        (graph_bundle['add_G'], 'additional_reactions_graph.png'),
        (graph_bundle['G'], 'reaction_graph.png'),
    ]:
        ax = visualize_reaction_graph(graph)
        ax.figure.savefig(save_dir / filename, dpi=300)
        close_figure(ax.figure)


def enumerate_library_dataframe(*,
                                cfg: dict,
                                graph_bundle: dict,
                                building_block_names: list[str],
                                counts_path: Path | None = None,
                                top_n: int | None = None,
                                debug: str = 'False',
                                errors: str = 'raise',
                                save_dir: Path | None = None) -> pd.DataFrame:
    """Enumerate a library and return it as a dataframe."""
    reactions = graph_bundle['reactions']
    products = graph_bundle['products']
    compounds = graph_bundle['compounds']
    add_G = graph_bundle['add_G']
    G = graph_bundle['G']

    combs = get_combinations(
        cfg=cfg,
        building_block_names=building_block_names,
        counts_path=counts_path,
        top_n=top_n,
    )

    library = []
    logger.info('Starting enumeration of library...')
    for i, comb in tqdm(enumerate(combs)):
        bb_edges = [(bb, c['reaction']) for bb, c in zip(building_block_names, comb)]
        bb_edges += [(c['reaction'], c['product']) for c in comb]
        bb_edges += [(c['educt'], c['reaction']) for c in comb]
        bb_nodes = {n for e in bb_edges for n in e}

        additional_edges = set()
        subgraphs = list(nx.weakly_connected_components(add_G))
        for n in bb_nodes:
            for sgn in subgraphs:
                sg = add_G.subgraph(sgn).copy()
                if n in sg and sg.out_degree(n) == 0:
                    additional_edges.update(sg.edges)

        edges = bb_edges + [tuple(e) for e in additional_edges]
        edges = [tuple(e) for e in edges]
        nodes = [n for e in edges for n in e]

        rnx = {r: cfg['catalog']['reactions'][r] for r in nodes if r in reactions}
        prods = {p: dict(smiles=None) for p in nodes if p in products}
        comps = {c: cfg['catalog']['compounds'][c] for c in nodes if c in compounds}
        bbs = {bbn: bb for bbn, bb in zip(building_block_names, comb) if not pd.isna(bb['smiles'])}

        node_attrs = {**comps, **prods, **bbs, **rnx}

        g = G.edge_subgraph(edges=edges).copy()
        sinks = [n for n, d in g.out_degree() if d == 0]
        is_valid = len(sinks) == 1

        if ((debug == 'all') or (debug == 'valid')) and is_valid and save_dir is not None:
            ax = visualize_reaction_graph(g)
            ax.figure.savefig(
                save_dir / f'reaction_graph_combination={i}_{"_".join(str(c["index"]) for c in comb)}.png',
                dpi=300,
            )
            plt.close(ax.figure)

        if is_valid:
            terminal = sinks[0]
        else:
            logger.warning(
                f'More than one terminal node for combination: {i}',
                f"Run with debug='all' to visualize reaction graphs with multiple terminal nodes.",
            )
            continue

        nx.set_node_attributes(g, node_attrs)

        try:
            g = complete_reaction_graph(g, errors=errors)
        except Exception as e:
            if debug == 'invalid' and save_dir is not None:
                ax = visualize_reaction_graph(g)
                ax.figure.savefig(
                    save_dir / f'reaction_graph_combination={i}_{"_".join(str(c["index"]) for c in comb)}.png',
                    dpi=300,
                )
                plt.close(ax.figure)

            if errors == 'raise':
                raise e
            if errors == 'ignore':
                continue

        smiles = g.nodes[terminal]['smiles']
        record = {f'code_{j}': c['index'] for j, c in enumerate(comb)}
        record['smiles'] = smiles
        library.append(record)

    df = pd.DataFrame(library)
    if 'smiles' in df.columns:
        df = df[df.smiles.notna()]
    return df


def visualize_reaction_graph(G: nx.DiGraph) -> plt.Axes:
    """Render a reaction graph with typed node coloring.

    Args:
        G: Reaction graph.

    Returns:
        Matplotlib Axes with the graph visualization.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    compounds = [n for n, d in G.nodes(data=True) if d.get("type") == "compound"]
    reactions = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    products = [n for n, d in G.nodes(data=True) if d.get("type") == "product"]
    building_blocks = [n for n, d in G.nodes(data=True) if d.get("type") == "building_block"]

    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")

    nx.draw_networkx_nodes(G, pos, nodelist=compounds, node_color="lightblue", node_shape="o", node_size=500, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=building_blocks, node_color="mediumorchid", node_shape="o", node_size=500, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=products, node_color="salmon", node_shape="o", node_size=500, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=reactions, node_color="lightgreen", node_shape="s", node_size=600, ax=ax)

    labels = {node: node for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax, font_size=8)
    nx.draw_networkx_edges(G, pos, ax=ax, arrows=True)
    ax.set_axis_off()
    ax.figure.tight_layout()

    return ax


def visualize_reaction_schemes(reactions: dict, save_dir: Path) -> None:
    """Save one PNG per reaction directly from config reaction data.

    Args:
        reactions: Mapping of reaction name to reaction metadata (must include 'smirks').
        save_dir: Directory to write PNG files into.
    """
    for name, data in reactions.items():
        smirks = data.get("smirks")
        fig, ax = plt.subplots(1, 1, figsize=(10, 3))
        ax.set_axis_off()

        if not smirks or (isinstance(smirks, float) and pd.isna(smirks)):
            ax.text(0.5, 0.5, name, ha='center', va='center', fontsize=14)
        else:
            rxn = rdChemReactions.ReactionFromSmarts(smirks)
            img = Draw.ReactionToImage(rxn, subImgSize=(300, 150))
            ax.imshow(img)
            ax.set_title(name, fontsize=12)

        fig.tight_layout()
        fig.savefig(save_dir / f"{name}.png", dpi=300)
        close_figure(fig)


def close_figure(fig) -> None:
    """Close a matplotlib figure when possible."""
    try:
        plt.close(fig)
    except TypeError:
        return


def find_next_reaction(G: nx.DiGraph):
    """Find the next reaction with all inputs defined.

    Args:
        G: Reaction graph with ``smiles`` node attributes.

    Returns:
        A dict describing the next reaction step, or None.
    """
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    for node in reaction_nodes:
        preds = sorted(G.predecessors(node))
        succ, = sorted(G.successors(node))  # note: currently only one product per reaction

        if G.nodes[succ]['smiles'] is None and all([G.nodes[i]['smiles'] is not None for i in preds]):
            return {'reactants': preds, 'reaction': node, 'product': succ}

    return None


def perform_reaction(smirks: str, reactants: list[str], use_smiles: bool = False) -> list[str]:
    """Run an RDKit reaction and return unique products.

    Args:
        smirks: Reaction SMARTS/SMIRKS string.
        reactants: List of reactant SMILES.
        use_smiles: Whether to treat the reaction as SMILES.

    Returns:
        A sorted list of product SMILES.
    """
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
    """Iteratively fill in missing product SMILES in a graph.

    Args:
        G: Reaction graph with ``smiles`` node attributes.
        errors: Error handling mode ('raise' or 'ignore').

    Returns:
        The updated reaction graph.
    """
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
            logger.error(f'Reaction graph at error: {G.nodes(data=True)}\n')

            if errors == 'raise':
                exit()
            elif errors == 'ignore':
                break

    return G


def visualize_smiles(smiles: list[str], nrow: int = 25, legends: list[str] | None = None,
                     title: str = 'Product Structures'):
    """Create a grid image of molecules.

    Args:
        smiles: List of SMILES strings.
        nrow: Maximum number of molecules per row.
        legends: Optional per-molecule legends.
        title: Figure title.

    Returns:
        Matplotlib Axes containing the image.
    """
    if len(smiles) == 0:
        fig, ax = plt.subplots(1, 1, figsize=(8, 3))
        ax.set_axis_off()
        ax.text(0.5, 0.5, "No structures to display", ha='center', va='center')
        ax.set_title(title)
        fig.tight_layout()
        return ax

    mols = [Chem.MolFromSmiles(s) for s in smiles]
    legends = legends or None
    nrow = min(nrow, len(mols))
    img = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=nrow, subImgSize=(200, 200))
    plt.figure(figsize=(10, 6))
    ax = plt.imshow(img)
    ax.axes.set_axis_off()
    ax.axes.set_title(title)
    ax.figure.tight_layout()
    return ax
