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

from delt_hit.cli.helper import (
    get_experiment_dir,
    get_library_dir,
    get_library_path,
    get_named_library_path,
)
from delt_hit.utils import read_yaml


class Library:
    def enumerate(self, *, config_path: Path, debug: str = 'False', overwrite: bool = False,
                  errors: str = 'raise', counts_path: Path | None = None,
                  top_n: int | None = None, library_name: str | None = None):
        """Enumerate the combinatorial library from a config.

        Args:
            config_path: Path to the YAML config file.
            debug: Debug output mode. Use 'all' or 'valid' to save
                per-combination reaction graphs for combinations with a
                single terminal node, or 'invalid' to save graphs for
                combinations that fail during reaction execution.
            overwrite: Whether to overwrite an existing library file.
            errors: Error handling mode. Use 'raise' to stop on the first
                reaction error or 'ignore' to skip failing combinations and
                continue enumeration.
            counts_path: Optional path to a file with observed combinations.
            top_n: Optional cap on the number of input combinations to enumerate.
            library_name: Optional output parquet base name for filtered mode.
        """
        if counts_path is not None:
            assert library_name, "Filtered enumeration requires `library_name`"
            lib_path = get_named_library_path(config_path, library_name)
        else:
            lib_path = get_library_path(config_path)
        if lib_path.exists() and not overwrite:
            logger.info(f'Library {lib_path} exists')
            return

        if top_n is not None:
            assert top_n > 0, "`top_n` must be a positive integer"

        cfg = read_yaml(config_path)
        graph_bundle = prepare_graph_bundle(cfg=cfg)

        lib_path.parent.mkdir(parents=True, exist_ok=True)

        building_block_names = sorted(graph_bundle['building_blocks'])

        if is_dual_display_library(cfg, building_block_names):
            df = enumerate_dual_display(
                cfg=cfg,
                building_block_names=building_block_names,
                counts_path=counts_path,
                top_n=top_n,
                debug=debug,
                errors=errors,
                save_dir=lib_path.parent,
            )
        else:
            df = enumerate_single_strand(
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
        elif library_name is not None:
            lib_path = get_named_library_path(config_path, library_name)
        else:
            lib_path = get_library_path(config_path)

        assert lib_path.exists(), f"Library file not found at {lib_path}"

        save_dir = lib_path.parent / 'properties' / lib_path.stem
        save_dir.mkdir(parents=True, exist_ok=True)

        df = pd.read_parquet(lib_path)
        assert 'smiles' in df.columns, "Dual-display libraries with `smiles_a`/`smiles_b` are not supported by `library properties`"
        df = self.compute_properties(data=df)
        df.to_parquet(save_dir / 'properties.parquet', index=False)

        prop_names = [col for col in df.columns if col.startswith('prop_')]
        plt.close('all')
        for name in tqdm(prop_names, desc="Property histograms"):
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
        for smiles in tqdm(data['smiles'], desc="Molecular properties"):
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
        exp_dir = get_experiment_dir(config_path)

        save_dir = exp_dir / 'representations'
        save_dir.mkdir(parents=True, exist_ok=True)

        lib_path = library_path or get_library_path(config_path)
        df = pd.read_parquet(lib_path)
        assert 'smiles' in df.columns, "Dual-display libraries with `smiles_a`/`smiles_b` are not supported by `library represent`"
        smiles = df.smiles

        match method:
            case 'morgan':
                run_morgan(smiles, save_path=save_dir / 'morgan.npz')
            case 'bert':
                run_morgan(smiles, save_path=save_dir / 'bert.npz')

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


def get_code_column_name(building_block_name: str) -> str:
    """Return the canonical code column for a building-block family."""
    assert building_block_name.startswith('B'), f"Building block name must start with 'B', got {building_block_name}"
    return f"code_{building_block_name[1:]}"


def load_filtered_combinations(*,
                               cfg: dict,
                               building_block_names: list[str],
                               counts_path: Path,
                               top_n: int | None = None) -> list[tuple[dict, ...]]:
    """Load explicit building-block combinations from a demultiplex-style counts file."""
    counts_path = Path(counts_path).expanduser().resolve()
    assert counts_path.exists(), f"Combination file not found at {counts_path}"

    df = pd.read_csv(counts_path, sep=None, engine='python')

    code_columns = [get_code_column_name(name) for name in building_block_names]
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
            raw_value = row[get_code_column_name(bb_name)]

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


def is_dual_display_library(cfg: dict, building_block_names: list[str]) -> bool:
    """Return whether the config declares an explicit dual-display library."""
    strand_by_building_block = get_building_block_strands(cfg)
    return bool(strand_by_building_block) and all(name in strand_by_building_block for name in building_block_names)


def get_building_block_strands(cfg: dict) -> dict[str, str]:
    """Extract explicit strand assignments for building-block regions."""
    structure = cfg['structure']
    return {
        entry['name']: entry['strand']
        for entry in structure
        if entry.get('type') == 'building_block' and entry.get('strand') is not None
    }


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


def save_figure_outputs(fig, output_path: Path, *, dpi: int = 300) -> None:
    """Save a figure as PNG with compact margins."""
    tight_layout = getattr(fig, "tight_layout", None)
    if callable(tight_layout):
        tight_layout()
    fig.savefig(output_path.with_suffix(".png"), dpi=dpi, bbox_inches="tight")



def enumerate_single_strand(*,
                            cfg: dict,
                            graph_bundle: dict,
                            building_block_names: list[str],
                            combinations: list[tuple[dict, ...]] | None = None,
                            counts_path: Path | None = None,
                            top_n: int | None = None,
                            debug: str = 'False',
                            errors: str = 'raise',
                            save_dir: Path | None = None,
                            label: str | None = None) -> pd.DataFrame:
    """Enumerate a single-display library and return it as a dataframe."""
    reactions = graph_bundle['reactions']
    products = graph_bundle['products']
    compounds = graph_bundle['compounds']
    add_G = graph_bundle['add_G']
    G = graph_bundle['G']
    if combinations is None:
        combs = get_combinations(
            cfg=cfg,
            building_block_names=building_block_names,
            counts_path=counts_path,
            top_n=top_n,
        )
    else:
        combs = combinations

    library = []
    suffix = f' ({label})' if label else ''
    logger.info(f'Starting enumeration of library{suffix}...')
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
                if n in sg and (sg.out_degree(n) == 0 or sg.in_degree(n) == 0):
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

        record = {get_code_column_name(bb_name): c['index'] for bb_name, c in zip(building_block_names, comb)}
        record['smiles'] = g.nodes[terminal]['smiles']
        library.append(record)

    df = pd.DataFrame(library)
    if 'smiles' in df.columns:
        df = df[df.smiles.notna()]
    return df


def enumerate_dual_display(*,
                           cfg: dict,
                           building_block_names: list[str],
                           counts_path: Path | None = None,
                           top_n: int | None = None,
                           debug: str = 'False',
                           errors: str = 'raise',
                           save_dir: Path | None = None) -> pd.DataFrame:
    """Enumerate both strands independently, then pair their products."""
    strand_by_building_block = get_building_block_strands(cfg)
    strand_a_names = [name for name in building_block_names if strand_by_building_block[name] == 'A']
    strand_b_names = [name for name in building_block_names if strand_by_building_block[name] == 'B']
    strand_a_cfg = get_strand_cfg(cfg=cfg, building_block_names=strand_a_names)
    strand_b_cfg = get_strand_cfg(cfg=cfg, building_block_names=strand_b_names)
    strand_a_graph_bundle = prepare_graph_bundle(cfg=strand_a_cfg)
    strand_b_graph_bundle = prepare_graph_bundle(cfg=strand_b_cfg)

    strand_a_counts_path = None
    strand_b_counts_path = None
    if counts_path is not None:
        dual_display_combinations = load_filtered_combinations(
            cfg=cfg,
            building_block_names=building_block_names,
            counts_path=counts_path,
            top_n=top_n,
        )
        strand_a_combinations = project_strand_combinations(
            combinations=dual_display_combinations,
            building_block_names=building_block_names,
            strand_building_block_names=strand_a_names,
        )
        strand_b_combinations = project_strand_combinations(
            combinations=dual_display_combinations,
            building_block_names=building_block_names,
            strand_building_block_names=strand_b_names,
        )
    else:
        strand_a_combinations = None
        strand_b_combinations = None

    strand_a_df = enumerate_single_strand(
        cfg=strand_a_cfg,
        graph_bundle=strand_a_graph_bundle,
        building_block_names=strand_a_names,
        combinations=strand_a_combinations,
        debug=debug,
        errors=errors,
        save_dir=save_dir,
        label='strand A',
    )
    strand_b_df = enumerate_single_strand(
        cfg=strand_b_cfg,
        graph_bundle=strand_b_graph_bundle,
        building_block_names=strand_b_names,
        combinations=strand_b_combinations,
        debug=debug,
        errors=errors,
        save_dir=save_dir,
        label='strand B',
    )

    strand_a_df = strand_a_df.rename(columns={'smiles': 'smiles_a'})
    strand_b_df = strand_b_df.rename(columns={'smiles': 'smiles_b'})

    strand_a_records = strand_a_df.to_dict('records')
    strand_b_records = strand_b_df.to_dict('records')

    combined_records = []
    for strand_a_record, strand_b_record in product(strand_a_records, strand_b_records):
        record = {**strand_a_record, **strand_b_record}
        combined_records.append(record)

    code_columns = [get_code_column_name(name) for name in building_block_names]
    if not combined_records:
        return pd.DataFrame(columns=[*code_columns, 'smiles_a', 'smiles_b'])
    df = pd.DataFrame(combined_records)[[*code_columns, 'smiles_a', 'smiles_b']]
    df = df.sort_values(code_columns).reset_index(drop=True)
    return df


def get_strand_cfg(*, cfg: dict, building_block_names: list[str]) -> dict:
    """Return a strand-specific config with only the requested building blocks."""
    strand_educts = {
        entry['educt']
        for name in building_block_names
        for entry in cfg['whitelists'][name]
    }
    strand_products = {
        entry['product']
        for name in building_block_names
        for entry in cfg['whitelists'][name]
    }
    strand_reactions = {
        entry['reaction']
        for name in building_block_names
        for entry in cfg['whitelists'][name]
    }
    strand_compounds = set(cfg['catalog']['compounds']) & strand_educts
    strand_nodes = set(building_block_names) | strand_products | strand_reactions | strand_compounds
    strand_library = {
        **cfg['library'],
        'bb_edges': [edge for edge in cfg['library']['bb_edges'] if edge[0] in strand_nodes and edge[1] in strand_nodes],
        'other_edges': [
            edge for edge in cfg['library']['other_edges']
            if edge[0] in strand_nodes and edge[1] in strand_nodes
        ],
        'building_blocks': building_block_names,
        'products': sorted(strand_products),
        'reactions': sorted(strand_reactions),
    }
    strand_catalog = {
        **cfg['catalog'],
        'reactions': {
            name: data for name, data in cfg['catalog']['reactions'].items()
            if name in strand_reactions
        },
    }
    strand_structure = [entry for entry in cfg['structure'] if entry['name'] in building_block_names]
    strand_whitelists = {name: cfg['whitelists'][name] for name in building_block_names}

    return {
        **cfg,
        'library': strand_library,
        'catalog': strand_catalog,
        'structure': strand_structure,
        'whitelists': strand_whitelists,
    }


def project_strand_combinations(*,
                                combinations: list[tuple[dict, ...]],
                                building_block_names: list[str],
                                strand_building_block_names: list[str]) -> list[tuple[dict, ...]]:
    """Project dual-display combinations onto one strand and drop duplicates."""
    strand_indices = [building_block_names.index(name) for name in strand_building_block_names]
    unique = {}
    for combination in combinations:
        strand_combination = tuple(combination[idx] for idx in strand_indices)
        unique[tuple(entry['index'] for entry in strand_combination)] = strand_combination
    return sorted(unique.values(), key=lambda combination: tuple(entry['index'] for entry in combination))


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


def visualize_reaction_schemes(reactions: dict, save_dir: Path, *, dpi: int = 300) -> None:
    """Save one reaction scheme per reaction directly from config reaction data.

    Args:
        reactions: Mapping of reaction name to reaction metadata (must include 'smirks').
        save_dir: Directory to write PNG files into.
    """
    for name, data in reactions.items():
        smirks = data.get("smirks")
        fig, ax = plt.subplots(1, 1, figsize=(10, 3))
        ax.set_axis_off()

        if not smirks or (isinstance(smirks, float) and pd.isna(smirks)):
            pass
        else:
            rxn = rdChemReactions.ReactionFromSmarts(smirks)
            img = Draw.ReactionToImage(rxn, subImgSize=(300, 150))
            ax.imshow(img)

        save_figure_outputs(fig, save_dir / name, dpi=dpi)
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

        if G.nodes[succ].get('smiles') is None and all([G.nodes[i]['smiles'] is not None for i in preds]):
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
            products.add(Chem.MolToSmiles(pmol, canonical=True, kekuleSmiles=False, isomericSmiles=True))

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
        next_reaction = None
        products = None

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


def visualize_smiles(
    smiles: list[str],
    nrow: int = 25,
    legends: list[str] | None = None,
    title: str = 'Structures',
    *,
    sub_img_size: tuple[int, int] = (300, 300),
):
    """Create a grid image of molecules.

    Args:
        smiles: List of SMILES strings.
        nrow: Maximum number of molecules per row.
        legends: Optional per-molecule legends.
        title: Figure title.
        sub_img_size: Pixel size used by RDKit for each molecule tile.

    Returns:
        Matplotlib Axes containing the image.
    """
    if len(smiles) == 0:
        fig, ax = plt.subplots(1, 1, figsize=(8, 3))
        ax.set_axis_off()
        ax.set_title(title)
        return ax

    mols = [Chem.MolFromSmiles(s) for s in smiles]
    legends = legends or None
    nrow = min(nrow, len(mols))
    img = Draw.MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=nrow,
        subImgSize=sub_img_size,
    )
    width_px, height_px = img.size
    fig, ax = plt.subplots(
        1,
        1,
        figsize=(width_px / 100, height_px / 100),
        dpi=100,
    )
    ax.imshow(img, interpolation="none")
    ax.axes.set_axis_off()
    ax.axes.set_title(title)
    return ax
