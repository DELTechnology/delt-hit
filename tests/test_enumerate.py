from itertools import product
from pathlib import Path
import networkx as nx
from tqdm import tqdm
from loguru import logger
from delt_core.cli.library.api import get_reaction_graph, visualize_reaction_graph, complete_reaction_graph
from delt_core.utils import read_yaml

lib_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/tests')
config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-rechecked-enum/config.yaml')
config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-rechecked-enum-v1/config.yaml')
cfg = read_yaml(config_path)

building_block_edges = cfg['library']['bb_edges']
other_edges = cfg['library']['other_edges']
steps = building_block_edges + other_edges
building_blocks = {k: dict(smiles=None) for k in cfg['library']['building_blocks']}
products = {k: dict(smiles=None) for k in cfg['library']['products']}

catalog = cfg['catalog']
reactions = catalog['reactions']
compounds = catalog['compounds']

G = get_reaction_graph(steps=building_block_edges,
                       reactions=reactions,
                       building_blocks=building_blocks,
                       compounds=compounds,
                       products=products)

ax = visualize_reaction_graph(G)
ax.figure.show()
ax.figure.savefig(lib_path.parent / 'building_block_reactions_graph.png', dpi=300)

add_G = get_reaction_graph(steps=other_edges,
                           reactions=reactions,
                           building_blocks=building_blocks,
                           compounds=compounds,
                           products=products)
ax = visualize_reaction_graph(add_G)
ax.figure.show()
ax.figure.savefig(lib_path.parent / 'additional_reactions_graph.png', dpi=300)

G = get_reaction_graph(steps=steps,
                       reactions=reactions,
                       building_blocks=building_blocks,
                       compounds=compounds,
                       products=products)

ax = visualize_reaction_graph(G)
ax.figure.show()
ax.figure.savefig(lib_path.parent / 'reaction_graph.png', dpi=300)



# %%
import pandas as pd
import matplotlib.pyplot as plt
debug = 'valid'
debug = 'all'
building_block_names = sorted(building_blocks)
lists = [cfg['whitelists'][bbn] for bbn in building_block_names]
combs = list(product(*lists))
library = []
# for comb in tqdm(combs[-1:]):
comb = combs[-1]
for i, comb in tqdm(enumerate(combs)):
    # break
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

    # additional_edges = list(filter(lambda step: set(step).intersection(bb_nodes), other_edges))
    additional_edges = [tuple(e) for e in additional_edges]
    additional_nodes = set([n for e in additional_edges for n in e])

    # NOTE: could be solved with nx.weakly_connected_components(add_G) as well
    # ancs = set()
    # for n in additional_nodes:
    #     ancs.update(nx.ancestors(add_G, n))
    # ancs_edges = list(filter(lambda e: ancs.intersection(e), other_edges))
    # edges = bb_edges + additional_edges + ancs_edges

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
        ax.figure.savefig(lib_path / f'reaction_graph_{i}_{"_".join(str(c["index"]) for c in comb)}.png', dpi=300)
        plt.close('all')
        ax.figure.show()

    if is_valid:
        assert len(sinks) == 1, f"Expected exactly one sink node, found {len(sinks)}"
        terminal = sinks[0]
    else:
        logger.warning(f'More than one sink node for combination: {i}')
        # NOTE: this means the combination is invalid,
        #   the participating reactions and compounds are not linearly connected. This happens for example if a
        #   certain product_1 is only used in a subset of subsequent reactions.
        continue

    nx.set_node_attributes(g, nodes)
    g = complete_reaction_graph(g)
    smiles = g.nodes[terminal]['smiles']
    record = {f'code_{i}': c['index'] for i, c in enumerate(comb)}
    record['smiles'] = smiles
    library.append(record)

# from delt_core.demultiplex.parser import library_and_catalog_from_excel
# from pathlib import Path
# path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/251021_NF2_library_rechecked_enumerate.xlsx')
# library_and_catalog_from_excel(path=path)
