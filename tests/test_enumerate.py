from itertools import product
from pathlib import Path
import networkx as nx
from tqdm import tqdm
from loguru import logger
from delt_core.cli.library.api import get_reaction_graph, visualize_reaction_graph, complete_reaction_graph
from delt_core.utils import read_yaml

lib_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/tests')
config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-rechecked-enum/config.yaml')
cfg = read_yaml(config_path)

steps = cfg['library']['steps']
catalog = cfg['catalog']
reactions = catalog['reactions']
compounds = catalog['compounds']
building_blocks = {k: dict(smiles=None) for k in cfg['library']['building_blocks']}
products = {k: dict(smiles=None) for k in cfg['library']['products']}

# steps = [
#     ('A', 'R1'),
#     ('B', 'R1'),
#     ('R1', 'P1'),
#     ('P1', 'R2'),
#     ('C', 'R2'),
#     ('Precursor', 'R3'),
#     ('R3', 'C'),
# ]
#
# reactions = {'R1': dict(smirks='[A].[B]>>[P1]')}
# compounds = {i: dict(smiles=None) for i in ['A', 'B']}
# products = {'P1': dict(smiles=None)}

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
building_block_names = sorted(building_blocks)
lists = [cfg['whitelists'][bbn] for bbn in building_block_names]
combs = list(product(*lists))
library = []
# for comb in tqdm(combs[-1:]):
for i, comb in tqdm(enumerate(combs)):
    # break
    rnx = set(c['reaction'] for c in comb)
    prods = set([c['product'] for c in comb])

    reacts = set([c['reactant'] for c in comb]) - set(prods) - set(reactions)
    ancs = set(reacts)
    for r in reacts:
        ancs.update(nx.ancestors(G, r))

    rnx = rnx | ancs & set(reactions)
    prods = prods | ancs & set(products)
    comps = ancs & set(compounds)

    bbs = {bbn: bb
           for bbn, bb in zip(building_block_names, comb)
           if not pd.isna(bb['smiles'])}

    rnx = {r: cfg['catalog']['reactions'][r] for r in rnx}
    prods = {p: dict(smiles=None) for p in prods}
    comps = {c: cfg['catalog']['compounds'][c] for c in comps}

    nodes = {**comps, **prods, **bbs, **rnx}
    g = G.subgraph(nodes).copy()

    sinks = [n for n, d in g.out_degree() if d == 0]
    is_valid = len(sinks) == 1

    if (debug == 'all') or (debug == 'valid') and is_valid:
        ax = visualize_reaction_graph(g)
        ax.figure.savefig(lib_path / f'reaction_graph_{"_".join(str(c["index"]) for c in comb)}.png', dpi=300)
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
