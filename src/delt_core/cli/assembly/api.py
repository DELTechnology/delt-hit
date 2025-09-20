from pathlib import Path

import networkx as nx
from loguru import logger
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from delt_core.demultiplex.parser import catalog_from_excel
from delt_core.utils import read_yaml

def perform_reaction(smirks: str, reactants: list[str]) -> list[str]:

    rxn = rdChemReactions.ReactionFromSmarts(smirks)
    mols = [Chem.MolFromSmiles(r) for r in reactants]

    # the order of mols matters, determined by rxn!
    product_sets = []
    product_sets.extend(rxn.RunReactants(mols))
    product_sets.extend(rxn.RunReactants(mols[::-1]))

    prods = set()
    for tup in product_sets:
        for pmol in tup:
            # sanitize and canonicalize
            Chem.SanitizeMol(pmol, catchErrors=True)
            prods.add(Chem.MolToSmiles(pmol, canonical=True))
    return list(prods)[0]


class Assembly:

    def init(self, *, excel_path: Path):
        pass

    def reaction_graph(self, *, config_path: Path):
        from delt_core.assembly.reaction_graph import get_reaction_graph, visualize_reaction_graph

        config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/scribble/test-1/config.yaml')
        config = read_yaml(config_path)

        G = get_reaction_graph(config)

        ax = visualize_reaction_graph(G)
        save_dir = Path(config['experiment']['save_dir']).expanduser().resolve()
        name = Path(config['experiment']['name'])
        save_path = save_dir / name / 'reaction_graph.pdf'
        save_path.parent.mkdir(parents=True, exist_ok=True)
        ax.figure.tight_layout()
        ax.figure.savefig(save_path)
        logger.info(f"Saved reaction graph to {save_path}")

    def perform_reactions(self, *, config: dict, G: nx.DiGraph):
        from itertools import product
        UNIMOLECULAR = ['SR', 'DH', 'DT', 'SN2-1']

        bb_names = list(map(lambda x: x['name'], filter(lambda x: x['type'] == 'building_block', config['structure'])))
        bb_lists = [config['whitelists'][name] for name in bb_names]

        combinations = list(product(*bb_lists))
        for combo in combinations:
            break
            g = G.copy()

            attrs = {}
            nodes = set()
            for name, bb in zip(bb_names, combo):
                attrs[name] = bb
                reactant = bb['reactant']
                product = bb['product']
                reaction = bb['reaction']
                nodes = nodes | {name, reactant, product, reaction}
                if not reactant.startswith('product_'):
                    attrs[reactant] = scaffolds[reactant]

            g = g.subgraph(nodes)
            visualize_reaction_graph(g)
            nx.set_node_attributes(g, attrs)

            products = sorted([i for i in g.nodes if i.startswith('product_')])
            p = 'product_2'
            for p in products:
                # break
                # TODO: how to ensure order of reactants matches the smirk definition?
                reaction, = list(g.predecessors(p))

                if reaction == 'PASSTHROUGH':
                    p = g.predecessors(reaction)
                    product = g.nodes[p]['smiles']

                smirks = reactions[reaction]['smirks']
                reactants = [g.nodes[p]['smiles'] for p in sorted(g.predecessors(reaction))]
                product = perform_reaction(smirks, reactants)
                g.nodes[p]['smiles'] = product



    def prepare(self, *, config_path: Path):
        config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/scribble/test-1/config.yaml')
        config = read_yaml(config_path)


    def run(self, config_path: Path):
        pass