from pathlib import Path
import networkx as nx
from delt_core.demultiplex.parser import catalog_from_excel
from delt_core.utils import read_yaml

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from pathlib import Path

import networkx as nx
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from delt_core.demultiplex.parser import catalog_from_excel
from delt_core.utils import read_yaml


class Assembly:

    def init(self, *, excel_path: Path):
        pass

    def reaction_graph(self, *, config_path: Path):
        from delt_core.assembly.reaction_graph import get_reaction_graph, visualize_reaction_graph

        # config_path = Path('/Users/adrianomartinelli/projects/delt/delt-core/scribble/test-1/config.yaml')
        config = read_yaml(config_path)

        G = get_reaction_graph(config)

        ax = visualize_reaction_graph(G)
        save_dir = Path(config['experiment']['save_dir'])
        name = Path(config['experiment']['name'])
        save_path = save_dir / name / 'reaction_graph.pdf'
        save_path.parent.mkdir(parents=True, exist_ok=True)
        ax.figure.tight_layout()
        ax.figure.savefig(save_path)

    def perform_reactions(self, *, config: dict, G: nx.DiGraph):
        from itertools import product
        # NOTE: this should not be needed anymore
        UNIMOLECULAR = ['SR', 'DH', 'DT', 'SN2-1']

        bb_names = list(map(lambda x: x['name'], filter(lambda x: x['type'] == 'building_block', config['structure'])))
        bb_lists = [config['whitelists'][name] for name in bb_names]

        catalog = config['catalog']
        scaffolds = catalog['scaffolds']
        reactions = catalog['reactions']

        combinations = list(product(*bb_lists))
        for combo in combinations:
            g = G.copy()

            attrs = {}
            if 'scaffold' in g.nodes:
                name = combo[0]['scaffold']
                attrs['scaffold'] = scaffolds[name]

            for name, bb in zip(bb_names, combo):
                attrs[name] = bb

            nx.set_node_attributes(g, attrs)

            products = sorted([i for i in g.nodes if i.startswith('P')])
            # p = 'P2'
            for p in products:
                # TODO: how to ensure order of reactants matches the smirk definition?
                predecessors = sorted(G.predecessors(p))

                reactants = [g.nodes[p]['smiles'] for p in predecessors]
                mols = [Chem.MolFromSmiles(r) for r in reactants]

                # NOTE: the building block carries the reaction info
                # reaction = {g.nodes[p].get('reaction') for p in predecessors if g.nodes[p].get('reaction') is not None}
                reaction = {g.nodes[p].get('reaction') for p in predecessors if p.startswith('B')}
                assert len(reaction) == 1
                reaction = reaction.pop()
                smirks = reactions[reaction]['smirks']
                rxn = rdChemReactions.ReactionFromSmarts(smirks)

                # the order of mols matters, determined by rxn!
                product_sets = []
                product_sets.extend(rxn.RunReactants(mols))
                # TODO: remove after fixing the smirks
                product_sets.extend(rxn.RunReactants(mols[::-1]))

                prods = set()
                for tup in product_sets:
                    for pmol in tup:
                        # sanitize and canonicalize
                        Chem.SanitizeMol(pmol, catchErrors=False)
                        prods.add(Chem.MolToSmiles(pmol, canonical=True))

                assert len(prods) == 1, f"Expected one product, got {prods} for reaction {reaction} ({smirks}) with reactants {reactants}"

                g.nodes[p]['smiles'] = prods.pop()



    def prepare(self, *, config_path: Path):
        config_path = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/test-1/config.yaml')
        config = read_yaml(config_path)

        save_dir = Path(config['experiment']['save_dir'])
        name = config['experiment']['name']

        save_path = save_dir / name / 'library'
        path = Path('/Users/adrianomartinelli/projects/delt/delt-core/templates/library.xlsx')
        catalog = catalog_from_excel(path)

        catalog = config['catalog']
        scaffolds = catalog['scaffolds']
        building_blocks = list(filter(lambda x: x['type'] == 'building_block', config['structure']))
        reactions = catalog['reactions']




    def run(self, config_path: Path):
        pass