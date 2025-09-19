import networkx as nx
from matplotlib import pyplot as plt


def get_reaction_graph(config: dict) -> nx.DiGraph:
    G = nx.DiGraph()
    for step in config['assembly']['steps']:
        name = step["name"]
        product = step["product"]
        for reactant in step["reactants"]:
            G.add_edge(reactant, product, label=name)
    return G


def visualize_reaction_graph(G: nx.DiGraph) -> plt.Axes:
    pos = nx.spring_layout(G, seed=42)  # layout for consistent placement
    fig, ax= plt.subplots(1, 1, figsize=(6, 4))
    nx.draw(G, pos, with_labels=True, node_size=2000, node_color="lightblue", arrows=True, ax=ax)
    return ax
