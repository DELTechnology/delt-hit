import networkx as nx
from matplotlib import pyplot as plt
import yaml

config = yaml.safe_load(open("/Users/adrianomartinelli/projects/delt/delt-core/scribble/test-1/config.yaml"))


def get_reaction_graph(config: dict) -> nx.DiGraph:
    G = nx.DiGraph()

    reactions = set()
    compounds = {'B0', 'B1'}
    for step in ['B0', 'B1']:
        reactions.update(i['reaction'] for i in config['whitelists'][step])
        compounds.update(i['reactant'] for i in config['whitelists'][step])
        compounds.update(i['product'] for i in config['whitelists'][step])

        edges = set((step, i['reaction']) for i in config['whitelists'][step])
        G.add_edges_from(edges)

        edges = set((i['reactant'], i['reaction']) for i in config['whitelists'][step])
        G.add_edges_from(edges)

        edges = set((i['reaction'], i['product']) for i in config['whitelists'][step])
        G.add_edges_from(edges)

    attrs = {i: {'type': 'reaction'} for i in reactions}
    nx.set_node_attributes(G, attrs)
    attrs = {i: {'type': 'compound'} for i in compounds}
    nx.set_node_attributes(G, attrs)

    return G


def visualize_reaction_graph(G: nx.DiGraph) -> plt.Axes:
    fig, ax= plt.subplots(1, 1, figsize=(12, 8))

    compounds = [n for n, d in G.nodes(data=True) if d.get("type") == "compound"]
    reactions = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]

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
