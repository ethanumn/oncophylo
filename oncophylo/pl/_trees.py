import networkx as nx 
import pandas as pd 
from IPython.display import Image, display
from networkx.drawing.nx_pydot import graphviz_layout
import pydot 

from oncophylo.ul import root_id 

def show_tree(
    tree,
    muts_as_number=False,
    cells_as_number=False,
    show_id=False,
    cell_info=None,
    save_path=None,
    color_attr=None,
    vertical=True,
    dpi=150,
):
    """Draw the tree
    
    Plots the tree in which edges are mutations and nodes are cells.

    Parameters
    ----------
    tree : :class:`networkx.DiGraph`
        The input tree.
    muts_as_number : :obj:`bool`, optional
        Change the mutation list to a number at edges, by default False
    cells_as_number : :obj:`bool`, optional
        Change the cell list to a number at edges, by default False
    show_id : :obj:`bool`, optional
        Whether to show IDs of nodes and edges or not, by default True
    cell_info : :class:`pandas.DataFrame`, optional
        Information of cells for coloring the nodes by a pie chart, by default None
    save_path : :obj:`str`, optional
        Path to a file for saving the tree in (extension should be .pdf or .png), by default None
    color_attr : :obj:`str`, optional
        Attributes in the `cell_info` dataframe for coloring the nodes, by default None
    dpi : :obj:`int`, optional
        Resolution of rendered figures – this influences the size of
        figures in notebooks, by default 150

    Returns
    -------
    :obj:`None`
    """
    tc = tree.copy()
    root = root_id(tree)
    tc.nodes[root]["label"] = tc.graph["root_name"]
    tc.nodes[root]["fontname"] = "Helvetica"
    tc.nodes[root]["style"] = "rounded"
    tc.nodes[root]["shape"] = "box"
    tc.nodes[root]["margin"] = 0.05
    tc.nodes[root]["pad"] = 0
    tc.nodes[root]["width"] = 0
    tc.nodes[root]["height"] = 0

    if muts_as_number:
        for u, v, label in tc.edges.data("label"):
            if label == "":
                ll = []
            else:
                ll = label.split(tc.graph["splitter_mut"])
            tc.add_edge(u, v, label=f"  {len(ll)}  ")

    if cells_as_number:
        for n in tc.nodes:
            if n != root:
                ll = tc.nodes[n]["label"].split(tc.graph["splitter_cell"])
                if "––" in tc.nodes[n]["label"]:
                    tc.nodes[n]["label"] = "0"
                else:
                    tc.nodes[n]["label"] = f"{len(ll)}"

    if cell_info is not None:
        tc.nodes[root]["label"] = tree.graph["splitter_cell"].join(
            tc.graph["normal_cells"]
        )
        mapping = cell_info[color_attr].to_dict()
        for node in tc:
            num = 0
            paths = nx.shortest_path(tc, source=root, target=node)
            for i in range(len(paths) - 1):
                x = paths[i]
                y = paths[i + 1]
                num += len(tc[x][y]["label"].split(tc.graph["splitter_mut"]))
            try:
                freq = [
                    mapping[x]
                    for x in tc.nodes[node]["label"].split(tc.graph["splitter_cell"])
                ]
            except KeyError:
                freq = ["#FFFFFF"]
            freq = pd.DataFrame(freq)[0].value_counts(normalize=True)
            fillcolor = ""
            for index, value in freq.items():
                fillcolor += f"{index};{value}:"
            tc.nodes[node]["fontsize"] = 14
            tc.nodes[node]["shape"] = "circle"
            tc.nodes[node]["fontname"] = "Helvetica"
            tc.nodes[node]["style"] = "wedged"
            tc.nodes[node]["margin"] = 0.05
            tc.nodes[node]["pad"] = 0
            if "––" in tc.nodes[node]["label"]:
                tc.nodes[node]["width"] = 0
                tc.nodes[node]["height"] = 0
            else:
                tc.nodes[node]["width"] = 0.8
                tc.nodes[node]["height"] = 0.8
            tc.nodes[node]["fillcolor"] = fillcolor 
            tc.nodes[node]["label"] = len(
                tc.nodes[node]["label"].split(tc.graph["splitter_cell"])
            )
            tc.nodes[node]["fontcolor"] = "white"
            tc.nodes[node]["color"] = "gray"

    if show_id:
        for u, v, label in tc.edges.data("label"):
            tc.add_edge(u, v, label=label + f"\n[{v}]")
            tc.nodes[v]["label"] = tc.nodes[v]["label"] + f"\n[{v}]"

    if vertical:
        tc.graph["graph"] = {"fontname": "Helvetica"}
    else:
        tc.graph["graph"] = {"fontname": "Helvetica", "rankdir": "LR"}
    tc.graph["node"] = {"fontname": "Helvetica", "fontsize": 14}
    tc.graph["edge"] = {"fontname": "Helvetica", "fontsize": 14}

    pydot_graph = nx.drawing.nx_pydot.to_pydot(tc)

    if "rank_same" in tc.graph:
        for group in tc.graph["rank_same"]:
            rank_group = pydot.Subgraph(rank='same')
            for node in group:
                rank_group.add_node(pydot.Node(node))
            pydot_graph.add_subgraph(rank_group)

    if save_path:
        if save_path.endswith(".pdf"):
            pydot_graph.write_pdf(save_path)
        elif save_path.endswith(".png"):
            pydot_graph.write_png(save_path)
        else:
            raise ValueError("Unsupported file extension. Use '.pdf' or '.png'.")

    return display(
        Image(pydot_graph.create_png(), embed=True, retina=True)
    )