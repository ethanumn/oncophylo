import re 
import networkx as nx 
from oncophylo.ul import DATA

# load dot file 
def load_dot(fn, 
             mutations = [], 
             cells = [], 
             mapping = {},
             _type="", 
             loss_prefix="-", 
             gain_prefix="+",
             set_id_to_label = False):
    """Load dot file. This function is heavily customized to be able to load extra data from the dot file stored
    in the 'graph' attribute.
    
    Input
    ------
    fn: str
        The path to the dot file
    mutations: list
        A list of mutation names. If it's empty, the mutation names will attempt to be read from the graph attribute of the tree.
    cells: list
        A list of cell names. If it's empty, the cell names will attempt to be read from the graph attribute of the tree.
    _type: str
        The type of tree being read from the file (clonal_tree, mutation_tree, cell_tree). These are all defined in the DATA dataclass (i.e., DATA.CLONAL_TREE, DATA.MUTATION_TREE, DATA.CELL_TREE)
    loss_prefix: str
        The prefix that designates a loss in the tree. Default = '-'
    gain_prefix: str
        The prefix that designates a gain in the tree. This is specific to recurrent mutations. Default = '+'
    set_id_to_label: bool
        Copies each node's ID to it's label attribute. This is used when a tree's nodes don't have a human readable label attribute, but their ID is human readable. Default = False

    Returns
    --------
    Networkx.DiGraph
        A clonal/cell/mutation tree read from the file
    Networkx.DiGraph
        A mutation tree. This will only be defined if the input is a cell tree. Otherwise, this is will be None.
    """
    T = nx.drawing.nx_pydot.read_dot(fn)

    # Extract the raw DOT source
    with open(fn, 'r') as file:
        dot_source = file.read()

    # Regular expression to match the graph attributes section
    graph_attr_pattern = re.compile(r'graph \[(.*?)\];', re.DOTALL)
    match = graph_attr_pattern.search(dot_source)

    if match:
        # Extract the attributes string
        attr_string = match.group(1)

        # Split the attributes and evaluate them
        attr_pairs = re.split(r',\s*(?=[a-zA-Z_]+\s*=)', attr_string)
        for pair in attr_pairs:
            key, value = pair.split('=', 1)
            key = key.strip()
            value = value.strip().strip('"')

            # Evaluate the value if it looks like a list or dict
            if value.startswith('[') and value.endswith(']'):
                value = eval(value)

            # Set the attribute in the NetworkX graph
            T.graph[key] = value
            
    # if node ids are integers, change it so the ids are their names
    if set_id_to_label:
        id_to_label = {}
        for n in T.nodes():
            id_to_label[n] = T.nodes()[n]["label"].replace("\"", "").replace("\'", "")
        T = nx.relabel_nodes(T, id_to_label)
            
    # if we need to relabel nodes
    if len(mapping) > 0:
        
        # set node label and actual name using mapping, label attribute needs to be changed first
        nx.set_node_attributes(T, mapping, 'label')
        T = nx.relabel_nodes(T, mapping)
        
        # and relabel edges
        new_edge_labels = {}
        for u, v, data in T.edges(data=True):
            old_label = data.get('label', None)
            if old_label:
                # Update edge labels if needed, e.g., keeping the same label or modifying it
                new_edge_labels[(u, v, 'label')] = old_label

        # Apply the new edge labels
        nx.set_edge_attributes(T, new_edge_labels)
                
    # do this after renaming things
    roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
    assert 1 == len(roots)
            
    # use provided attributes if they aren't present in the dot file
    if "root_name" not in T.graph:
        if roots[0][0] not in mutations:
            T.graph["root_name"] = roots[0][0]
    if len(_type) > 0:
        T.graph["type"] = _type
    if len(cells) > 0:
        T.graph["cells"] = cells
    if len(mutations) > 0:
        T.graph["mutations"] = mutations
    if "losses" not in T.graph and len(loss_prefix) > 0:
        T.graph["loss_prefix"] = loss_prefix
        T.graph["losses"] = sorted(filter(lambda x: x.startswith(loss_prefix), T.nodes()))
    if "gains" not in T.graph and len(gain_prefix) > 0:
        T.graph["gain_prefix"] = gain_prefix
        T.graph["gains"] = sorted(filter(lambda x: x.startswith(gain_prefix), T.nodes()))
        
    T_prime = None
    # mutation tree is recovered by removing all cells
    if _type == DATA.CELL_TREE:
        T_prime = T.copy()
        T_prime.graph["type"] = DATA.MUTATION_TREE
        for cell in T_prime.graph["cells"]:
            if cell in T_prime.nodes():
                T_prime.remove_node(cell)
                

    return T, T_prime