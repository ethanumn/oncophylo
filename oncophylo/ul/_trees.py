import numpy as np 
import pandas as pd 
import networkx as nx 

import oncophylo as op
from oncophylo.ul import CONST 

def resolve_genotypes(T, input_df=None):
    """Resolves genotypes given a cell tree
    
    Input
    -------
    T : Networkx.DiGraph
        A cell tree represented with a Networkx Digraph. All leaves and mutations in the input_df must be present in the tree.
    input_df : pd.DataFrame, optional
        The character matrix used as input to reconstruct the tree. Rows must be cells and columns mutations. If not provided,
        then the graph attribute will create its own character matrix containing only the characters and no other data that might
        have been provided with the input.

    Returns
    --------
    Networkx.DiGraph
        The input tree T with an updated graph attribute containing the implied genotypes. Accessible via T.graph["data"]
    pd.DataFrame
        A pandas dataframe containing the genotypes for each cell implied by the input tree T. The index/column names of this matrix 
        will be the same as input_df.
    """
    assert T.graph["type"] == CONST.CELL_TREE, "Tree must be of type %s not %s" % (CONST.CELL_TREE, T.graph["type"])
    
    if isinstance(input_df, pd.DataFrame):
        output_df = input_df.copy()
    else:
        output_df = pd.DataFrame(0, index=T.graph["cells"], columns=T.graph["mutations"])
        
    for col in output_df.columns:
        if col in T.graph["mutations"]:
            output_df[col] = 0
        
    # find all cells with the same parent, they have the same genotypes
    cell_parents = {}
    for s in T.graph["cells"]:
        am = list(T.predecessors(s))[0]
        if am in cell_parents:
            cell_parents[am].append(s)
        else:
            cell_parents[am] = [s]
        
    # go through each group of cells, find the genotype for one of them and you find it for all
    for am, cells in cell_parents.items():

        seen = set() # keep track of mutations we've already seen

        while am != T.graph["root_name"]: # going from leaves up
            if am in T.graph["losses"]: # detect mutation losses
                lost_mutation = am.lstrip(T.graph["loss_prefix"]).split(" ")[0]
                if lost_mutation in seen:
                    pass
                else:
                    output_df.loc[cells, lost_mutation] = 2
                    seen.add(lost_mutation)
            elif am in T.graph["gains"]: # detect homomplasy and recurrence
                gained_mutation = am.lstrip(T.graph["gain_prefix"]).split(" ")[0]
                if gained_mutation in seen:
                    pass
                else:
                    output_df.loc[cells, gained_mutation] = 4
                    seen.add(gained_mutation)

            elif am in T.graph["mutations"]:
                if am in seen:
                    pass
                else:
                    output_df.loc[cells, am] = 1
                    seen.add(am)
            else:
                pass
            am = list(T.predecessors(am))[0]

    T.graph["data"] = output_df
    
    return T, output_df

def to_clonal_tree(T, df):
    """Converts a mutation tree with cell attachments (i.e., a cell tree) to a clonal tree."""
     # copy because we're first going to remove any mutations that are leaves without cells attached
    if not isinstance(T, nx.DiGraph):
        print("Input tree is not a Networkx DiGraph, returning an empty tree...")
        return nx.DiGraph()
    
    T_ = T.copy()
    leaf_nodes = [x for x in T_.nodes() if T_.out_degree(x)==0 and T_.in_degree(x)==1]
    to_remove = list(filter(lambda x: x in list(df.columns), leaf_nodes))
    T_.remove_nodes_from(to_remove)
    
    tree = nx.DiGraph()
    tree.graph["data"] = df
    tree.graph["splitter_mut"] = "\n"
    tree.graph["splitter_cell"] = "\n"
    tree.graph["become_germline"] = list(df.columns[(df == 0).all(axis=0)])
    tree.graph["type"] = CONST.CLONAL_TREE
    node_id = 0
    
    # store list of cells and nodes for future reference
    node_list = list(T_.nodes())
    tree.graph["cells"] = list(df.index)
    tree.graph["mutations"] = list(df.columns)
    tree.graph["losses"] = sorted(filter(lambda n: n.startswith(T_.graph["loss_prefix"]), node_list))
    tree.graph["gains"] = sorted(filter(lambda n: n.startswith(T_.graph["gain_prefix"]), node_list))
    tree.graph["normal_cells"] = list(df[df.sum(axis=1) == 0].index)
    tree.graph["root_id"] = 0
    tree.graph["loss_prefix"] = T_.graph["loss_prefix"]
    tree.graph["gain_prefix"] = T_.graph["gain_prefix"]

    tree.add_node(node_id)
    tree.nodes[node_id]["label"] = "root"
    
    if "root_name" in T_.graph:
        root_name = T_.graph["root_name"]
    else:
        roots = list(filter(lambda p: p[1] == 0, T_.in_degree()))
        assert 1 == len(roots)
        root_name = roots[0][0]

    queue = [(None,root_name,"")]
    node_id = 1
    parent_id = 0
    edge_name = ""
    node_name = ""

    node_id_list = [0]
    node_names = ["root"]
    edges = []
    edge_names = []

    while len(queue) > 0:

        # grab next node off of queue
        parent_id, node, edge_name = queue.pop(0) 
        edges = list(nx.bfs_tree(T_, source=node, depth_limit=1).edges)
        node_name = ""
        make_node = False

        # initialize edge name 
        if node != "root":
            if len(edge_name) == 0:
                edge_name = node
            else:
                edge_name = edge_name + tree.graph["splitter_mut"] + node

        if parent_id == None:
            node_id = node_id_list[0]
        elif len(edges) > 1 or any([e[1] in tree.graph["cells"] for e in edges]):
            node_id = np.max(node_id_list) + 1
            node_id_list.append(node_id)
            make_node = True
        else:
            node_id = parent_id

        # add all mutations to queue and label this node by all cells attached
        for e in edges:

            # collect all cells attached to this mutation to label this node
            if e[1] in tree.graph["cells"]:
                if len(node_name) == 0:
                    node_name = e[1]
                else:
                    node_name = node_name + tree.graph["splitter_cell"] + e[1]

            # add all mutations to queue
            elif (e[1] in tree.graph["mutations"]) or (e[1] in tree.graph["losses"]) or (e[1] in tree.graph["gains"]):
                if len(edges) == 1:
                    queue.append((node_id,e[1],edge_name))
                else:
                    queue.append((node_id,e[1],""))

        # make node
        if make_node:
            if len(node_name) == 0:
                node_name = "––"
            tree.add_node(node_id)
            tree.nodes[node_id]["label"] = node_name
            tree.add_edge(parent_id, node_id, label=edge_name)

    
    return tree


def conflict_free_matrix_to_clonal_tree(df):
    """Convert a CONFLICT FREE cell x mutation matrix to a Networkx.DiGraph.

    This function is a modification of https://github.com/faridrashidi/scphylo-tools/blob/main/scphylo/ul/_trees.py
    and is primarily used to process HUNTRESS's output.

    Parameters
    ----------
    df: pandas.DataFrame
        A cell x mutation matrix from which a perfect phylogeny can be resolved via Gusfield's algorithm

    Returns
    -------
    networkx.DiGraph
        A perfect phylogeny
    """
    if not op.ul.is_conflict_free_gusfield(df):
        print("Matrix is not conflict free!")

    def _contains(col1, col2):
        for i in range(len(col1)):
            if not col1[i] >= col2[i]:
                return False
        return True

    tree = nx.DiGraph()
    tree.graph["data"] = df
    tree.graph["splitter_mut"] = "\n"
    tree.graph["splitter_cell"] = "\n"
    tree.graph["become_germline"] = list(df.columns[(df == 0).all(axis=0)])
    tree.graph["type"] = CONST.CLONAL_TREE
    tree.graph["cells"] = list(df.index.values)
    tree.graph["mutations"] = list(df.columns.values)

    # these are just dummy values
    tree.graph["gain_prefix"] = "+"
    tree.graph["loss_prefix"] = "-"
    tree.graph["losses"] = []
    tree.graph["gains"] = []
    
    matrix = df.values
    names_mut = list(df.columns)

    i = 0
    while i < matrix.shape[1]:
        j = i + 1
        while j < matrix.shape[1]:
            if np.array_equal(matrix[:, i], matrix[:, j]):
                matrix = np.delete(matrix, j, 1)
                x = names_mut.pop(j)
                names_mut[i] += tree.graph["splitter_mut"] + x
                j -= 1
            j += 1
        i += 1

    cols = matrix.shape[1]
    dimensions = np.sum(matrix, axis=0)
    indices = np.argsort(dimensions)
    dimensions = np.sort(dimensions)
    names_mut = [names_mut[indices[i]] for i in range(cols)]

    tree.add_node(cols)
    tree.add_node(cols - 1)
    tree.add_edge(cols, cols - 1, label=names_mut[cols - 1])
    node_mud = {}
    node_mud[names_mut[cols - 1]] = cols - 1

    i = cols - 2
    while i >= 0:
        if dimensions[i] == 0:
            break
        attached = False
        for j in range(i + 1, cols):
            if _contains(matrix[:, indices[j]], matrix[:, indices[i]]):
                tree.add_node(i)
                tree.add_edge(node_mud[names_mut[j]], i, label=names_mut[i])
                node_mud[names_mut[i]] = i
                attached = True
                break
        if not attached:
            tree.add_node(i)
            tree.add_edge(cols, i, label=names_mut[i])
            node_mud[names_mut[i]] = i
        i -= 1

    tumor_cells = []
    clusters = {cols: "root"}
    for node in tree:
        if node == cols:
            tree.nodes[node]["label"] = "root"
            continue
        untilnow_mut = []
        sp = nx.shortest_path(tree, cols, node)
        for i in range(len(sp) - 1):
            untilnow_mut += tree.get_edge_data(sp[i], sp[i + 1])["label"].split(
                tree.graph["splitter_mut"]
            )
        untilnow_cell = df.loc[
            (df[untilnow_mut] == 1).all(axis=1)
            & (df[[x for x in df.columns if x not in untilnow_mut]] == 0).all(axis=1)
        ].index
        if len(untilnow_cell) > 0:
            clusters[node] = f"{tree.graph['splitter_cell'].join(untilnow_cell)}"
            tumor_cells += list(
                y for y in tree.graph["splitter_cell"].join(untilnow_cell)
            )
        else:
            clusters[node] = "––"

        tree.nodes[node]["label"] = clusters[node]

    tree.graph["root_id"] = cols

    i = 1
    for k, v in clusters.items():
        if v == "––":
            clusters[k] = i * "––"
            i += 1
    return tree

def is_conflict_free_gusfield(df):
    """Check conflict-free criteria via Gusfield algorithm.

    Parameters
    ----------
    df: pd.DataFrame
        A cell x mutation genotype matrix

    Returns
    -------
    bool
        True if the matrix is conflict free, False otherwise
    """
    I_mtr = df.astype(int).values
    if not np.array_equal(np.unique(I_mtr), [0, 1]):
        return False

    def _sort_bin(a):
        b = np.transpose(a)
        b_view = np.ascontiguousarray(b).view(
            np.dtype((np.void, b.dtype.itemsize * b.shape[1]))
        )
        idx = np.argsort(b_view.ravel())[::-1]
        c = b[idx]
        return np.transpose(c), idx

    Ip = I_mtr.copy()
    O_mtr, _ = _sort_bin(Ip)
    Lij = np.zeros(O_mtr.shape, dtype=int)
    for i in range(O_mtr.shape[0]):
        maxK = 0
        for j in range(O_mtr.shape[1]):
            if O_mtr[i, j] == 1:
                Lij[i, j] = maxK
                maxK = j + 1
    Lj = np.amax(Lij, axis=0)
    for i in range(O_mtr.shape[0]):
        for j in range(O_mtr.shape[1]):
            if O_mtr[i, j] == 1:
                if Lij[i, j] != Lj[j]:
                    return False
    return True

def clonal_to_cell_tree(T_clonal):
    """Converts a clonal tree to a cell tree. 
    
    TODO:
        This is working well, but there's some issue with how Networkx sets the node ID's.
        This is more of a minor technical issue that can be reproduced by performing the following:

            cell_tree = solution[op.ul.CONST.CELL_TREE]
            clonal_tree = solution[op.ul.CONST.CLONAL_TREE]
            cell_tree_duplicate = op.ul.clonal_to_cell_tree(clonal_tree)
            clonal_tree_duplicate = op.ul.to_clonal_tree(cell_tree_duplicate)

            nx.utils.graphs_equal(cell_tree, cell_tree_duplicate) # should be True, but returns False
            nx.utils.graphs_equal(clonal_tree, clonal_tree_duplicate) # should be True, but returns False

        This occurs only sometimes (with larger trees especially), and I believe it's because translating between
        cell/clonal trees with these functions changes the numeric ID values for nodes. Consequently, nx.utils.graphs_equal 
        will find that the .edges attributes are not the same, and therefore returns False.

        However, I do not think this impacts comparing trees because the labels match exactly in all the tests I've done.

    Parameters
    ----------
    T_clonal: networkx.DiGraph
        A clonal tree where nodes are labeled with cells and edges are labeled with mutations

    Returns
    -------
    networkx.DiGraph
        A cell tree where 
    """
    T_cell = nx.DiGraph()
    T_cell.add_node("root", label="root")
    T_cell.graph["type"] = CONST.CELL_TREE
    
    # copy all relevant data
    T_cell.graph["losses"] = T_clonal.graph["losses"]
    T_cell.graph["gains"] = T_clonal.graph["gains"]
    T_cell.graph["cells"] = T_clonal.graph["cells"]
    T_cell.graph["mutations"] = T_clonal.graph["mutations"] 
    T_cell.graph["become_germline"] = T_clonal.graph["become_germline"]
    T_cell.graph["data"] =  T_clonal.graph["data"] 
    T_cell.graph["loss_prefix"] = T_clonal.graph["loss_prefix"]
    T_cell.graph["gain_prefix"] = T_clonal.graph["gain_prefix"]

    # go through all of the edges
    for u, v, l in T_clonal.edges.data("label"):
        
        # split apart all mutations on this edge
        muts = l.split(T_clonal.graph["splitter_mut"])
        
        # base case where u is the root
        if T_clonal.in_degree(u) == 0:
            T_cell.add_node(muts[0], label=muts[0])
            T_cell.add_edge("root", muts[0])
                
        # go through and add edges between mutations
        last_mut = muts[0]
        if last_mut not in T_cell.nodes:
            T_cell.add_node(last_mut, label=last_mut)
        
        for _mut in muts[1:]:
            T_cell.add_node(_mut, label=_mut)
            T_cell.add_edge(last_mut, _mut)
            last_mut = _mut
           
        # now attach all cells to last mutation
        for cell in T_clonal.nodes[v]["label"].split(T_clonal.graph["splitter_cell"]):
            if cell == "––":
                continue
            T_cell.add_node(cell, label=cell)
            T_cell.add_edge(last_mut, cell)
                
        # find edges that are children of this edge (u,v) and add an edge from the last mutation on this edge 
        # to the first mutation on all child edges
        child_edges = nx.dfs_edges(T_clonal, source=v, depth_limit=1)
        for (w,y) in child_edges:
            l = T_clonal[w][y].get('label', None)
            # go through each edge, and all an edge from the last mutation in the edge (u,v) to
            # the first mutation on the edge e
            first_mut = l.split(T_clonal.graph["splitter_mut"])[0]
            T_cell.add_node(first_mut, label=first_mut)
            T_cell.add_edge(last_mut, first_mut)

    T_mut = T_cell.copy()
    T_mut.graph["type"] = CONST.MUTATION_TREE
    for cell in T_mut.graph["cells"]:
        if cell in T_mut.nodes():
            T_mut.remove_node(cell)
                
            
    return T_cell, T_mut

def root_id(tree):
    """Finds the root in a Networkx tree
    
    Inputs
    -------
    tree: Networkx.DiGraph
        A mutation, cell, or clonal tree

    Returns
    --------
    object
        The root id. This could be a string or integer.
    """
    for x in tree.nodes:
        if tree.in_degree(x) == 0:
            return x
    return None