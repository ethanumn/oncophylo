import numpy as np 
import pandas as pd 
import networkx as nx 

from oncophylo.ul import CONST 

# conflict_free.py
def resolve_genotypes(T, input_df):
    """Resolves genotypes given a cell tree
    
    Input
    -------
    T : Networkx.DiGraph
        A cell tree represented with a Networkx Digraph. All leaves and mutations in the input_df must be present in the tree.
    input_df : pd.DataFrame
        The character matrix used as input to reconstruct the tree. Rows must be cells and columns mutations.

    Returns
    --------
    Networkx.DiGraph
        The input tree T with an updated graph attribute containing the implied genotypes. Accessible via T.graph["data"]
    pd.DataFrame
        A pandas dataframe containing the genotypes for each cell implied by the input tree T. The index/column names of this matrix 
        will be the same as input_df.
    """
    assert T.graph["type"] == CONST.CELL_TREE, "Tree must be of type %s not %s" % (CONST.CELL_TREE, T.graph["type"])
    
    pred_df = input_df.copy()
        
    for col in pred_df.columns:
        if col in T.graph["mutations"]:
            pred_df[col] = 0
        
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

        lost = []
        gained = []

        while am != T.graph["root_name"]: # going from leaves up
            if am in T.graph["losses"]:
                lost_mutation = am.lstrip(T.graph["loss_prefix"]).split(" ")[0]
                if lost_mutation in gained: # recurrence
                    pred_df.loc[cells, lost_mutation] = 5
                else:
                    pred_df.loc[cells, lost_mutation] = 2
                lost.append(lost_mutation)
            elif am in T.graph["gains"]: # homoplasy
                gained_mutation = am.lstrip(T.graph["gain_prefix"]).split(" ")[0]
                if gained_mutation in lost:
                    pass
                else:
                    pred_df.loc[cells, gained_mutation] = 4
                    gained.append(gained_mutation)

            elif am in T.graph["mutations"] and (am not in lost) and (am not in gained):
                pred_df.loc[cells, am] = 1
            else:
                pass
            am = list(T.predecessors(am))[0]

    T.graph["data"] = pred_df
    
    return T, pred_df

def to_clonal_tree(T, df):
    """Converts a mutation tree with cell attachments to a clonal tree."""
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
    tree.graph["losses"] = list(filter(lambda n: n.startswith(T_.graph["loss_prefix"]), node_list))
    tree.graph["gains"] = list(filter(lambda n: n.startswith(T_.graph["gain_prefix"]), node_list))
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


def matrix_to_clonal_tree(df):
    """Convert a conflict-free matrix to a tree object.

    This function converts a conflict-free matrix to a tree object in which
    nodes are labeled with cells and edges are lables with mutations. The root is
    labled by 'root'. Mutations are seperated by `.graph['splitter_mut']` and cells
    are seperated by `.graph['splitter_cell']`. Those mutations that are not present
    in any cell are stored in `.graph['become_germline']`. Mutations happed once
    during the evolution so there is no repetitive mutation.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        A genotype dataframe in which rows are cells and columns are mutations.
        Note that this dataframe must be conflict-free.

    Returns
    -------
    :class:`networkx.DiGraph`
        A perfect phylogenetic tree.
    """
    if not is_conflict_free_gusfield(df):
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
    tree.graph["losses"] = []
    
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

    # rows = matrix.shape[0]
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

    tree.graph["normal_cells"] = df[df.sum(axis=1) == 0].index
    tree.graph["root_id"] = cols

    i = 1
    for k, v in clusters.items():
        if v == "––":
            clusters[k] = i * "––"
            i += 1
    return tree

def is_conflict_free_gusfield(df_in):
    """Check conflict-free criteria via Gusfield algorithm.

    This is an implementation of algorithm 1.1 in :cite:`Gusfield_1991`.

    The order of this algorithm is :math:`O(nm)`
    where n is the number of cells and m is the number of mutations.

    Parameters
    ----------
    df_in : :class:`pandas.DataFrame`
        Input genotype matrix.

    Returns
    -------
    :obj:`bool`
        A Boolean checking if the input conflict-free or not.

    Examples
    --------
    >>> sc = scp.datasets.test()
    >>> scp.ul.is_conflict_free_gusfield(sc)
    False

    See Also
    --------
    :func:`scphylo.ul.is_conflict_free`.
    """
    I_mtr = df_in.astype(int).values
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


def clonal_to_mutation_tree(tree):
    """Converts a clonal tree to a mutation tree.

    Parameters
    ----------
    tree : :class:`networkx.DiGraph`
        The phylogenetic tree in which cells are in nodes and
        mutations are at edges.

    Returns
    -------
    :class:`networkx.DiGraph`
        The mutation tree in which mutations are in nodes.
    """
    mutation_tree = nx.DiGraph()
    mutation_tree.graph["type"] = CONST.CELL_TREE
    for u, v, l in tree.edges.data("label"):
        if tree.in_degree(u) == 0:
            mutation_tree.add_node(u, label="root")
        muts = l.split(tree.graph["splitter_mut"])
        mutation_tree.add_node(v, label=muts)
        mutation_tree.add_edge(u, v)
    return mutation_tree

    def _clonal_cell_mutation_list(tree):
        muts_list = []
        cells_list = []
        for _, v, l in tree.edges.data("label"):
            muts = l.split(tree.graph["splitter_mut"])
            if "––" not in tree.nodes[v]["label"]:
                cells = tree.nodes[v]["label"].split(tree.graph["splitter_cell"])
                for mut in muts:
                    muts_list.append({"mut": mut, "Node": f"[{v}]"})
                for cell in cells:
                    cells_list.append({"cell": cell, "Node": f"[{v}]"})
            else:
                for mut in muts:
                    muts_list.append({"mut": mut, "Node": f"[{v}]"})
        cells_list = pd.DataFrame(cells_list).set_index("cell")
        muts_list = pd.DataFrame(muts_list).set_index("mut")
        return cells_list, muts_list


def root_id(tree):
    """Finds the root in a Networkx tree
    
    Inputs
    -------
    tree : Networkx.DiGraph
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