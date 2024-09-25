# _lt_utils.py

import pandas as pd 
import numpy as np
import cassiopeia as cas
import networkx as nx

def load_lineage_tracing_file(fn, min_cells=1):
    """Loads a lineage tracing file in the same format that Cassiopeia uses"""
    raw_df = pd.read_csv(fn, sep="\t", index_col=0, header=0)
    return prep_lineage_tracing_data(raw_df, min_cells)

def prep_lineage_tracing_data(character_matrix, min_cells, missing_character="-"):
    """Prepares a lineage tracing indel file by converting it to a binary representation"""
    mapping = {barcode:"s%d" % i for i,barcode in enumerate(list(character_matrix.index))}
    character_matrix = character_matrix.rename(index=mapping,
              columns={r:"indel%d" % (i+1) for i,r in enumerate(list(character_matrix.columns))})
    character_matrix.index.names = ["cells"]
    
    binary_character_matrix = to_binary_df(character_matrix, missing_character)

    binary_character_matrix[binary_character_matrix == "0"] = 0
    binary_character_matrix[binary_character_matrix == "-"] = -1
    binary_character_matrix = binary_character_matrix.astype(int)
    
    input_df = binary_character_matrix.loc[:,(binary_character_matrix > 0).sum(axis=0) > min_cells]
    return input_df, binary_character_matrix,  mapping

def to_binary_df(df, missing_character="-", unknown_character=3):
    """Converts a lineage tracing indel file to a binary representation
    
    Parameters
    -----------
    df: pandas.DataFrame    
        A data frame where rows are cells and columns are indels
    missing_character: str, optional
        The character which represents a missing value for a cutsite
    missing_character_numeric: int, optional
        The numeric value to convert the missing characters in the indel data frame to

    Returns
    --------
    pandas.DataFrame    
        A binary representation of the original input data where rows are cells and columns are indel-cutsite pairs. A 1 means the cell has a particular
        indel at that cutsite, and a 0 means the particular indel did not occur at that cutsite. Each cutsite can only have a single "1" in one of the possible indels.
    """
    characters_to_skip = ["0", 0, missing_character]
    columns = list(df.columns)   
    index = df.index
    size = len(index)
    data = {}
    
    for col in columns:
        unique_indels = np.unique(df[col])
        unknown = df[col] == "-"
        for i,indel in enumerate(unique_indels):
            if indel in characters_to_skip: # no need to add unknown or absent as a column
                continue
            indel_col = "%s_%s" % (col, indel)
            data[indel_col] = np.zeros(size)
            data[indel_col][unknown] = unknown_character
            data[indel_col][df[col] == indel] = 1
            
    df_b = pd.DataFrame(data, index=df.index)

    return df_b.astype(int)    

def post_process_celltree(T, character_matrix, muts_to_add, mapping):
    """
    Post processes a cell tree output by a single-cell phylogeny reconstruction method into a Cassiopeia tree

    Parameters
    -----------
    T: Networkx.DiGraph
        A directed tree where internal nodes represent indels at particular cutsites, and leaf nodes are cells
    character_matrix: pandas.DataFrame
        The original character matrix where rows are cells and column are cutsites, and entries at the type of indel that occurred 
        at the cutsite
    muts_to_add: pandas.DataFrame
        A data frame where rows are cells and columns are indel-cutsite pairs. These indel-cutsite pairs must occur in only a single cell, and
        are added to the tree as a direct parent of the single cells they occur in.
    mapping: dict
        A dictionary for relabeling nodes in the tree

    Returns
    --------
    cas.data.CassiopeiaTree
        A CassiopeiaTree object containing the phylogeny and character matrix
    """
    for c in list(muts_to_add.index):
        muts = list(muts_to_add.loc[c, muts_to_add.loc[c] == 1].index)

        # each mutation only present in this cell c
        for m in muts:

            # find parent of cell
            am = list(T.predecessors(c))[0]

            # remove edge from parent to cell
            T.remove_edge(am, c)

            # add edge from parent to mutation
            T.add_edge(am, m)

            # add edge from mutation to cell
            T.add_edge(m, c)
            
    Tnew = nx.DiGraph(T)
    Tnew = nx.relabel_nodes(Tnew, {j:i for i,j in mapping.items()})
    
    # remove leaf nodes that have no cells attached
    leaf_nodes = [x for x in Tnew.nodes() if Tnew.out_degree(x)==0 and Tnew.in_degree(x)==1]
    for n in leaf_nodes:
        if n.startswith("indel"):
            Tnew.remove_node(n)
            
    return cas.data.CassiopeiaTree(tree=Tnew, character_matrix = character_matrix)