import os, subprocess, shutil, copy, time
import pandas as pd
import numpy as np
import networkx as nx
from networkx.drawing.nx_pydot import write_dot

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import ListVector

import oncophylo as op 

def ancestry_to_adjacency(ancestry):
    """Simple function to convert an ancestry matrix into an adjaceny matrix"""
    n = ancestry.shape[0]
    adjacency = ancestry.copy()

    for k in range(n):
        for i in range(n):
            for j in range(n):
                if ancestry[i, k] and ancestry[k, j] and i != k and k != j and i != j:
                    adjacency[i, j] = 0

    # Remove self-edges (diagonal)
    np.fill_diagonal(adjacency, 0)

    return adjacency

def LACE(character_matrix,
         cell_samples,
         fp=0.02,
         fn=0.05,
         keep_equivalent=True,
         num_rs=5,
         num_iter=10,
         n_try_bs=5,
         seed=0,
         remove_temp_dir=True,        
         destination_dir=""):

    """
    LACE Python Wrapper (https://github.com/BIMIB-DISCo/LACE)
    
    Parameters
    ----------
    character_matrix: pd.DataFrame
        A dataframe containing the called genotypes (0,1,2,-1) for each SNV/SNP in each cell
    cell_samples: np.array
        A NumPy array containing an integer value indicating which sample each cell came from
    fp: float
        The estimated false positive rate for the sequencing data (SNVs)
    fn: float
        The estimated false negative rate for the sequencing data (SNVs)

    Returns
    --------
    list or dictionary
        A dictionary or a list of dictionaries containing solutions and their corresponding meta data
    """
    temp_path = os.path.join(os.path.abspath(""),"lace_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)

    dot_fn = os.path.join(temp_path, "tree.gv")
    log_file = os.path.join(temp_path, "log.txt")

    # Replace -1/3's with NA, replace 2 with 1 
    modified_character_matrix = character_matrix.copy()
    modified_character_matrix = modified_character_matrix.replace(-1, pd.NA).replace(3, pd.NA).replace(2,1)

    # Create a list of sample names 
    unique_samples, cell_counts = np.unique(cell_samples, return_counts=True)

    # Create a dictionary of per-sample dataframes
    sample_dfs = {
        sample: modified_character_matrix.loc[cell_samples == sample]
        for sample in unique_samples
    }

    r_matrices = []
    for sample, sdf in sample_dfs.items():
        sdf = sdf.astype("Int64")  # nullable integer for NA compatibility
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_df = ro.conversion.py2rpy(sdf)
        # Convert R data.frame to matrix
        r_matrix = ro.r("as.matrix")(r_df)
        # Force mode to integer
        ro.r("storage.mode")(r_matrix).ro = ro.StrVector(["integer"])
        r_matrices.append(r_matrix)

    # data sets for each sample
    r_D = ListVector({f"sample_{i}": mat for i, mat in enumerate(r_matrices)})

    # fraction of cells from each sample
    lik_weights = cell_counts / cell_counts.sum()
    r_lik_weights = ro.FloatVector(lik_weights)

    # false positive/negative rates per sample
    alpha = [[fp]*len(unique_samples)]
    beta = [[fn]*len(unique_samples)]
    r_alpha = ro.ListVector({str(i+1): ro.FloatVector(a) for i, a in enumerate(alpha)})
    r_beta = ro.ListVector({str(i+1): ro.FloatVector(b) for i, b in enumerate(beta)})

    # run LACE
    lace = importr("LACE")

    start_time = time.time() # collect run time
    inference = lace.LACE(D=r_D,
                          lik_w=r_lik_weights,
                          alpha=r_alpha,
                          beta=r_beta,
                          keep_equivalent=keep_equivalent,
                          num_rs=num_rs,
                          num_iter=num_iter,
                          n_try_bs=n_try_bs,
                          num_processes=ro.rinterface.NA_Integer,
                          seed=seed,
                          verbose=True,
                          log_file=log_file)
    runtime = time.time() - start_time

    # collect re-ordering of mutations
    columns = list(inference[inference.names.index("B")].colnames)

    # grab ancestry/adjacency matrix and cell clone assignments
    variants = columns[1:] # variants is all variants without the root
    ancestry_matrix = np.array(inference[inference.names.index("B")]).T
    cell_assignments = np.array(inference[inference.names.index("C")]).reshape(-1).astype(int)
    cells = modified_character_matrix.index.tolist()
    adjacency_matrix = ancestry_to_adjacency(ancestry_matrix)

    # create mutation tree
    T_mut = nx.from_numpy_array(adjacency_matrix, create_using=nx.DiGraph)  
    T_mut.nodes[0]["label"] = "root"
    T_mut.graph["root_name"] = "root"
    T_mut.graph["root"] = "root"
    T_mut.graph["cell_assignments"] = cell_assignments
    T_mut.graph["mutations"] = variants
    T_mut.graph["cells"] = cells
    T_mut.graph["type"] = op.ul.DATA.MUTATION_TREE
    T_mut.graph["losses"] = []
    T_mut.graph["gains"] = []
    T_mut.graph["loss_prefix"] = "-"
    T_mut.graph["gain_prefix"] = "+"
    T_mut.graph["corrected_genotypes"] = np.array(inference[inference.names.index("corrected_genotypes")])
    T_mut.graph["B"] = np.array(inference[inference.names.index("B")])
    T_mut.graph["B_columns"] = list(inference[inference.names.index("B")].colnames)
    T_mut.graph["B_rows"] = list(inference[inference.names.index("B")].rownames)
    T_mut.graph["relative_likelihoods"] = list(inference[inference.names.index("relative_likelihoods")])
    T_mut.graph["C"] = np.array(inference[inference.names.index("C")])
    T_mut.graph["clones_summary"] = [list(l) for l in inference[inference.names.index("clones_summary")]]

    # Relabel the nodes
    mapping = {i+1:var for i,var in enumerate(variants)}
    mapping[0] = "root"

    # create cell tree
    T_mut = nx.relabel_nodes(T_mut, mapping)
    T_cell = copy.deepcopy(T_mut)
    T_cell.graph["type"] = op.ul.DATA.CELL_TREE

    # add cell attachments
    nodes = ["root"] + variants
    for i, node_id in enumerate(cell_assignments):
        T_cell.add_edge(nodes[node_id], cells[i])

    # get inferred genotypes per cell
    (T_cell, corrected_character_matrix) = op.ul.resolve_genotypes(T_cell, character_matrix)

    # save dot file
    write_dot(T_cell, dot_fn)
    op.ul.save_output_files(destination_dir, [dot_fn])

    if remove_temp_dir:
        shutil.rmtree(temp_path)

    return op.ul.solution(T_cell, 
                          T_mut, 
                          character_matrix,
                          corrected_character_matrix, 
                          "",
                          runtime)
