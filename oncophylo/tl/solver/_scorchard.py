import os, subprocess, shutil
import pandas as pd 
import numpy as np 

import oncophylo as op 

def scOrchard(input_df, 
              k=10, 
              e=10, 
              fp=0.001, 
              fn=0.2, 
              K=0,
              J=0,
              minCellsL=1,
              minCellsG=1,
              seed=None, 
              homoplasy_only=False,
              only_add_leaves=False,
              constrained=False,
              greedy=False,
              n_solutions=1, 
              remove_temp_dir=True):
    """
    scOrchard Python Wrapper
    
    Parameters
    ----------
    input_df: pd.DataFrame
        A Pandas dataframe where the row are cells and the columns are mutations
    k: int
        The beam width -- which controls the number of unique tree returned
    e: int 
        The branching factor -- which controls how many solutions are expanded from the queue at each iteration
    fp: float
        The false positive rate to use for likelihood calculations
    fn: float
        The false negative rate to use for likelihood calculations
    K: int
        The maximum number of losses per mutation to consider. Default = 0
    J: int
        The maximum number of gains per mutation to consider. Default = 0
    minCellsL: int
        The minimum number of cells that benefit from a mutation loss in order to consider the loss. 
        Can be a whole number or a percentage.
    minCellsG: int
        The minimum number of cells that benefit from a mutation gain in order to consider the gain.
        Can be a whole number or a percentage.
    seed: int
        The random seed to use
    greedy: bool
        Flag to run scOrchard in its greedy search mode. Default = False
    homoplasy_only: bool
        Flag to only consider mutation gains due to homoplasy. Default = False
    n_solutions: int
        The number of solutions to return. Default = 1
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files. Default = True

    Returns
    --------
    list or dictionary
        A dictionary or a list of dictionaries containing solutions and their corresponding meta data
    """
    
    temp_path = os.path.join(os.path.abspath(""),"orchard_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "input.txt")
    cell_names_fn = os.path.join(temp_path, "cell_names.txt")
    gene_names_fn = os.path.join(temp_path, "gene_names.txt")
    clusters_fn = os.path.join(temp_path, "clusters.txt")
    loss_pairs_fn = os.path.join(temp_path, "loss_pairs.txt")

    output_path = os.path.join(temp_path)
    output_prefix = "out"
    
    clusters, loss_pairs = [], []
    if constrained:
        if op.ul.CONST.CLUSTER_ID not in input_df.columns:
            raise Exception("%s isn't in the input data frame columns!" % op.ul.CONST.CLUSTER_ID) 
        clusters = input_df[op.ul.CONST.CLUSTER_ID].values
        input_df[op.ul.CONST.CLUSTER_ID].to_csv(clusters_fn, index=False, header=False)
        input_df = input_df.drop(columns=op.ul.CONST.CLUSTER_ID)
        loss_pairs = op.ul.find_loss_pairs(input_df, list(input_df.columns), np.unique(clusters).size)

    else:
        if op.ul.CONST.CLUSTER_ID in input_df.columns:
            input_df = input_df.drop(columns=op.ul.CONST.CLUSTER_ID)
        
    input_df.T.to_csv(input_fn, sep=" ", index=False, header=False)
    input_df.index.to_series().to_csv(cell_names_fn, index=False, header=False)
    input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)
    pd.DataFrame(loss_pairs).to_csv(loss_pairs_fn, index=False, header=False, sep=" ")

    n, m = input_df.shape
    
    # process minCellsL and minCellsG
    if minCellsL < 1:
        minCellsL = np.maximum(1, np.ceil(minCellsL*n))

    if minCellsG < 1:
        minCellsG = np.maximum(1, np.ceil(minCellsG*n))    

    args = [
            "-k", "%s" % k,
            "-e", "%d" % e,
            "-K", "%d" % K,
            "-J", "%d" % J,
            "-minCellsL", "%d" % minCellsL,
            "-minCellsG", "%d" % minCellsG,
            "-n", "%d" % n,
            "-m", "%d" % m,
            "-i", "%s" % input_fn,
            "-o", "%s" % output_path,
            "-p", "%s" % output_prefix,
            "-g", "%s" % gene_names_fn,
            "-c", "%s" % cell_names_fn,
            "-fd", "%.6f" % fp,
            "-ad", "%.6f" % fn,
            "-seed", "%d" % seed if seed is not None else "",
            "-greedy" if greedy else "",
            "-homoplasyOnly" if homoplasy_only else "",
            "-onlyAddLeaves" if only_add_leaves else "",
            "-s", "%s" % clusters_fn if constrained else "",
            "-l", "%s" % loss_pairs_fn if constrained else ""
    ]

    # run scOrchard
    output, time = op.ul.subprocess([os.path.join(os.path.abspath(""), "../sgts/bin/sgts")] + args)

    solutions = []
    for i in range(n_solutions):  
        dot_fn = os.path.join(output_path, output_prefix + "_ml%d.gv" % i)
        T_cell, T_mut = op.io.load_dot(dot_fn, 
                                       mutations = list(input_df.columns), 
                                       cells = list(input_df.index), 
                                       _type="cell_tree")
        (T_cell, output_df) = op.ul.resolve_genotypes(T_cell, input_df)
        solutions.append(op.ul.solution(T_cell, 
                                        T_mut, 
                                        input_df,
                                        output_df, 
                                        fp,
                                        fn,
                                        output,
                                        time))
    if remove_temp_dir:
        shutil.rmtree(temp_path)
    if n_solutions == 1:
         return solutions[0]
    else:
        return solutions