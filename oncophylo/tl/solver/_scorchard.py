import os, subprocess, shutil
import pandas as pd 
import numpy as np 

import oncophylo as op 

def scOrchard(input_df, 
              variant_reads_df=None,
              total_reads_df=None,
              k=10, 
              e=10, 
              fp=0.001, 
              fn=0.2, 
              K=0,
              R=0,
              ado_precision=15.0,
              minCellsL=1,
              minCellsG=1,
              mcmc_iterations=0,
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
        The maximum number of losses per mutation to consider (Default = 0)
    R: int
        The maximum number of gains per mutation to consider (Default = 0)
    ado_precision: float
        The allelic dropout precision for the beta binomial likelihood model. This is only used if read counts are provided.
    minCellsL: int
        The minimum number of cells that benefit from a mutation loss in order to consider the loss. 
        Can be a whole number or a percentage.
    minCellsG: int
        The minimum number of cells that benefit from a mutation gain in order to consider the gain.
        Can be a whole number or a percentage.
    mcmc_iterations: int, optional
        The number of MCMC iterations to perform for each sampled tree (Default = 0)
    seed: int
        The random seed to use
    greedy: bool
        Flag to run scOrchard in its greedy search mode. (Default = False)
    homoplasy_only: bool
        Flag to only consider mutation gains due to homoplasy. (Default = False)
    n_solutions: int
        The number of solutions to return (Default = 1)
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files (Default = True)

    Returns
    --------
    list or dictionary
        A dictionary or a list of dictionaries containing solutions and their corresponding meta data
    """
    
    temp_path = os.path.join(os.path.abspath(""),"orchard_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "input.txt")
    variant_reads_fn = os.path.join(temp_path, "variant_reads.csv")
    total_reads_fn = os.path.join(temp_path, "total_reads.csv")
    cell_names_fn = os.path.join(temp_path, "cell_names.txt")
    gene_names_fn = os.path.join(temp_path, "gene_names.txt")
    clusters_fn = os.path.join(temp_path, "clusters.txt")
    mutation_clusters_fn = os.path.join(temp_path, "mutation_clusters.txt")

    output_path = os.path.join(temp_path)
    output_prefix = "out"
    
    clusters, pairs = [], []
    if constrained:
        if op.ul.CONST.CLUSTER_ID not in input_df.columns:
            raise Exception("%s isn't in the input data frame columns!" % op.ul.CONST.CLUSTER_ID) 
        clusters = input_df[op.ul.CONST.CLUSTER_ID].values
        input_df[op.ul.CONST.CLUSTER_ID].to_csv(clusters_fn, index=False, header=False)
        input_df = input_df.drop(columns=op.ul.CONST.CLUSTER_ID)
        clusterings, _ = op.tl.cluster.KModes(input_df.T, 
                                              list(input_df.index), 
                                              k=len(np.unique(clusters)))
        clusters = [list(input_df.columns[indices]) for indices in clusterings[0]]
        with open(mutation_clusters_fn, "w") as f:
            for cluster in clusters:
                f.write(" ".join(map(str, cluster)) + "\n")

    else:
        if op.ul.CONST.CLUSTER_ID in input_df.columns:
            input_df = input_df.drop(columns=op.ul.CONST.CLUSTER_ID)
        
    input_df.T.to_csv(input_fn, sep=" ", index=False, header=False)
    input_df.index.to_series().to_csv(cell_names_fn, index=False, header=False)
    input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)

    if variant_reads_df is not None and total_reads_df is not None:
        assert variant_reads_df.shape == input_df.shape, "variant read count matrix and character matrix shape mismatch!"
        assert total_reads_df.shape == input_df.shape, "variant read count matrix and character matrix shape mismatch!"
        variant_reads_df.to_csv(variant_reads_fn)
        total_reads_df.to_csv(total_reads_fn)

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
            "-R", "%d" % R,
            "-precision", "%.2f" % ado_precision,
            "-minCellsL", "%d" % minCellsL,
            "-minCellsG", "%d" % minCellsG,
            "-mcmcIterations", "%d" % mcmc_iterations,
            "-n", "%d" % n,
            "-m", "%d" % m,
            "-i", "%s" % input_fn,
            "-o", "%s" % output_path,
            "-p", "%s" % output_prefix,
            "-g", "%s" % gene_names_fn,
            "-c", "%s" % cell_names_fn,
            "-fd", "%.6f" % fp,
            "-ad", "%.6f" % fn,
            "-v", "%s" % variant_reads_fn if variant_reads_df is not None else "",
            "-t", "%s" % total_reads_fn if total_reads_df is not None else "",
            "-seed", "%d" % seed if seed is not None else "",
            "-greedy" if greedy else "",
            "-homoplasyOnly" if homoplasy_only else "",
            "-onlyAddLeaves" if only_add_leaves else "",
            "-s", "%s" % clusters_fn if constrained else "",
            "-l", "%s" % mutation_clusters_fn if constrained else "",
            "-constraint", "%d" % (3 if constrained else 0)
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
                                        time,
                                        var_reads=variant_reads_df,
                                        total_reads=total_reads_df,
                                        ado_precision=ado_precision))
    if remove_temp_dir:
        shutil.rmtree(temp_path)
    if n_solutions == 1:
         return solutions[0]
    else:
        return solutions