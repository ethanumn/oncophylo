import os, subprocess, shutil
import pandas as pd 
import numpy as np 

import oncophylo as op 

def scOrchard(input_df, 
              variant_reads_df=None,
              total_reads_df=None,
              region_reads_df=None,
              meta_df=None,
              region_weights_df=None,
              k=10, 
              e=10, 
              fp=0.001, 
              fn=0.2, 
              K=0,
              R=0,
              ado_precision=15.0,
              min_clone_fractions=[0.03],
              hc_iterations=0,
              max_restarts=1,
              patience=0,
              seed=None, 
              hill_climbing_only=False,
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
        The estimated false positive rate for the sequencing data
    fn: float
        The estimated false negative rate for the sequencing data
    K: int
        The maximum number of losses per mutation to consider (Default = 0)
    R: int
        The maximum number of gains per mutation to consider (Default = 0)
    ado_precision: float
        The allelic dropout precision for the beta binomial likelihood model. This is only used if read counts are provided.
    min_clone_fractions: list
        A list of floats where each describes the minimum fraction of cells that can be used to define a new clonal population. Default is [0.03]
    hc_iterations: int, optional
        The number of hill climbing iterations to perform for each sampled tree (Default = 0)
    max_restarts: int, optional
        The maximum number of times to run scOrchard in an effort to learn a better mutation order (Default = 1)
    patience: int, optional
        The number of consecutive restarts scOrchard will perform unless the log-likelihood of its best tree improves. Only is used when max_restarts > 1.
    seed: int
        The random seed to use
    hill_climbing_only: bool
        Flag to only run hill climbing. (Default = False)
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
    region_reads_fn = os.path.join(temp_path, "region_reads.csv")
    mutation_location_fn = os.path.join(temp_path, "mutation_location.txt")
    region_weights_fn = os.path.join(temp_path, "region_weights.txt")
    cell_names_fn = os.path.join(temp_path, "cell_names.txt")
    gene_names_fn = os.path.join(temp_path, "gene_names.txt")
    clusters_fn = os.path.join(temp_path, "clusters.txt")
    mutation_clusters_fn = os.path.join(temp_path, "mutation_clusters.txt")

    output_path = os.path.join(temp_path)
    output_prefix = "out"
    
    clusters, pairs = [], []
    if constrained:
        if op.ul.DATA.CLUSTER_ID not in input_df.columns:
            raise Exception("%s isn't in the input data frame columns!" % op.ul.DATA.CLUSTER_ID) 
        clusters = input_df[op.ul.DATA.CLUSTER_ID].values
        input_df[op.ul.DATA.CLUSTER_ID].to_csv(clusters_fn, index=False, header=False)
        input_df = input_df.drop(columns=op.ul.DATA.CLUSTER_ID)
        clusterings, _ = op.tl.cluster.KModes(input_df.T, 
                                              list(input_df.index), 
                                              k=len(np.unique(clusters)))
        clusters = [list(input_df.columns[indices]) for indices in clusterings[0]]
        with open(mutation_clusters_fn, "w") as f:
            for cluster in clusters:
                f.write(" ".join(map(str, cluster)) + "\n")

    else:
        if op.ul.DATA.CLUSTER_ID in input_df.columns:
            input_df = input_df.drop(columns=op.ul.DATA.CLUSTER_ID)
        
    input_df.loc[:,input_df.columns != op.ul.DATA.CLUSTER_ID] = input_df.loc[:,input_df.columns != op.ul.DATA.CLUSTER_ID].replace(-1,3)
    input_df.T.to_csv(input_fn, sep=" ", index=False, header=False)
    input_df.index.to_series().to_csv(cell_names_fn, index=False, header=False)
    input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)

    if variant_reads_df is not None and total_reads_df is not None:
        assert variant_reads_df.shape == input_df.shape, "variant read count matrix and character matrix shape mismatch!"
        assert total_reads_df.shape == input_df.shape, "variant read count matrix and character matrix shape mismatch!"
        variant_reads_df.to_csv(variant_reads_fn)
        total_reads_df.to_csv(total_reads_fn)

    if region_reads_df is not None and meta_df is not None:
        assert np.array_equal(meta_df["NAME"], variant_reads_df.columns), "meta_df and variant_reads_df do not contain the same set of mutations"
        regions_index = [x.split("_")[1] for x in region_reads_df.index]
        region_index_series = pd.Series([regions_index.index(x) for x in meta_df["REGION"]])
        region_index_series.to_csv(mutation_location_fn, header=False, index=False)
        region_reads_df.to_csv(region_reads_fn, header=False)

        if region_weights_df is not None: 
            assert np.array_equal(regions_index, list(region_weights_df.index)), "region_reads_df and region_weights_df do not have the same index"
            region_weights_df.to_csv(region_weights_fn, index=False, header=False)

    n, m = input_df.shape
    
    # check minimum clone fractions
    assert len(min_clone_fractions) > 0 and all(isinstance(x, float) and 0 <= x <= 1 for x in min_clone_fractions), "min_clone_fractions should be a list of doubles (0,1)"


    args = [
            "-k", "%s" % k,
            "-e", "%d" % e,
            "-K", "%d" % K,
            "-R", "%d" % R,
            "-precision", "%.2f" % ado_precision,
            "-minCloneFractions", *[str(x) for x in min_clone_fractions],
            "-hcIterations", "%d" % hc_iterations,
            "-max_restarts", "%d" % max_restarts,
            "-patience", "%d" % patience,
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
            "-rr", "%s" % region_reads_fn if region_reads_df is not None else "",
            "-ri", "%s" % mutation_location_fn if meta_df is not None else "",
            "-rw", "%s" % region_weights_fn if region_weights_df is not None else "",
            "-seed" if seed is not None else "", "%d" % seed if seed is not None else "",
            "-greedy" if greedy else "",
            "-hcOnly" if hill_climbing_only else "",
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
                                       _type="cell_tree")
                                       
        # resolve genotypes for each cell
        clone_genotypes = np.array(T_cell.graph["genotypes"], dtype=int)
        cell_assignments = np.array(T_cell.graph["cell_assignments"], dtype=int)
        output_df = pd.DataFrame(clone_genotypes[cell_assignments], index=input_df.index, columns=input_df.columns)
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
