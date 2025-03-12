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
              hom_precision=15.0,
              het_precision=4.0,
              dropout_concentration=100.0,
              dropout_rate_prior=0.05,
              min_clone_fraction=0.02,
              theta=6.0,
              seed=None, 
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
    ado_precision: float
        The allelic dropout precision for the beta binomial likelihood model. This is only used if read counts are provided.
    min_clone_fraction: float
        The minimum fraction of cells that can be used to define a new clonal population. Default is 0.02
    seed: int
        The random seed to use
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
    meta_fn = os.path.join(temp_path, "meta.csv")
    variant_reads_fn = os.path.join(temp_path, "variant_reads.csv")
    total_reads_fn = os.path.join(temp_path, "total_reads.csv")
    region_reads_fn = os.path.join(temp_path, "region_reads.csv")
    mutation_location_fn = os.path.join(temp_path, "mutation_location.txt")
    region_weights_fn = os.path.join(temp_path, "region_weights.txt")
    snps_fn = os.path.join(temp_path, "snps.txt")
    cell_names_fn = os.path.join(temp_path, "cell_names.txt")
    gene_names_fn = os.path.join(temp_path, "gene_names.txt")
    clusters_fn = os.path.join(temp_path, "clusters.txt")
    mutation_clusters_fn = os.path.join(temp_path, "mutation_clusters.txt")

    output_path = os.path.join(temp_path)
    output_prefix = "out"
        
    input_df.loc[:,input_df.columns != op.ul.DATA.CLUSTER_ID] = input_df.loc[:,input_df.columns != op.ul.DATA.CLUSTER_ID].replace(-1,3)
    input_df.T.to_csv(input_fn, sep=" ", index=False, header=False)
    input_df.index.to_series().to_csv(cell_names_fn, index=False, header=False)
    # input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)

    if variant_reads_df is not None and total_reads_df is not None:
        assert variant_reads_df.shape == input_df.shape, "variant read count matrix and character matrix shape mismatch!"
        assert total_reads_df.shape == input_df.shape, "variant read count matrix and character matrix shape mismatch!"
        variant_reads_df.to_csv(variant_reads_fn)
        total_reads_df.to_csv(total_reads_fn)

    is_snp = [0]*meta_df.shape[0]
    meta = None
    if region_reads_df is not None and meta_df is not None:
        assert np.array_equal(meta_df["NAME"], variant_reads_df.columns), "meta_df and variant_reads_df do not contain the same set of mutations"
        meta = meta_df.copy()
        regions_ids = [x.split("_")[1] for x in region_reads_df.index]
        region_index = [regions_ids.index(x) for x in meta["REGION"]]
        meta["REGION_INDEX"] = region_index
        if "FREQ" in meta.columns:
            is_snp = [int(x > 0.0) for x in meta["FREQ"]]
        # snps_series.to_csv(snps_fn, header=False, index=False) # file with 0/1 indicating if mutation is SNP
        meta["SNP"] = is_snp
        #region_index_series.to_csv(mutation_location_fn, header=False, index=False) # file with each mutation location
        region_reads_df.to_csv(region_reads_fn, header=False) # file with region read counts

        if region_weights_df is not None: 
            assert np.array_equal(regions_index, list(region_weights_df.index)), "region_reads_df and region_weights_df do not have the same index"
            region_weights_df.to_csv(region_weights_fn, index=False, header=False)

        meta.to_csv(meta_fn, header=True)

    n, m = input_df.shape
    
    # check minimum clone fraction
    assert (min_clone_fraction >= 0.0) and (min_clone_fraction <= 1), "min_clone_fraction must be in range (0,1)"


    args = [
            "-k", "%s" % k,
            "-e", "%d" % e,
            "-homp", "%.2f" % hom_precision,
            "-hetp", "%.2f" % het_precision,
            "-dropoutc", "%.2f" % dropout_concentration,
            "-dropoutrp", "%.6f" % dropout_rate_prior,
            "-mincf", "%.6f" % min_clone_fraction,
            "-theta", "%.2f" % theta,
            "-n", "%d" % n,
            "-m", "%d" % m,
            "-i", "%s" % input_fn,
            "-meta", "%s" % meta_fn if meta is not None else "",
            "-o", "%s" % output_path,
            "-p", "%s" % output_prefix,
            "-c", "%s" % cell_names_fn,
            "-fd", "%.6f" % fp,
            "-ad", "%.6f" % fn,
            "-v", "%s" % variant_reads_fn if variant_reads_df is not None else "",
            "-t", "%s" % total_reads_fn if total_reads_df is not None else "",
            "-rr", "%s" % region_reads_fn if region_reads_df is not None else "",
            "-rw", "%s" % region_weights_fn if region_weights_df is not None else "",
            "-seed" if seed is not None else "", "%d" % seed if seed is not None else "",
    ]

    # run scOrchard
    output, time = op.ul.subprocess([os.path.join(os.path.expanduser("~"), "scaOrchard/bin/scaOrchard")] + args)

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
                                        total_reads=total_reads_df))
    if remove_temp_dir:
        shutil.rmtree(temp_path)
    if n_solutions == 1:
         return solutions[0]
    else:
        return solutions
