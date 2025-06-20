import os, subprocess, shutil
import pandas as pd 
import numpy as np 

import oncophylo as op 

def scaOrchard(character_matrix, 
               variant_reads_df=None,
               total_reads_df=None,
               region_reads_df=None,
               meta_df=None,
               cell_samples=None,
               fp=0.02, 
               fn=0.2, 
               hom_precision=15.0,
               het_precision=4.0,
               dropout_concentration=100.0,
               dropout_rate_prior=0.05,
               theta=6.0,
               mcmc_iters=200,
               seed=None, 
               remove_temp_dir=True,
               destination_dir=""):
    """
    scaOrchard Python Wrapper
    
    Parameters
    ----------
    character_matrix: pd.DataFrame
        A dataframe containing the called genotypes (0,1,2,-1) for each SNV/SNP in each cell
    variant_reads_df: pd.DataFrame
        A dataframe containing the count of variant reads for each SNV/SNP in each cell
    total_reads_df: pd.DataFrame
        A dataframe containing the count of total reads for each SNV/SNP in each cell
    region_reads_df: pd.DataFrame
        A dataframe containing the count of reads that fell into each genomic region in each cell
    meta_df: pd.DataFrame
        A dataframe containing meta data for each variant
    cell_samples: np.array
        A NumPy array containing an integer value indicating which sample each cell came from
    fp: float
        The estimated false positive rate for the sequencing data
    fn: float
        The estimated false negative rate for the sequencing data
    hom_precision: float
        The precision parameter for the beta binomial likelihood model when the loci is homozygous
    het_precision: float
        The precision parameter for the beta binomial likelihood model when the loci is heterozygous
    dropout_concentration: float
        The dropout concentration parameter used when updating the dropout rates for each loci
    dropout_rate_prior: float
        The prior dropout rate used when updating the dropout rates for each loci
    mcmc_iters: int
        The maximum number of Markov Chain Monte Carlo moves at each iteration for sampling CNAs
    seed: int
        The random seed to use
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files (Default = True)
    destination_dir: str
        A directory to save all output files to

    Returns
    --------
    list or dictionary
        A dictionary or a list of dictionaries containing solutions and their corresponding meta data
    """
    
    temp_path = os.path.join(os.path.abspath(""),"orchard_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    character_matrix_fn = os.path.join(temp_path, "character_matrix.csv")
    meta_fn = os.path.join(temp_path, "meta.csv")
    variant_reads_fn = os.path.join(temp_path, "variant_reads.csv")
    total_reads_fn = os.path.join(temp_path, "total_reads.csv")
    region_reads_fn = os.path.join(temp_path, "region_reads.csv")
    cell_samples_fn = os.path.join(temp_path, "cell_samples.txt")

    output_path = os.path.join(temp_path)
    output_prefix = "out"
        
    character_matrix.to_csv(character_matrix_fn, index=True, header=True)

    if variant_reads_df is not None and total_reads_df is not None:
        assert variant_reads_df.shape == character_matrix.shape, "variant read count matrix and character matrix shape mismatch!"
        assert total_reads_df.shape == character_matrix.shape, "variant read count matrix and character matrix shape mismatch!"
        variant_reads_df.to_csv(variant_reads_fn, index=True, header=True)
        total_reads_df.to_csv(total_reads_fn, index=True, header=True)

    is_snp = [0]*meta_df.shape[0]
    meta = None
    if region_reads_df is not None and meta_df is not None:
        assert np.all(np.sort(meta_df["NAME"]) == np.sort(variant_reads_df.columns)), "meta_df and variant_reads_df do not contain the same set of mutations"
        meta = meta_df.copy()
        regions_ids = [x.split("_")[1] for x in region_reads_df.index]
        region_index = [regions_ids.index(x) for x in meta["REGION"]]
        meta["REGION_INDEX"] = region_index
        if "FREQ" in meta.columns:
            is_snp = [int(x > 0.0) for x in meta["FREQ"]]

        meta["SNP"] = is_snp
        region_reads_df.T.to_csv(region_reads_fn, index=True, header=True) # file with region read counts

    # cell samples is a label for which sample each cell belongs to
    # if no cell_samples are provided and no sample label is provided to SNVs, then we assume
    # everything is from the same sample
    if cell_samples is None:
        cell_samples = pd.Series([0]*character_matrix.shape[0])
    if "SAMPLE" not in meta.columns:
        meta["SAMPLE"] = 0
    cell_samples.to_csv(cell_samples_fn, header=False, index=False)
    meta.to_csv(meta_fn, index=True, header=True)

    args = [
            "-homp", "%.2f" % hom_precision,
            "-hetp", "%.2f" % het_precision,
            "-dropoutc", "%.2f" % dropout_concentration,
            "-dropoutrp", "%.6f" % dropout_rate_prior,
            "-theta", "%.2f" % theta,
            "-c", "%s" % character_matrix_fn,
            "-meta", "%s" % meta_fn if meta is not None else "",
            "-o", "%s" % output_path,
            "-p", "%s" % output_prefix,
            "-s", "%s" % cell_samples_fn,
            "-fp", "%.6f" % fp,
            "-fn", "%.6f" % fn,
            "-v", "%s" % variant_reads_fn if variant_reads_df is not None else "",
            "-t", "%s" % total_reads_fn if total_reads_df is not None else "",
            "-r", "%s" % region_reads_fn if region_reads_df is not None else "",
            "-iters", "%d" % mcmc_iters, 
            "-seed" if seed is not None else "", "%d" % seed if seed is not None else "",
    ]

    # run scOrchard
    output, time = op.ul.subprocess([os.path.join(os.path.expanduser("~"), "scaOrchard/bin/scaOrchard")] + args)
 
    dot_fn = os.path.join(output_path, output_prefix + "_ml0.gv" % i)
    T_cell, T_mut = op.io.load_dot(dot_fn, 
                                    _type="cell_tree")
                                    
    # resolve genotypes for each cell
    clone_genotypes = np.array(T_cell.graph["genotypes"], dtype=int)
    cell_assignments = np.array(T_cell.graph["cell_assignments"], dtype=int)
    out_character_matrix = pd.DataFrame(clone_genotypes[cell_assignments], index= character_matrix.index, columns=character_matrix.columns)

    # save output files into a directory if provided a valid directory
    op.ul.save_output_files(destination_dir, [dot_fn])

    if remove_temp_dir:
        shutil.rmtree(temp_path)

    return op.ul.solution(T_cell, 
                          T_mut, 
                          character_matrix,
                          out_character_matrix, 
                          fp,
                          fn,
                          output,
                          time,
                          var_reads=variant_reads_df,
                          total_reads=total_reads_df)
