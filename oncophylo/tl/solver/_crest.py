import os, subprocess, shutil
import pandas as pd 
import numpy as np 

import oncophylo as op 

def CREST(character_matrix, 
          variant_reads_df=None,
          total_reads_df=None,
          meta_df=None,
          fp=0.02, 
          fn=0.02, 
          iters=200,
          seed=None, 
          log_stderr=True,
          log_stdout=True,
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
    fp: float
        The estimated false positive rate for the sequencing data (SNVs)
    fn: float
        The estimated false negative rate for the sequencing data (SNVs)
    iters: int
        The maximum number of hill climbing iterations to perform for each subtree
    seed: int
        The random seed to use
    log_stderr: bool
        Flag to write the stderr to a text file
    log_stdout: bool
        Flag to write the stdout to a text file
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files (Default = True)
    destination_dir: str
        A directory to save all output files to

    Returns
    --------
    list or dictionary
        A dictionary or a list of dictionaries containing solutions and their corresponding meta data
    """
    
    temp_path = os.path.join(os.path.abspath(""),"crest_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    character_matrix_fn = os.path.join(temp_path, "character_matrix.csv")
    variant_reads_fn = os.path.join(temp_path, "variant_reads.csv")
    total_reads_fn = os.path.join(temp_path, "total_reads.csv")
    meta_fn = os.path.join(temp_path, "meta.csv")
    clusters_fn = os.path.join(temp_path, "clusters.txt")

    stderr_fn = os.path.join(temp_path, "stderr.txt")
    stdout_fn = os.path.join(temp_path, "stdout.txt")

    output_path = os.path.join(temp_path)
    output_prefix = "out"
        
    # process cluster_id
    if op.ul.DATA.CLUSTER_ID in character_matrix:
        clusters = character_matrix[op.ul.DATA.CLUSTER_ID]
        character_matrix = character_matrix.drop(op.ul.DATA.CLUSTER_ID, axis=1)
    else:
        clusters = pd.Series([0]*character_matrix.shape[0])

    # write character matrix and clusters to files
    character_matrix.to_csv(character_matrix_fn, index=True, header=True)
    clusters.to_csv(clusters_fn, header=False, index=False)

    # write variant and total reads to file
    if variant_reads_df is not None and total_reads_df is not None:
        assert variant_reads_df.shape == character_matrix.shape, "variant read count matrix and character matrix shape mismatch!"
        assert total_reads_df.shape == character_matrix.shape, "variant read count matrix and character matrix shape mismatch!"
        variant_reads_df.to_csv(variant_reads_fn, index=True, header=True)
        total_reads_df.to_csv(total_reads_fn, index=True, header=True)

    # write meta data to file
    if meta_df == None:
        meta_df = pd.DataFrame({"NAME": list(character_matrix.columns), "CLUSTER": np.arange(0, len(character_matrix.columns))})
    meta_df.to_csv(meta_fn, index=True, header=True)

    args_dict = {
        "-c": character_matrix_fn,
        "-v": variant_reads_fn if variant_reads_df is not None else "",
        "-t": total_reads_fn if total_reads_df is not None else "",
        "-m": meta_fn if meta_df is not None else "",
        "-s": clusters_fn,
        "-o": output_path,
        "-p": output_prefix,
        "-fp": f"{fp:.6f}",
        "-fn": f"{fn:.6f}",
        "-iters": str(iters),
    }

    # if seed is not None:
    #     args_dict["-seed"] = str(seed)

    # args = op.ul.convert_args(args_dict) # convert to list for subprocess

    # # run scOrchard
    # output, time = op.ul.subprocess([os.path.join(os.path.expanduser("~"), "crest/crest/bin/crest")] + args)
 
    # # log stderr
    # if log_stderr:
    #     with open(stderr_fn, "w") as f:
    #         f.write(output.stderr)

    # # log stdout
    # if log_stdout:
    #     with open(stdout_fn, "w") as f:
    #         f.write(output.stdout)

    # dot_fn = os.path.join(output_path, output_prefix + "_ml0.gv")
    # T_cell, T_mut = op.io.load_dot(dot_fn, 
    #                                _type="cell_tree")
                                    
    # # resolve genotypes for each cell
    # clone_genotypes = np.array(T_cell.graph["genotypes"], dtype=int)
    # cell_assignments = np.array(T_cell.graph["cell_assignments"], dtype=int)
    # corrected_character_matrix = pd.DataFrame(clone_genotypes[cell_assignments], index= character_matrix.index, columns=character_matrix.columns)

    # # save output files into a directory if provided a valid directory
    # op.ul.save_output_files(destination_dir, [dot_fn])

    # if remove_temp_dir:
    #     shutil.rmtree(temp_path)

    # return op.ul.solution(T_cell, 
    #                       T_mut, 
    #                       character_matrix,
    #                       corrected_character_matrix, 
    #                       output,
    #                       time)
