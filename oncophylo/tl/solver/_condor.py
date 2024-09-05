import os, subprocess, shutil
import pandas as pd 
import numpy as np 

import oncophylo as op 

def ConDoR(input_df, 
           alt_reads_df,
           total_reads_df,
           k=0, 
           fp=0.001, 
           fn=0.1, 
           ado=15,
           remove_temp_dir = True):
    
    """Python wrapper for ConDoR
    
    
    Parameters
    ----------
    input_df: pd.DataFrame
        A Pandas dataframe where rows are cells and columns are mutations and the entries at the present/absense of each SNV in each cell
    alt_reads_df: pd.DataFrame
        A Pandas dataframe where rows are cells and columns are mutations and the entries at the variant reads mapped per cell/SNV
    total_reads_df: pd.DataFrame
        A Pandas dataframe where rows are cells and columns are mutations and the entries at the total reads mapped per cell/SNV
    k: int
        The number of mutation losses to allow in the k-Dollo model
    fp: float, optional
        The false positive rate. Default = 0.001
    fn: float, optional
        The false negative rate. Default = 0.1 
    ado: int, optional
        The precision parameter for the beta distribution modeling allelic dropout. Default = 15
    remove_temp_dir: bool, optional
        Flag which when true removes the temporary directory that contains the output files. Default = True
        
    Returns
    -------
    dictionary
        A dictionary of results
    """
    # 2 is a loss -> 0, -1 is missing -> 3
    ConDoR_mutation_types = {2:0, -1:3}
    
    assert op.ul.CONST.CLUSTER_ID in input_df.columns, "Input must contain column: %s" % op.ul.CONST.CLUSTER_ID
    
    temp_path = os.path.join(os.path.abspath(""),"condor_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "character_matrix.csv")
    alt_reads_fn = os.path.join(temp_path, "alt_reads.csv")
    total_reads_fn = os.path.join(temp_path, "total_reads.csv")
    output_fn = os.path.join(temp_path, "out")

    n, m = input_df.shape
    
    input_df.loc[:,input_df.columns != op.ul.CONST.CLUSTER_ID] = input_df.loc[:,input_df.columns != op.ul.CONST.CLUSTER_ID].replace(3,-1)
    input_df.to_csv(input_fn, index=True, header=True)
    alt_reads_df.to_csv(alt_reads_fn, index=True, header=True)
    total_reads_df.to_csv(total_reads_fn, index=True, header=True)

    args = ["-i", "%s" % input_fn,
            "-v", "%s" % alt_reads_fn,
            "-r", "%s" % total_reads_fn,
            "-k", "%d" % k,
            "-a", "%.6f" % fp,
            "-b", "%.6f" % fn,
            "--ado", "%d" % ado,
            "-o", "%s" % output_fn
    ]

    # run SCITE
    output, time = op.ul.subprocess(["python", op.ul.script_path("ConDoR", "src/condor.py")] + args)
        
    mapping = {}
    mutations = list(input_df.columns) # don't really care if cluster_id column is in here
    for mut in mutations:
        for i in range(1,k+2):
            name = mut
            if i > 1:
                name = "-" + name
                if k > 1:
                    name = name + " (%d)" % i
                
            mapping[mut + "_%d" % i] = name
            
    # prepare output dataframe
    output_df = pd.read_csv(output_fn + "_B.csv", index_col=0, header=0)
    output_df.rename(columns=mapping, inplace=True)

    T_cell, T_mut = op.io.load_dot(output_fn + "_tree.dot", 
                                   mapping=mapping,
                                   mutations = list(input_df.columns),
                                   cells = list(input_df.index),
                                   _type="cell_tree",
                                   set_id_to_label = True)    

    solution  = op.ul.solution(T_cell, 
                               T_mut, 
                               input_df,
                               output_df, 
                               fp,
                               fn,
                               output,
                               time)
    
    return solution