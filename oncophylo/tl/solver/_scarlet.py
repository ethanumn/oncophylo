# _scarlet.py
import os, subprocess, shutil
import pandas as pd 
import numpy as np 

import oncophylo as op 

def scarlet(input_df,
            alt_reads_df,
            total_reads_df,
            cn_tree,
            ado_precision=15,
            fp=0.001,
            fn=0.2,
            remove_temp_dir = True):
    
    """Python wrapper for SCARLET
    
    
    Parameters
    ----------
    input_df: pd.DataFrame
        A Pandas dataframe where rows are cells and columns are mutations and the entries at the present/absense of each SNV in each cell
    alt_reads_df: pd.DataFrame
        A Pandas dataframe where rows are cells and columns are mutations and the entries at the variant reads mapped per cell/SNV
    total_reads_df: pd.DataFrame
        A Pandas dataframe where rows are cells and columns are mutations and the entries at the total reads mapped per cell/SNV
    cn_tree: Networkx.DiGraph
        A copy number tree
    ado_precision: float, optional
        scarlet does not use this parameter, it's only used to compute the beta-binomial likelihood
    fp: float, optional
        The false positive rate. This is only for scoring the genotype matrix it is not used by SCARLET. Default = 0.001
    fn: float, optional
        The false negative rate. This is only for scoring the genotype matrix it is not used by SCARLET.Default = 0.1 
    remove_temp_dir: bool, optional
        Flag which when true removes the temporary directory that contains the output files. Default = True
        
    Returns
    -------
    dictionary
        A dictionary of results
    """
    
    temp_path = os.path.join(os.path.abspath(""),"scarlet_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    reads_fn = os.path.join(temp_path, "reads.csv")
    cn_tree_fn = os.path.join(temp_path, "tree.csv")
    output_fn = os.path.join(temp_path, "out")

    # make a dataframe containing alt reads, total reads, and copy number cluster assignment
    reads_df = pd.DataFrame(index=input_df.index)
    reads_df["c"] = input_df[op.ul.DATA.CLUSTER_ID]
    reads_df.index.name = "cell_id"
        
    # make reads.csv from alt and total reads df
    for mut in alt_reads_df.columns:
        if mut == op.ul.DATA.CLUSTER_ID:
            continue
        reads_df[mut + "_v"] = alt_reads_df[mut]
        reads_df[mut + "_t"] = total_reads_df[mut]

    cn_tree_lines = []

    # now prepare copy number tree file
    for e in cn_tree.edges():
        line = [e[0],e[1]]

        # get all cell ids in cluster
        cells_in_cluster = reads_df.index[reads_df["c"] == e[1]].values.tolist() 
        line += cells_in_cluster
        cn_tree_lines.append(line)
    
    # write reads and cn tree file 
    reads_df.to_csv(reads_fn, index=True, header=True)
    pd.DataFrame(cn_tree_lines).to_csv(cn_tree_fn, header=False, index=False)

    args = ["%s" % reads_fn,
            "%s" % cn_tree_fn,
            "%s" % output_fn,
            "ALL"
    ]

    # capture runtime for subprocess
    output, time = op.ul.subprocess(["bash", op.ul.script_path("scarlet", "code/scarlet.sh")] + args)
    
    # prepare output dataframe
    output_df = pd.read_csv(output_fn + ".T", index_col=0, header=0)
    output_df = output_df[input_df.columns[:-1]] # scarlet rearranges the columns, reverse this
    
    mapping = {"ROOT:0":"root"}
    
    # add roots for each copy clone
    for c in np.unique(reads_df["c"]):
        if c == 0:
            continue
        mapping["ROOT:%d"] = "cn_clone%d" % c
        
    mutations = list(input_df.columns)
    cells = list(input_df.index)
    for mut in mutations:   
        mapping["MUT:"+mut] = mut
    for cell in cells:
        mapping["CELL:"+cell] = cell
            
    T_cell, T_mut = op.io.load_dot(output_fn + ".dot", 
                                   mapping=mapping,
                                   mutations = list(input_df.columns),
                                   cells = list(input_df.index),
                                   _type="cell_tree")    

    solution = op.ul.solution(T_cell, 
                              T_mut, 
                              input_df.drop(columns=op.ul.DATA.CLUSTER_ID),
                              output_df, 
                              fp,
                              fn,
                              output,
                              time,
                              var_reads=alt_reads_df,
                              total_reads=total_reads_df,
                              ado_precision=ado_precision)
    
    return solution