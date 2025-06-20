import os, shutil
import pandas as pd 
import numpy as np 
import networkx as nx

import oncophylo as op 

def BitSC2(input_df,
            alt_reads_df,
            total_reads_df,
            regions_df,
            meta_df,
            sex="female",
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
    ado_precision: int, optional
        The precision parameter for the beta distribution modeling allelic dropout. Default = 15
    remove_temp_dir: bool, optional
        Flag which when true removes the temporary directory that contains the output files. Default = True
        
    Returns
    -------
    dictionary
        A dictionary of results
    """
        
    temp_path = os.path.join(os.path.abspath(""),"bitsc2_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)

    basename = os.path.join(temp_path, "input")
    # variants_fn = os.path.join(temp_path, "input_variants.csv")
    # regions_fn = os.path.join(temp_path, "input_regions.csv")

    # # output file names
    # tree_fn = os.path.join(temp_path, "out_tree.gv")
    # genotypes_fn = os.path.join(temp_path, "out_nodes_genotypes.tsv")
    # cell_assignments_fn = os.path.join(temp_path, "out_cellAssignments.tsv")
    # copy_numbers_fn = os.path.join(temp_path, "out_nodes_copynumbers.tsv")

    n, m = input_df.shape
    SNV_to_region = meta_df["REGION"].values

    # Add regions which do not have any variants
    region2loci={}
    for x in regions_df.index:
        region2loci[x[x.find("_")+1:]] = []
    for i in range(meta_df.shape[0]):
        region2loci[SNV_to_region[i]].append(i)
    regions_depth=[]
    for x in regions_df.index:
        region = x[1+x.find("_"):]
        if len(region2loci[region])==0:
            regions_depth.append(regions_df.loc[x,:])
    if len(regions_depth)>0:
        regions_depth = np.array(regions_depth)
        DP_matrix = np.concatenate([total_reads_df.T.values,regions_depth],axis=0)
        AD_matrix = np.concatenate([alt_reads_df.T.values,np.zeros(regions_depth.shape)],axis=0)

    
    np.savetxt(basename+"_DP.csv",DP_matrix.astype(int),delimiter=",",fmt='%i')
    np.savetxt(basename+"_AD.csv",AD_matrix.astype(int),delimiter=",",fmt='%i')

    # Create genomic segments for BITSC2
    segments = []
    start_region = 0
    end_region = 0
    while end_region < m: 
        if end_region==m-1 or SNV_to_region[start_region] !=SNV_to_region[end_region+1]:
            segments.append((start_region+1,end_region+1))
            start_region = end_region+1
            end_region = start_region
        else:
            end_region+=1
    for i in range(m+1,AD_matrix.shape[0]+1):
        segments.append((i,i))
    np.savetxt(basename+"_segments.csv",np.array(segments).astype(int),delimiter=",",fmt='%i')


    # # write input files
    # variants_df = (total_reads_df - alt_reads_df).astype(str) + ":" + alt_reads_df.astype(str) + ":" + input_df.astype(str)
    # variants_df = pd.concat([meta_df.reset_index(), variants_df.T.reset_index(drop=True)], axis=1)
    # variants_df.to_csv(variants_fn, index=False, header=True)
    # regions_df.to_csv(regions_fn, index=True, header=False)

    # args = ["-i", "%s" % os.path.join(temp_path, "input"),
    #         "-o", "%s" % os.path.join(temp_path, "out"),
    #         "--nchains", "%d" % n_chains,
    #         "--chainlength", "%d" % chain_length,
    #         "--CNA", "%d" % int(infer_cnvs),
    #         "--sex", "%s" % sex,
    #         "-d", "%d" % int(infer_doublets),
    #         "--doubletrate", "%.6f" % double_rate,
    #         "--dropoutrate", "%.6f" % dropout_rate_concentration,
    #         "--seqerror", "%.6f" % seq_error,
    #         "--nodecost", "%d" % node_cost,
    #         "--cnacost", "%d" % cna_cost,
    #         "--lohcost", "%d" % loh_cost
    # ]

    # # run COMPASS
    # output, time = op.ul.subprocess([op.ul.binary_path("COMPASS")] + args)

    # T_cell, T_mut, output_genotypes = postprocess_COMPASS(tree_fn, cell_assignments_fn, genotypes_fn, copy_numbers_fn)
    # output_df = pd.DataFrame(output_genotypes, index=input_df.index, columns=input_df.columns)

    # solution = op.ul.solution(T_cell, 
    #                           T_mut, 
    #                           input_df,
    #                           output_df, 
    #                           seq_error,
    #                           dropout_rate,
    #                           output,
    #                           time,
    #                           T_clonal=nx.DiGraph(), # empty clonal tree 
    #                           var_reads=alt_reads_df,
    #                           total_reads=total_reads_df,
    #                           ado_precision=ado_precision)
    
    # # return solution
    # return solution

        
    # return T_cell, T_mut, genotypes


