import os, shutil, copy
import pandas as pd 
import numpy as np 
import networkx as nx

import oncophylo as op 

def COMPASS(input_df,
            alt_reads_df,
            total_reads_df,
            regions_df,
            meta_df,
            n_chains=4,
            chain_length=10000,
            infer_cnvs=True,
            infer_doublets=False,
            double_rate=0.08,
            dropout_rate=0.05,
            dropout_rate_concentration=100,
            seq_error=0.02,
            node_cost=1,
            cna_cost=85,
            loh_cost=85,
            ado_precision=15.0,
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
        
    temp_path = os.path.join(os.path.abspath(""),"compass_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    variants_fn = os.path.join(temp_path, "input_variants.csv")
    regions_fn = os.path.join(temp_path, "input_regions.csv")

    # output file names
    tree_fn = os.path.join(temp_path, "out_tree.gv")
    genotypes_fn = os.path.join(temp_path, "out_nodes_genotypes.tsv")
    cell_assignments_fn = os.path.join(temp_path, "out_cellAssignments.tsv")
    copy_numbers_fn = os.path.join(temp_path, "out_nodes_copynumbers.tsv")

    # write input files
    variants_df = (total_reads_df - alt_reads_df).astype(str) + ":" + alt_reads_df.astype(str) + ":" + input_df.astype(str)
    variants_df = pd.concat([meta_df.reset_index(drop=True), variants_df.T.reset_index(drop=True)], axis=1)
    variants_df.to_csv(variants_fn, index=False, header=True)
    regions_df.to_csv(regions_fn, index=True, header=False)

    args = ["-i", "%s" % os.path.join(temp_path, "input"),
            "-o", "%s" % os.path.join(temp_path, "out"),
            # "--nchains", "%d" % n_chains,
            "--chainlength", "%d" % chain_length,
            # "--CNA", "%d" % int(infer_cnvs) if not infer_cnvs else "",
            # "--sex", "%s" % sex,
            "-d", "%d" % int(infer_doublets),
            # "--doubletrate", "%.6f" % double_rate,
            # "--dropoutrate", "%.6f" % dropout_rate_concentration,
            # "--seqerror", "%.6f" % seq_error,
            # "--nodecost", "%d" % node_cost,
            # "--cnacost", "%d" % cna_cost,
            # "--lohcost", "%d" % loh_cost
    ]

    # run COMPASS
    output, time = op.ul.subprocess([op.ul.binary_path("COMPASS")] + args)

    T_cell, T_mut, T, output_genotypes = postprocess_COMPASS(tree_fn, cell_assignments_fn, genotypes_fn, copy_numbers_fn)
    output_df = pd.DataFrame(output_genotypes, index=input_df.index, columns=input_df.columns)

    solution = op.ul.solution(T_cell, 
                              T_mut, 
                              input_df,
                              output_df, 
                              seq_error,
                              dropout_rate,
                              output,
                              time,
                              T_clonal=nx.DiGraph(), # empty clonal tree 
                              var_reads=alt_reads_df,
                              total_reads=total_reads_df,
                              ado_precision=ado_precision)

    solution["tree_output"] = T
    
    # return solution
    return solution

def postprocess_COMPASS(tree_fn, cell_assignments_fn, genotypes_fn, copy_numbers_fn):
    
    # read dot file
    T_mut = nx.nx_pydot.read_dot(tree_fn)

    cell_assignments_df = pd.read_csv(cell_assignments_fn, header=0, index_col=0, sep="\t")
    node_genotypes_df = pd.read_csv(genotypes_fn, index_col=0, sep="\t")
    copy_numbers_df = pd.DataFrame()
    if os.path.isfile(copy_numbers_fn): # only load if we ran with CNA = 1
        copy_numbers_df = pd.read_csv(copy_numbers_fn, header=0, index_col=0, sep="\t")

    # process labels so they match scOrchard
    for node in T_mut.nodes:
        if "label" in T_mut.nodes[node]:
            label_parts = []
            label = T_mut.nodes[node]["label"]
            label_split = label.split("<br/>")
            for l in label_split:
                l = l.replace("<B>","").replace("</B>","").replace("_", " ")
                l = l.split(":")
                if len(l) > 1:
                    split = l[1].split(" ")
                    allele = split[0]
                    chrom = split[1].replace("(","").replace(")","")
                    label_parts.append(f"{l[0]} ({chrom}, {allele})")
                else:
                    label_parts.append(l[0].replace("(", " (").replace("  ", " "))
            if len(label_parts) > 1:
                T_mut.nodes[node]["label"] = "<br/>".join(label_parts)
            else:
                T_mut.nodes[node]["label"] = label_parts[0]

    T = copy.deepcopy(T_mut.copy())

    # remove nodes from tree that aren't relevant
    node_indices = [n.split()[1] for n in node_genotypes_df.index]
    nodes_to_remove = np.setdiff1d(T_mut.nodes(), node_indices).tolist()
    T_mut.remove_nodes_from(nodes_to_remove)
    node_names = [f"Node {node}" for node in T_mut.nodes]
    cell_assignments = cell_assignments_df["node"].values

    # add tree graph
    T_mut.graph[op.ul.DATA.TOTAL_COPY_NUMBERS] = copy_numbers_df.T.values
    T_mut.graph["mutations"] = node_names
    T_mut.graph["root_id"] = str(0)
    T_mut.graph["root_name"] = str(0)
    T_mut.graph["cell_assignments"] = cell_assignments

    T_cell = T_mut.copy() # create cell tree

    # add cell attachments
    for i, row in cell_assignments_df.iterrows():
        T_cell.add_edge(str(row["node"]), i)

    T_cell.graph["cells"] = cell_assignments_df.index.tolist()
    T_cell.graph["type"] = op.ul.DATA.CELL_TREE


    # process genotypes
    genotypes = node_genotypes_df.values[cell_assignments]
        
    return T_cell, T_mut, T, genotypes