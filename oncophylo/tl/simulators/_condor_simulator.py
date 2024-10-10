import itertools, math
import pandas as pd
import numpy as np
from scipy.stats import betabinom
import networkx as nx
import anndata as ad

import oncophylo as op 

def simulate(num_cells = 25,
             num_mutations = 25,
             num_clusters = 1,
             max_losses = 0,
             seed = 0,
             missing_rate = 0.0,
             mean_coverage = 50,
             fp_rate = 0.001,
             fn_rate = 0.001,
             ado_precision = 15,
             max_cn = 8,
             vaf_threshold = 0.1,
             mutation_rate = 0.2,
             prefix = "",
             save_path = ""):
    """This is a modified version of ConDoR's simulation framework 
    (https://github.com/raphael-group/ConDoR/blob/main/src/simulation_reads.py).
    This simulator will generate trees with losses that occur due to copy number loss. 
    
    Parameters
    -----------
    num_cells: int
        The number of cells in the simulation
    num_mutations: int
        The number of mutations in the simulation
    num_clusters: int
        The number of clusters in the simulation
    n_loss_per_character: int
        The number of losses per mutation in the simulation
    seed: int
        The random seed
    missing_rate: float
        The rate of missing mutations due to sequencing errors (introduces ? or 3 in mutation matrix)
    mean_coverage: int
        The mean coverage of sequencing
    fp_rate: float
        The false positive rate during sequencing
    fn_rate: float
        The false negative rate during sequencing
    ado_precision: int
        The precision parameter for the Beta distribution used to model allelic dropout
    max_cn: int 
        Maximum number of copies of genomic section
    vaf_threshold: float 
        The threshold of variant allele frequency (VAF) to call a mutation present in a cell
    mutation_rate: float
        The rate at which a locus is mutated
    prefix: string
        The prefix to give the simulation data when saving the file
    save_path: string
        The path to save the simulation data to
    
    Returns
    --------
    dictionary
        A dictionary containing the requested simulation data. The keys for this dictionary are 
        as follows: 
            'obs_data' = observed binary genotypes of each cell, last column is the copy cluster ID
            'true_data' = true cell genotypes, last column is the copy cluster ID
            'cell_tree' = true mutation tree with cell attachments used to generate the data
            'mutation_tree' = true mutation tree used to generate the data
            'cn_tree' = the true copy number tree
            'noise_free_read_counts' = dictionary containing noise free read counts in two dataframes 
                                       accessible with the 'variant' and 'total' keys
            'noisy_read_counts' = dictionary containing noisy read counts in two dataframes
                                  accessible with the 'variant' and 'total' keys
                                  
        
    """
    
    assert fp_rate >= 0 and fp_rate <= 1, "False positive rate must be in the range (0,1)"
    assert fn_rate >= 0 and fn_rate <= 1, "False negative rate must be in the range (0,1)"
    assert vaf_threshold >= 0 and vaf_threshold <= 1, "VAF threshold must be in the range (0,1)"
    assert mutation_rate >= 0 and mutation_rate <= 1, "Mutation rate must be in the range (0,1)"
    assert missing_rate >= 0 and missing_rate <= 1, "Missing rate must be in the range (0,1)"

    np.random.seed(seed)

    T = nx.DiGraph() # mutation tree
    Tc = nx.DiGraph() # copy number tree

    # add root nodes
    T.add_node('root')
    Tc.add_node(0)
    Tc.graph["root_id"] = 0
    Tc.graph["root_name"] = "root"
    
    # set types
    T.graph["type"] = op.ul.DATA.MUTATION_TREE
    Tc.graph["type"] = op.ul.DATA.CN_TREE
    T.graph["root_id"] = 0
    T.graph["root_name"] = "root"
    T.graph["loss_prefix"] = "-"
    T.graph["gain_prefix"] = "+"
    
    character_list = [f'm{character_index}' for character_index in range(num_mutations)]
    T.graph["mutations"] = character_list
    cluster_list = [f'{T.graph["loss_prefix"]}m{cluster_index}' for cluster_index in range(1, num_clusters)]
    event_order = np.random.permutation(character_list + cluster_list).tolist()

    loss_counter = np.zeros((num_mutations, 1))
    loss_dictionary = {f'{T.graph["loss_prefix"]}m{cluster_index}': [] for cluster_index in range(1, num_clusters)}
    losses = []
    
    # ground truth 
    B = np.zeros((num_mutations + num_clusters, num_mutations + 1), dtype=int)
    
    # number of copies of each allele in each cell
    R = np.zeros((num_clusters, num_mutations), dtype=int)    
    R[0, :] = np.random.randint(max_cn - max_losses - 1, size = num_mutations) + max_losses + 1
        
    # add losses 
    for node_index, event in enumerate(event_order):
        node_index += 1

        if event.startswith(f'{T.graph["loss_prefix"]}m'):
            parent_node_index = np.random.randint(node_index)
            parent_node = list(T.nodes)[parent_node_index]
            while parent_node.startswith(f'{T.graph["loss_prefix"]}m') or parent_node == "root":
                parent_node_index = np.random.randint(node_index)
                parent_node = list(T.nodes)[parent_node_index]
        else:
            parent_node_index = np.random.randint(node_index)
            parent_node = list(T.nodes)[parent_node_index]
                        
        # DO NOT ADD A MUTATION LOSS THAT'S ATTACHED TO THE ROOT
        # But, continue to do house keeping to complete the simulation 
        if event.startswith(f'{T.graph["loss_prefix"]}m'):
            if parent_node == "root":
                pass
            else:
                T.add_edge(parent_node, event)
        else:
            T.add_edge(parent_node, event)
        B[node_index, :] = B[parent_node_index, :]
        if event.startswith(f'{T.graph["loss_prefix"]}m'):
            cluster_id = int(event.lstrip(f'{T.graph["loss_prefix"]}m'))
            B[node_index, -1] = cluster_id
            parent_cluster_id = B[parent_node_index, -1]
            Tc.add_edge(parent_cluster_id, cluster_id)
            R[cluster_id, :] = R[parent_cluster_id, :]

            for mutation in range(num_mutations):
                if B[parent_node_index, mutation] == 1 and loss_counter[mutation] < max_losses:
                    if np.random.rand() < (1 - mutation_rate):
                        B[node_index, mutation] = loss_counter[mutation] + 2
                        loss_counter[mutation] += 1
                        loss_dictionary[event].append(mutation)
                        R[cluster_id, mutation] -= 1
                        losses.append(event)
        elif event.startswith('m'):
            mutation = int(event.lstrip('m'))
            B[node_index, mutation] = 1

    # randomize the number of copies of each allele for the mutations not impacted by copy-losses
    for mutation in range(num_mutations):
        if loss_counter[mutation] == 0:
            for cluster_id in range(num_clusters):
                R[cluster_id, mutation] = np.random.randint(max_cn - 1) + 1

    # check that all copy number states are non-zero positive
    assert(len(np.where(R == 0)[0]) == 0), "Mutations have a copy number of 0!"
    
    # check all SNV losses are supported by CNVs
    for cn_edge in Tc.edges:
        for mutation in loss_dictionary[f'{T.graph["loss_prefix"]}m{cn_edge[1]}']:
            assert(R[cn_edge[0], mutation] > R[cn_edge[1], mutation])    

    # assign cells and generate character-state matrix
    leaf_indices = []
    for idx, node in enumerate(T.nodes):
        if len(T[node]) == 0:
            leaf_indices.append(idx)    
    nleaves = len(leaf_indices)
    
    assert(num_cells > num_clusters)
    cell_assignment = np.random.randint(num_mutations, size=num_cells-nleaves)
    complete_cell_assignment = list(cell_assignment) + leaf_indices
    Bcell = B[complete_cell_assignment, :]
        
    # observed matrix
    A = B.copy()
    for mutation in range(num_mutations):
        A[A[:,mutation] > 1, mutation] = 0
    Acell = A[complete_cell_assignment, :]
    
    # cell tree
    celltree = T.copy()
    celltree.graph["type"] = op.ul.DATA.CELL_TREE
    for cell_id, assigned_node_index in enumerate(complete_cell_assignment):
        celltree.add_edge(list(T.nodes)[assigned_node_index], f's{cell_id}')

    # generate read counts
    Rtotal = np.zeros((num_cells, num_mutations), dtype=int)
    Vcount = np.zeros((num_cells, num_mutations), dtype=int)
    
    for cell in range(num_cells):
        for mutation in range(num_mutations):
            cluster_id = Acell[cell, -1]
            nvariant = Acell[cell, mutation]
            ntotal = R[cluster_id, mutation]

            latent_vaf = nvariant / ntotal
            
            nreads = np.random.poisson(mean_coverage)
            Rtotal[cell, mutation] = int(nreads)

            post_error_vaf = fp_rate + (1 - fp_rate - fn_rate) * latent_vaf
            ado_alpha = post_error_vaf * ado_precision
            ado_beta = ado_precision * (1 - post_error_vaf)
            nvariant_reads = betabinom.rvs(nreads, ado_alpha, ado_beta)
            Vcount[cell, mutation] = int(nvariant_reads)
            
    # generate the binarized mutation matrix
    VAF_mat = Vcount / Rtotal
    mutation_mat = (VAF_mat >= vaf_threshold).astype(int)
    mutation_mat = np.hstack((mutation_mat, Acell[:,-1][:,np.newaxis]))
    
    # introduce missing entries
    Acell_missing = Acell.copy()
    Rtotal_missing = Rtotal.copy()
    Vcount_missing = Vcount.copy()
    Acell_noisy = mutation_mat.copy()

    n_entries = num_cells * num_mutations
    nmissing = math.floor(missing_rate * n_entries)
    selected_cell_indices = np.random.randint(num_cells, size=nmissing)
    selected_character_indices = np.random.randint(num_mutations, size=nmissing)
    Acell_missing[selected_cell_indices, selected_character_indices] = -1
    Rtotal_missing[selected_cell_indices, selected_character_indices] = 0
    Vcount_missing[selected_cell_indices, selected_character_indices] = 0
    Acell_noisy[selected_cell_indices, selected_character_indices] = -1
    
    df_B = pd.DataFrame(B, index=["root"] + event_order,
                        columns = [f'm{idx}' for idx in range(num_mutations)] + [op.ul.DATA.CLUSTER_ID], dtype=int)            
    df_Bcell = pd.DataFrame(Bcell, index=[f's{idx}' for idx in range(num_cells)],
                            columns = [f'm{idx}' for idx in range(num_mutations)] + [op.ul.DATA.CLUSTER_ID], dtype=int)            
    df_Acell = pd.DataFrame(Acell, index=[f's{idx}' for idx in range(num_cells)],
                            columns = [f'm{idx}' for idx in range(num_mutations)] + [op.ul.DATA.CLUSTER_ID], dtype=int)    
    df_Acell_noisy = pd.DataFrame(Acell_noisy, index=[f's{idx}' for idx in range(num_cells)],
                                  columns = [f'm{idx}' for idx in range(num_mutations)] + [op.ul.DATA.CLUSTER_ID], dtype=int)
    
    df_Rtotal = pd.DataFrame(Rtotal, index=[f's{idx}' for idx in range(num_cells)],
                        columns = [f'm{idx}' for idx in range(num_mutations)], dtype=int)
    df_Vcount = pd.DataFrame(Vcount, index=[f's{idx}' for idx in range(num_cells)],
                        columns = [f'm{idx}' for idx in range(num_mutations)], dtype=int)    
    df_Rtotal_missing = pd.DataFrame(Rtotal_missing, index=[f's{idx}' for idx in range(num_cells)],
                        columns = [f'm{idx}' for idx in range(num_mutations)], dtype=int)
    df_Vcount_missing = pd.DataFrame(Vcount_missing, index=[f's{idx}' for idx in range(num_cells)],
                        columns = [f'm{idx}' for idx in range(num_mutations)], dtype=int)    
        
    celltree.graph["cells"] = list(df_Acell.index.values)
    celltree.graph["mutations"] = list(df_Acell.columns.values[:-1]) # everything except the cluster_id column
    celltree.graph["losses"] = losses

    adata = ad.AnnData(df_Acell_noisy.drop(columns=op.ul.DATA.CLUSTER_ID))
    adata.var["mutation_type"] = ["SNV"] * len(df_Acell_noisy.columns[:-1])  

    # collect all trees
    adata.uns[op.ul.DATA.CELL_TREE] = celltree
    adata.uns[op.ul.DATA.MUTATION_TREE] = T
    adata.uns[op.ul.DATA.CN_TREE] = Tc
    adata.uns[op.ul.DATA.CLONAL_TREE] = op.ul.to_clonal_tree(celltree, df_Acell.drop(columns=op.ul.DATA.CLUSTER_ID))
    
    # collect all cell/mutation data 
    adata.layers[op.ul.DATA.TRUE_DATA] = df_Acell.drop(columns=op.ul.DATA.CLUSTER_ID)
    adata.layers[op.ul.DATA.OBS_DATA] = df_Acell_noisy.drop(columns=op.ul.DATA.CLUSTER_ID)
    adata.layers[op.ul.DATA.VARIANT_READS] = df_Vcount
    adata.layers[op.ul.DATA.TOTAL_READS] = df_Rtotal
    adata.layers[op.ul.DATA.VARIANT_READS_CORRUPT] = df_Vcount_missing
    adata.layers[op.ul.DATA.TOTAL_READS_CORRUPT] = df_Rtotal_missing
    
    # collect cell specific data
    adata.obs[op.ul.DATA.CLUSTER_ID] = df_Acell[op.ul.DATA.CLUSTER_ID]
    
    # compute FPR and FNR and missing rate
    obs_values = df_Acell_noisy.drop(columns=op.ul.DATA.CLUSTER_ID).values
    true_values = df_Acell.drop(columns=op.ul.DATA.CLUSTER_ID).values
    adata.uns[op.ul.DATA.MUTANT_COPY_NUMBERS] = None
    adata.uns[op.ul.DATA.TOTAL_COPY_NUMBERS] = None
    adata.uns[op.ul.SIM_KEYS.FPR] = np.maximum(1e-6, np.sum((obs_values == 1) & (true_values == 0)) / np.sum(true_values == 0)) # percentage of 0's flipped to 1's
    adata.uns[op.ul.SIM_KEYS.FNR] = np.maximum(1e-6, np.sum((obs_values == 0) & (true_values == 1)) / np.sum(true_values == 1)) # percentage of 1's flipped to 0's
    adata.uns[op.ul.SIM_KEYS.MISSING_RATE] = np.sum(obs_values == -1) / obs_values.size # percentage of entries flipped to -1
        
    return adata