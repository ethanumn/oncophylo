# _scorchard_simulator.py 

import os, sys
import pandas as pd
import numpy as np
import networkx as nx
import anndata as ad

import oncophylo as op 

# CONSTRAINTS
NONE = 0
LOSS_ONLY = 1
GAIN_ONLY = 2
LOSS_AND_GAIN = 3

def simulate_tree(num_mutations, 
                  num_cells, 
                  num_clusters, 
                  K, 
                  R,
                  constraint,
                  loss_prob=0.1,
                  gain_prob=0.05,
                  max_cn=2,
                  seed=None,
                  max_retries=3):
    """Simulates a mutation tree by generating a random undirected graph using Networkx,
    making it into a rooted directed tree, then add loss/gains. This only simulates
    trees that adhere to a limited version of the finite-sites model. It considers only
    the following per clone:
        mutation acquired 
        mutation acquired -> lost 
        
    A clone here means all ancestral states that can be reached from a leaf node in a mutation tree. Clones
    can share common ancestry, but each mutation is gained at most once and lost at most once in each clone."""
    import networkx as nx
    
    # here are a set of helper functions
    def process_event(e):
        """Identify the original mutation associated with this event"""
        if e.startswith("+"):
            return e.lstrip("+").split(" ")[0]
        elif e.startswith("-"):
            return e.lstrip("-").split(" ")[0]
        else:
            return e
        
    def check_events(gen, chosen_losses, chosen_gains, possible_losses, possible_gains, p_loss=0.5):
        """Checks and makes sure at least one event (gain or loss) has occurred"""
        if len(chosen_losses) == 0 and len(chosen_gains) == 0:
            if len(possible_losses) > 0 and len(possible_gains) > 0:
                if gen.binomial(1,p=0.5):
                    chosen_losses.append(gen.choice(possible_losses))
                else:
                    chosen_gains.append(gen.choice(possible_gains))
            elif len(possible_losses) > 0:
                chosen_losses = np.array([gen.choice(possible_losses)])
            elif len(possible_gains) > 0:
                chosen_gains = np.array([gen.choice(possible_gains)])
            else:
                print("No gains or losses possible, retrying...")
        return chosen_losses, chosen_gains
    
    def find_losses(nodes):
        """Return anything will a loss prefix"""
        losses = []
        for n in nodes:
            if n.startswith("-"):
                losses.append(n.lstrip("-").split(" ")[0])
        return set(losses)
    
    def find_gains(nodes):
        """Return anything will a gain prefix"""
        gains = []
        for n in nodes:
            if n.startswith("+"):
                gains.append(n.lstrip("+").split(" ")[0])
        return set(gains)
    
    def process_events(chosen_losses, 
                       chosen_gains, 
                       loss_count, 
                       gain_count, 
                       losses, 
                       gains,
                       K,
                       R):
        """Processes the events so that the proper nodes are added to the tree"""
        events_to_add = []
        # parse losses
        for e in chosen_losses:
            mut = process_event(e)
            if loss_count[mut] >= K:
                continue
            else:
                loss_count[mut] += 1
                if loss_count[mut] > 1:
                    new_loss = f"-{mut} ({loss_count[mut]})"
                else:
                    new_loss = f"-{mut}"
                losses.append(new_loss)
                events_to_add.append(new_loss)

        # parse gain
        for e in chosen_gains:
            mut = process_event(e)
            if gain_count[mut] >= R:
                continue
            else:
                gain_count[mut] += 1
                if gain_count[mut] > 1:
                    new_gain = f"+{mut} ({gain_count[mut]})"
                else:
                    new_gain = f"+{mut}"
                gains.append(new_gain)
                events_to_add.append(new_gain)
                
        return events_to_add
    
    assert num_clusters <= num_mutations, "The number of clusters %d must be greater than the number of mutations %d" % (num_clusters, num_mutations)
    assert max_cn > 0, "The maximum copy number must be non-zero, %d was provided" % max_cn
    ##################################
    # Step 1: simulate mutation tree #
    ##################################
    gen = np.random.default_rng(seed)

    mutations = [f"m{i}" for i in range(num_mutations)]
    cells = [f"s{i}" for i in range(num_cells)]

    # Step 1: Generate a random undirected tree
    undirected_tree = nx.random_tree(num_mutations + 1)

    # Step 2: Choose a root node
    root = gen.choice(list(undirected_tree.nodes))

    # Step 3: Convert to a directed tree using depth-first search (DFS) and relabel nodes
    mutation_tree = nx.DiGraph()
    mutation_tree.add_edges_from(nx.dfs_edges(undirected_tree, source=root))
    
    mapping = {mut:mutations[i] for i,mut in enumerate(np.setdiff1d(mutation_tree.nodes, [root]))}
    mapping[root] = "root"
    
    mutation_tree = nx.relabel_nodes(mutation_tree, mapping=mapping)
    mutation_tree.graph["root_name"] = "root"
    mutation_tree.graph["cells"] = cells
    mutation_tree.graph["mutations"] = mutations
    mutation_tree.graph["loss_prefix"] = "-"
    mutation_tree.graph["gain_prefix"] = "+"
    mutation_tree.graph["type"] = op.ul.DATA.CELL_TREE

    #########################################################
    # Step 2: add in gains/losses that adhere to clustering #
    #########################################################
    gains = []
    losses = []
    loss_count = {mut:0 for mut in mutations}
    gain_count = {mut:0 for mut in mutations}
    subtree_roots = ["root"]
    initiating_events = []
        
    # create clustering
    if num_clusters > 1 and (K > 0 or R > 0) and (constraint is not NONE):
        valid_clusters = 1
        retries = 0
        while valid_clusters < num_clusters:

            possible_gains, possible_losses, chosen_gains, chosen_losses = [], [], [], []
            break_point = gen.choice(mutations)
            ancestors = set(nx.ancestors(mutation_tree, break_point))
            descendants = set(nx.descendants(mutation_tree, break_point))

            if K > 0 and (constraint == LOSS_ONLY or constraint == LOSS_AND_GAIN):                
                # mutations acquired above breakpoint can be lost, but only those that haven't been lost on this lineage
                losses_in_lineage = find_losses(ancestors.union(descendants))
                gains_above = find_gains(ancestors)
                gains_above.add(break_point)
                possible_losses = ((ancestors.union(gains_above)).intersection(set(mutations))).difference(losses_in_lineage)
                possible_losses = np.array(list(possible_losses))
                chosen_losses = possible_losses[[bool(gen.binomial(1, p=loss_prob)) for i in range(len(possible_losses))]].tolist()

            if R > 0 and (constraint == GAIN_ONLY or constraint == LOSS_AND_GAIN):
                # gains can include any mutation that's hasn't occurred on this lineage, or has been lossed above the breakpoint
                gains_in_lineage = find_gains(ancestors.union(descendants))
                gains_in_lineage.add(break_point)
                possible_gains = np.array(list(set(mutations).difference(ancestors.union(descendants).union(gains_in_lineage))))
                chosen_gains = possible_gains[[bool(gen.binomial(1, p=gain_prob)) for i in range(len(possible_gains))]].tolist()
            
            chosen_losses, chosen_gains = check_events(gen, chosen_losses, chosen_gains, possible_losses, possible_gains, p_loss=loss_prob/(loss_prob+gain_prob))

            # if we can't create the desired cluster then try a new break point
            if len(chosen_losses) == 0 and len(chosen_gains) == 0:
                retries += 1
                if retries > max_retries:
                    print("Unable to create phylogeny, please check your parameters!")
                    return
                continue
            else:
                retries = 0
                valid_clusters += 1 # update if we have a new cluster
            
            events_to_add = process_events(chosen_losses, 
                                           chosen_gains, 
                                           loss_count, 
                                           gain_count, 
                                           losses, 
                                           gains,
                                           K,
                                           R)
            
            if len(events_to_add) == 0:
                retries += 1
                continue
            
            # randomly choose zero or more children to parent from the break point
            children = list(mutation_tree.successors(break_point))

            size = gen.integers(0, len(children) + 1)
            selected_children = gen.choice(children, size=size, replace=False)
            
            # remove edges from break_point -> selected_children
            for ch in selected_children:
                mutation_tree.remove_edge(break_point, ch)

            # add all new edges
            p = break_point
            for e in events_to_add:
                mutation_tree.add_edge(p, e)
                p = e 
            
            # add edge from last event -> selected_children
            for ch in selected_children:
                mutation_tree.add_edge(p, ch)
            subtree_roots.append(e) # last event is the subtree root for this subclone 
            initiating_events.append(events_to_add[0])

    mutation_tree.graph["losses"] = losses
    mutation_tree.graph["gains"] = gains
    
    ####################################
    # Step 3: resolve copy number tree #
    ####################################
    copy_number_tree = nx.DiGraph()
    
    # resolve ancestral relationships between the subtree roots, and this is the copy number tree
    for i,r in enumerate(subtree_roots):
        
        if r == "root":
            copy_number_tree.add_node(0)
            continue
        else:   
            ancestors = nx.ancestors(mutation_tree, r)
            relevant_ancestors = set(subtree_roots) & ancestors

            closest_ancestor = None
            min_distance = float('inf')

            for ancestor in relevant_ancestors:
                distance = nx.shortest_path_length(mutation_tree, ancestor, r)
                if distance < min_distance:
                    min_distance = distance
                    closest_ancestor = ancestor
            copy_number_tree.add_edge(subtree_roots.index(closest_ancestor),i)

        
    #####################
    # Step 4: add cells #
    #####################
    cell_attachments = []
    attachment_points = set()
    clusters = []

    # attach at least one cell to a mutation for each clone, and each leaf should have a cell attached
    for r in subtree_roots:
        descendants = nx.descendants(mutation_tree, r)
        cluster_mutations = set(descendants.copy())
        for d in descendants:
            if d in initiating_events:
                
                cluster_mutations = cluster_mutations.difference(set(nx.descendants(mutation_tree, d)).union({d}))
        cluster_mutations.add(r)
        clusters.append(list(cluster_mutations))
        cell_attachments.append(gen.choice(list(cluster_mutations)))
        attachment_points = attachment_points.union(cluster_mutations)
        
    cell_attachments += [node for node in mutation_tree.nodes if mutation_tree.out_degree(node) == 0]
    
    if len(cell_attachments) > num_cells:
        print("Not enough cells to make it so the entire tree can be resolve from the character matrix. Please \
        increase the number of cells or decrease the number of mutations.")
        return

    attachment_points = list(attachment_points)
    # sample the remaining cell attachment points
    cell_attachments += gen.choice(attachment_points, size=num_cells-len(cell_attachments), replace=True).tolist()
    gen.shuffle(cell_attachments) # shuffle

    # resolve cluster assignment and copy states for each cell
    cluster_assignments = []
    cell_tree = mutation_tree.copy()    
    # now add cells to tree 
    for i,c in enumerate(cells):
        cell_tree.add_edge(cell_attachments[i], c)  
        
        # record which cluster cell c is in 
        for j in range(len(clusters)):
            if cell_attachments[i] in clusters[j]:
                cluster_assignments.append(j)
                break

    # resolve the character matrix
    copy_states = pd.DataFrame(0.0, index=cells, columns=mutations)
    _, character_matrix = op.ul.resolve_genotypes(cell_tree.copy())
    character_matrix = character_matrix.replace(op.ul.mutation_types)
        
    # simulate copy number states for each clone
    total_copy_numbers = np.zeros((num_clusters, num_mutations), dtype=int)
    mutant_copy_numbers = np.zeros((num_clusters, num_mutations), dtype=int)
    for c in range(num_clusters):
        missing_in_clone = np.flatnonzero(np.sum(character_matrix.replace(op.ul.mutation_types).values[np.flatnonzero(np.array(cluster_assignments) == c)], axis=0) == 0)
        for j in range(num_mutations):
            # get all of the cells in cluster, and fin
            np.flatnonzero(np.array(cluster_assignments) == c)
            total_copy_numbers[c,j] = gen.integers(low=1, high=max_cn+1) 
            if j in missing_in_clone:
                mutant_copy_numbers[c,j] = 0
            else:
                mutant_copy_numbers[c,j] = gen.integers(low=1, high=total_copy_numbers[c,j]+1)

    # compute copy states using mutant and total copy numbers
    for i,cell in enumerate(character_matrix.index):
        for j,mut in enumerate(character_matrix.columns):
            if character_matrix.loc[cell,mut] == 1:
                c = cluster_assignments[i]
                copy_states.loc[cell,mut] = mutant_copy_numbers[c,j] / total_copy_numbers[c,j]

    cell_tree.graph["cluster_assignments"] = cluster_assignments
    mutation_tree.graph["cluster_assignments"] = cluster_assignments

    
    return cell_tree, mutation_tree, copy_number_tree, copy_states, mutant_copy_numbers, total_copy_numbers


from scipy.stats import betabinom

def simulate_reads(character_matrix, 
                   mean_coverage=50,
                   ado_precision=5,
                   fp_rate=0.001,
                   fn_rate=0.001,
                   missing_rate=0.0,
                   vaf_threshold=0.1,
                   unknown_value=3,
                   copy_states="diploid",
                   seed=None):
    """Simulates observed and ground truth variant and total reads for each cells
    
    Parameters
    -----------
    character_matrix: pd.DataFrame
        A binary character matrix where rows and cells and columns are mutations
    mean_coverage: int
        The average number of reads mapped to each locus
    ado_precision: int
        The precision parameter for the Beta distribution used to model allelic dropout
    fp_rate: float, optional
        The false positive rate used to corrupt the read counts and observed character matrix.
    fn_rate: float, optional
        The false negative rate used to corrupt the read counts and observed character matrix.
    missing_rate: float, optional
        The rate at which each entry in the observed character matrix will be missing.
    vaf_threshold: float, optional
        The threshold of variant allele frequency (VAF) to call a mutation present in a cell
    copy_states: str or pd.DataFrame, optional
        If a string, then it must be either 'diploid' or 'haploid'. If a Pandas dataframe, then the rows and columns
        should be the same as the character matrix, but each entry is the fraction of copies at the locus 
        containing the mutation that are mutant alleles. This is used as an expected variant allele frequency.
        
    Returns
    -------
    pd.DataFrame
        An 'observed' character matrix like one might obtain after calling variants from a sequencing experiment
    """
    assert np.all(np.unique(character_matrix.values) == np.array([0,1])), "Non-binary characters detected in character matrix, please ensure all characters are [0,1]"
    assert fp_rate >= 0 and fp_rate <= 1, "False positive rate must be in the range (0,1)"
    assert fn_rate >= 0 and fn_rate <= 1, "False negative rate must be in the range (0,1)"
    assert vaf_threshold >= 0 and vaf_threshold <= 1, "VAF threshold must be in the range (0,1)"
    assert missing_rate >= 0 and missing_rate <= 1, "Missing rate must be in the range (0,1)"
    
    # get the number of cells and mutations
    # we assume there are no other columns in the character_matrix besides mutations and possibly cluster_id's
    n, m = character_matrix.shape
    
    # compute copy states for each locus
    if isinstance(character_matrix, pd.DataFrame):
        assert np.all(character_matrix.index == copy_states.index) and np.all(character_matrix.columns == copy_states.columns), \
                "The index/columns between the character_matrix and copy_states do not match!"
        assert np.all(copy_states.values >= 0) and np.all(copy_states.values <= 1), "Copy states are not bounded between 0 and 1."
    else:
        if copy_states == "diploid":
            copy_states = character_matrix.values * (1/2)
        elif copy_states == "haploid":
            copy_states  = character_matrix.values
        else:
            print("copy_states is not a valid input!")
            return

    gen = np.random.default_rng(seed)

    ################################
    # Step 1: generate read counts #
    ################################
    total_reads = np.zeros((n, m), dtype=int)
    variant_reads = np.zeros((n, m), dtype=int)
    
    for i,cell in enumerate(character_matrix.index):
        for j,mutation in enumerate(character_matrix.columns):
            
            # same type of model for simulating reads as ConDoR's
            latent_vaf = copy_states.loc[cell,mutation]
            total_reads[i, j] = gen.poisson(mean_coverage)
            post_error_vaf = fp_rate + (1 - fp_rate - fn_rate) * latent_vaf
            ado_alpha = post_error_vaf * ado_precision
            ado_beta = ado_precision * (1 - post_error_vaf)
            variant_reads[i, j] = betabinom.rvs(total_reads[i, j], ado_alpha, ado_beta, random_state=gen)
            
    # generate the binarized mutation matrix
    VAFs = variant_reads / total_reads
    observed_characters = (VAFs >= vaf_threshold).astype(int)
        
    # corrupt read counts and observed character matrix
    variant_reads_corrupt = variant_reads.copy()
    total_reads_corrupt = total_reads.copy()
    
    estimated_unknown_entries = int(np.floor(missing_rate * (n*m)))

    # only corrupt if we have a non-zero missing rate
    if estimated_unknown_entries > 0:
        selected_cells = gen.integers(low=0, high=n, size=estimated_unknown_entries)
        selected_mutations = gen.integers(low=0, high=m, size=estimated_unknown_entries)

        observed_characters[selected_cells, selected_mutations] = unknown_value
        variant_reads_corrupt[selected_cells, selected_mutations] = 0
        total_reads_corrupt[selected_cells, selected_mutations] = 0
    
    # turn everything into a dataframe
    observed_character_matrix = pd.DataFrame(observed_characters, index=character_matrix.index, columns=character_matrix.columns)
    variant_reads = pd.DataFrame(variant_reads, index=character_matrix.index, columns=character_matrix.columns)
    total_reads = pd.DataFrame(total_reads, index=character_matrix.index, columns=character_matrix.columns)
    variant_reads_corrupt = pd.DataFrame(variant_reads_corrupt, index=character_matrix.index, columns=character_matrix.columns)
    total_reads_corrupt = pd.DataFrame(total_reads_corrupt, index=character_matrix.index, columns=character_matrix.columns)

    return observed_character_matrix, variant_reads, total_reads, variant_reads_corrupt, total_reads_corrupt


def simulate(num_mutations, 
             num_cells, 
             num_clusters,
             max_losses=0,
             max_gains=0,
             constraint=LOSS_AND_GAIN,
             loss_prob=0.1,
             gain_prob=0.05,
             mean_coverage=50,
             ado_precision=15,
             missing_rate=0.0,
             fp_rate=0.001,
             fn_rate=0.001,
             vaf_threshold=0.1,
             unknown_value=3,
             max_cn=8,
             seed=None):
    """Simulator for trees with parallel mutations and losses
    
    Parameters
    -----------
    num_mutations: int
        The number of mutations the simulated tree will have
    num_cells: int
        The number of cells that will be simulated in the output data
    num_clusters: int
        The number of distinct clonal populations in the simulation
    K: int, optional
        The number of times a mutation can be lost during evolution (corresponds to the k in k-Dollo). Default = 0
    R: int, optional
        The number of times mutation can be acquired in parallel. Default = 0
    constraints: int, optional
        The type of constraint to apply for mutation gains/losses. Default = LOSS_AND_GAIN
    loss_prob: float, optional
        The probability that a mutation is lost when a new clone is formed. Default = 0.1
    gain_prob: float, optional
        The probability that a mutation occurs in parallel when a new clone is formed. Default = 0.05
    mean_coverage: int, optional
        The mean coverage of sequencing. Default = 50
    fp_rate: float, optional
        The false positive rate during sequencing. Default = 0.001
    fn_rate: float, optional
        The false negative rate during sequencing. Default = 0.001
    ado_precision: int, optional
        The precision parameter for the Beta distribution used to model allelic dropout. Default = 5
    missing_rate: float, optional
        The percentage of mutations that are unknown after sequencing. Default = 0.0
    vaf_threshold: float, optional
        The threshold of variant allele frequency (VAF) to call a mutation present in a cell. Default = 0.01
    max_cn: int, optional
        The maximum number of allele copies at a genomic locus. Default = 2
    seed: int, optional
        The random seed for reproducibility, Default = None
        
    Returns
    --------
    AnnData
        An annotated data structure containing the simulation
    """
    
    # simulate tree
    (T_cell, 
     T_mut, 
     Tc, 
     copy_states, 
     mutant_copy_numbers, 
     total_copy_numbers) = simulate_tree(num_mutations=num_mutations, 
                                         num_cells=num_cells,
                                         num_clusters=num_clusters, 
                                         K=max_losses, 
                                         R=max_gains,
                                         constraint=constraint,
                                         loss_prob=loss_prob,
                                         gain_prob=gain_prob,
                                         max_cn=max_cn,
                                         seed=seed)
    
    # resolve genotypes
    T_cell, character_matrix = op.ul.resolve_genotypes(T_cell)
    
    # simulate reads
    (observed_character_matrix, 
     variant_reads, 
     total_reads, 
     variant_reads_corrupt, 
     total_reads_corrupt) = simulate_reads(character_matrix.replace(op.ul.mutation_types), 
                                           mean_coverage=mean_coverage,
                                           ado_precision=ado_precision,
                                           fp_rate=fp_rate,
                                           fn_rate=fn_rate,
                                           missing_rate=missing_rate,
                                           vaf_threshold=vaf_threshold,
                                           unknown_value=unknown_value,
                                           copy_states=copy_states,
                                           seed=seed)

    adata = ad.AnnData(observed_character_matrix)
    adata.var["mutation_type"] = ["SNV"] * len(observed_character_matrix.columns)  

    # collect all trees
    adata.uns[op.ul.DATA.CELL_TREE] = T_cell
    adata.uns[op.ul.DATA.MUTATION_TREE] = T_mut
    adata.uns[op.ul.DATA.CN_TREE] = Tc
    adata.uns[op.ul.DATA.CLONAL_TREE] = op.ul.to_clonal_tree(T_cell, character_matrix.replace(op.ul.mutation_types))
    
    # collect all cell/mutation data 
    adata.layers[op.ul.DATA.TRUE_DATA] = character_matrix
    adata.layers[op.ul.DATA.OBS_DATA] = observed_character_matrix
    adata.layers[op.ul.DATA.VARIANT_READS] = variant_reads
    adata.layers[op.ul.DATA.TOTAL_READS] = total_reads
    adata.layers[op.ul.DATA.VARIANT_READS_CORRUPT] = variant_reads_corrupt
    adata.layers[op.ul.DATA.TOTAL_READS_CORRUPT] = total_reads_corrupt
    adata.layers[op.ul.DATA.COPY_STATES] = copy_states
    
    # collect cell specific data
    adata.obs[op.ul.DATA.CLUSTER_ID] = T_cell.graph["cluster_assignments"]
    
    # compute FPR and FNR and missing rate
    obs_values = observed_character_matrix.values
    true_values = character_matrix.replace(op.ul.mutation_types).values
    adata.uns[op.ul.DATA.MUTANT_COPY_NUMBERS] = mutant_copy_numbers
    adata.uns[op.ul.DATA.TOTAL_COPY_NUMBERS] = total_copy_numbers
    adata.uns[op.ul.SIM_KEYS.FPR] = np.maximum(1e-6, np.sum((obs_values == 1) & (true_values == 0)) / np.sum(true_values == 0)) # percentage of 0's flipped to 1's
    adata.uns[op.ul.SIM_KEYS.FNR] = np.maximum(1e-6, np.sum((obs_values == 0) & (true_values == 1)) / np.sum(true_values == 1)) # percentage of 1's flipped to 0's
    adata.uns[op.ul.SIM_KEYS.MISSING_RATE] = np.sum(obs_values == 3) / obs_values.size # percentage of entries flipped to 3

    return adata
    