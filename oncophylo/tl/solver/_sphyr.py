# SPhyR 
import shutil, os
import oncophylo as op 

def SPhyR(input_df, 
          k=0, 
          N=10,
          fp=0.001, 
          fn=0.1, 
          lC=15,
          lT=10,
          seed=0,
          num_threads=1,
          remove_temp_dir=True):
    
    """Python wrapper for SPhyR
    
    
    Parameters
    ----------
    input_df: pd.DataFrame
        A Pandas dataframe where rows are cells and columns are mutations
    k: int
        The number of mutation losses to allow in the k-Dollo model
    N: int, optional
        Number of restarts. Default = 10
    fp: float, optional
        The false positive rate. Default = 0.001
    fn: float, optional
        The false negative rate. Default = 0.1 
    lC: int, optional
        The number of character clusters. Default = 15
    lT: int, optional
        The number of taxon clusters. Default = 10
    seed: int, optional
        The random seed. Default = 0
    num_threads: int, optional
        The number of threads to use. Default = 1
    remove_temp_dir: bool, optional
        Flag which when true removes the temporary directory that contains the output files. Default = True
        
    Returns
    -------
    dictionary
        A dictionary of results
    """
    
    temp_path = os.path.join(os.path.abspath(""),"sphyr_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "input.txt")
    gene_names_fn = os.path.join(temp_path, "genes.labels")
    cell_names_fn = os.path.join(temp_path, "cells.labels")
    output_fn = os.path.join(temp_path, "output.txt")
    dot_fn = os.path.join(temp_path, "output.dot")

    n, m = input_df.shape
    
    # replace 3 with -1
    input_df.replace(3, -1).to_csv(input_fn, sep=" ", index=False, header=False)
    
    # write first two lines to tell SPhyR how many mutations/cells there are
    with open(input_fn, 'r') as original:
        data = original.read()
        
    with open(input_fn, 'w') as modified:
        modified.write("%d #taxa\n%d #characters\n" % (n, m) + data)
    
    input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)
    input_df.index.to_series().to_csv(cell_names_fn, index=False, header=False)
    
    args = ["%s" % input_fn,
            "-k", "%s" % k,
            "-N", "%d" % N,
            "-a", "%.6f" % fp,
            "-b", "%.6f" % fn,
            "-lC", "%d" % lC,
            "-lT", "%d" % lT,
            "-t", "%d" % num_threads,
            "-s", "%d" % seed 
    ]

    # capture output from SPhyR
    f = open(output_fn, "w")

    # run SPhyR
    output, time = op.ul.subprocess([op.ul.binary_path("kDPFC")] + args, stdout=f)
    f.close()
    
    args = ["%s" % output_fn,
            "-c", "%s" % gene_names_fn,
            "-t", "%s" % cell_names_fn
    ]
    
    # write dot file
    f = open(dot_fn, "w")
    op.ul.subprocess([op.ul.binary_path("visualize")] + args, stdout=f)
    f.close() 
    
    # prepare output dataframe
    predicted_genotypes = pd.read_csv(output_fn, skiprows=2, names=list(input_df.columns), index_col=None, sep=" ")
    predicted_genotypes.index = list(input_df.index)

    T, _ = load_dot(dot_fn, 
                 mutations = list(input_df.columns), 
                 cells = list(input_df.index), 
                 _type="clonal_tree")
    
    T_clonal = post_process_SPhyR(T, predicted_genotypes.replace(mutation_types))
            
    solution  = {CONST.CELL_TREE: None, 
                 CONST.CLONAL_TREE:T_clonal,
                 CONST.MUTATION_TREE: None,
                 CONST.PRED_DATA:predicted_genotypes,
                 CONST.RUNTIME:end-start,
                 CONST.LLH: score_genotypes(predicted_genotypes.replace(mutation_types), input_df, fn, fp),
                 CONST.TERMINAL_OUTPUT:res}
    
    return solution

def post_process_SPhyR(T, df):
    """Perform post-processing on SPhyR's tree output. Translates it into a clonal_tree."""
    T_clonal = nx.DiGraph()
    T_clonal.graph["data"] = df
    T_clonal.graph["splitter_mut"] = "\n"
    T_clonal.graph["splitter_cell"] = "\n"
    T_clonal.graph["loss_prefix"] = "-"
    T_clonal.graph["become_germline"] = list(df.columns[(df == 0).all(axis=0)])
    T_clonal.graph["type"] = CONST.CLONAL_TREE
    T_clonal.graph["cells"] = list(df.index.values)
    T_clonal.graph["mutations"] = list(df.columns.values)
    T_clonal.graph["losses"] = []
        
    # get root for finding predecessor
    root = root_id(T)
            
    # process nodes
    new_nodes = {}
    to_remove = []
    
    # relabel nodes with mutations, and attach all cells as individual nodes to the
    # correct mutation
    for v,label in T.nodes().data("label"):
        label_pieces = []
        for s in label.replace('"', '').split("\\n"):
            label_pieces += s.split(" ")
        predecessors = nx.predecessor(T, source=root, target=v)
        if np.all(np.in1d(label_pieces,T.graph["cells"])):
            if len(predecessors) == 1:
                u = predecessors[0]
                T.nodes[u]["label"] = "\n".join(label_pieces)
                to_remove.append(v)
                
    # remove nodes that we don't need any more
    for v in to_remove:
        T.remove_node(v)
        
    # any node that doesn't have a label, label it with "––"
    for v,label in T.nodes().data("label"):
        if label == '""':
            T.nodes[v]["label"] = "––"
            
    # edge delimiter is wrong so we need to correct it
    for u,v,label in T.edges.data("label"):
        edge_muts = label.replace('"', "").replace("--", "-").split("\\n")
        for mut in edge_muts:
            if mut.startswith(T_clonal.graph["loss_prefix"]):
                T_clonal.graph["losses"].append(mut)
        T[u][v][0]["label"] = "\n".join(edge_muts)
           
    # copy everything into new tree
    for v,label in T.nodes().data("label"):
        T_clonal.add_node(v, label=label)

    for (u,v,data) in T.edges(data=True):
        T_clonal.add_edge(u,v,label=data["label"])      
        
    return T_clonal