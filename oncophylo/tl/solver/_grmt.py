# _grmt.py 
import shutil, os
import oncophylo as op 

def GRMT(input_df, 
          l=0.7,
          K=1.0,
          k=0,
          t=1,
          fp=0.001,
          fn=0.2,
          n=100,
          N=30,
          remove_temp_dir=True):
    """
    GRMT Python Wrapper
    
    Parameters
    ----------
    input_df: pd.DataFrame
        A Pandas dataframe where the row are cells and the columns are mutations
    l: float
        The lambda parameter. Should be between 0.5 and 1.0. (Default = 0.7)
    K: float
        The Kappa parameter. Should be non-negative. (Default = 1.0)
    k: int
        The maximum number of times a mutation can be lost. (Default = 1)
    t: int
        The number of threads to use. (Default = 1)
    fp: float
        The estimated false positive rate for the sequencing data
    fn: float
        The estimated false negative rate for the sequencing data
    n: int
        The desired number of data points to sample in the Bayesian Optimization algorithm. (Default = 100)
    N: int
        The desired number of data points to evaluate the BO algorithm. (Default = 30)
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files. Default = True

    Returns
    --------
    list or dictionary
        A dictionary or a list of dictionaries containing solutions and their corresponding meta data
    """

    # create folder to store temp files
    temp_path = os.path.join(os.path.abspath(""),"grmt_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "input.txt")
    output_fn = os.path.join(temp_path, "output")
    cell_names_fn = os.path.join(temp_path, "cell_names.txt")
    gene_names_fn = os.path.join(temp_path, "gene_names.txt")
    
    # create temp files
    input_df = input_df.replace(-1,3)
    input_df.to_csv(input_fn, sep="\t", index=False, header=False)
    input_df.index.to_series().to_csv(cell_names_fn, index=False, header=False)
    input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)

    cells = input_df.index.to_series().values
    
    args = [
            "-l", "%.6f" % l,
            "-K", "%.6f" % K,
            "-k", "%d" % k,
            "-t", "%d" % t,
            "-i", "%s" % input_fn,
            "-o", "%s" % output_fn,
            "-c", "%s" % cell_names_fn,
            "-m", "%s" % gene_names_fn,
            "-a", "%.6f" % fp,
            "-b", "%.6f" % fn,
            "-n", "%d" % n,
            "-N", "%d" % N,
    ]

    # run GRMT
    output, time = op.ul.subprocess([op.ul.binary_path("grmt")] + args)

    # collect solutions
    dot_fn = os.path.join(output_fn + ".dot")
    mapping = {'s%d' % i:cells[i] for i in range(len(cells))}
    mapping[""] = "root"
    T_cell, T_mut = op.io.load_dot(dot_fn, 
                                   mutations = list(input_df.columns), 
                                   cells = list(input_df.index), 
                                   mapping=mapping,
                                   _type="cell_tree",
                                   set_id_to_label=True,
                                   preprocessor=preprocess_GRMT)
    (T_cell, output_df) = op.ul.resolve_genotypes(T_cell, input_df)
    solution = op.ul.solution(T_cell, 
                              T_mut, 
                              input_df,
                              output_df, 
                              fp,
                              fn,
                              output,
                              time)
    if remove_temp_dir:
        shutil.rmtree(temp_path)
        
    return solution
    
def preprocess_GRMT(T):
    """Preprocesses the GRMT tree by splitting apart the cells into distinct nodes"""
    T_nodes = list(T.nodes())
    for n in T_nodes:
        cell_list = n.split(" ")
        if len(cell_list) > 1:
            p = list(T.predecessors(n))[0]
            for s in cell_list:
                T.add_edge(p,s)
            T.remove_node(n)
    return T
    