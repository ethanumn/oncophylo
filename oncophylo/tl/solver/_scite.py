# SCITE 
import shutil, os
import oncophylo as op 

def SCITE(input_df, 
          repetitions=3, 
          chain_length=90000, 
          fp=0.001, 
          fn=0.2, 
          e=0,
          n_solutions=1, 
          seed=None, 
          remove_temp_dir=True):
    """
    SCITE Python Wrapper
    
    Parameters
    ----------
    input_df: pd.DataFrame
        A Pandas dataframe where the row are cells and the columns are mutations
    repetitions: int
        The number of restarts to run (-r for SCITE)
    chain_length: int
        The number of trees to sample during MCMC (-l for SCITE)
    fp: float
        The false positive rate (-fd for SCITE). Default = 0.001
    fn: float
        The false negative rate (-ad for SCITE). Default = 0.2
    e: float
        Error rate jump probability (-e for SCITE). Default = 0
    n_solutions: int
        The number of solutions to return. Default = 1
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files. Default = True

    Returns
    --------
    list or dictionary
        A dictionary or a list of dictionaries containing solutions and their corresponding meta data
    """

    # create folder to store temp files
    temp_path = os.path.join(os.path.abspath(""),"scite_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "input.txt")
    output_fn = os.path.join(temp_path, "output")
    gene_names_fn = os.path.join(temp_path, "gene_names.txt")
    
    # create temp files
    m, n = input_df.shape
    input_df = input_df.replace(-1,3)
    input_df.T.to_csv(input_fn, sep=" ", index=False, header=False)
    input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)
    cells = input_df.index.to_series().values
    
    args = [
            "-r", "%s" % repetitions,
            "-l", "%d" % chain_length,
            "-n", "%d" % n,
            "-m", "%d" % m,
            "-i", "%s" % input_fn,
            "-o", "%s" % output_fn,
            "-names", "%s" % gene_names_fn,
            "-fd", "%.6f" % fp,
            "-ad", "%.6f" % fn,
            "-e", "%.4f" % e,
            "-seed", "%d" % seed if seed is not None else "",
            "-a"
    ]

    # run SCITE
    output, time = op.ul.subprocess([op.ul.binary_path("SCITE")] + args)

    # collect solutions
    solutions = []
    for i in range(n_solutions):  
        dot_fn = os.path.join(output_fn + "_ml%d.gv" % i,)
        mapping = {'s%d' % i:cells[i] for i in range(len(cells))}
        T_cell, T_mut = op.io.load_dot(dot_fn, 
                                       mutations = list(input_df.columns), 
                                       cells = list(input_df.index), 
                                       mapping=mapping,
                                       _type="cell_tree")
        (T_cell, output_df) = op.ul.resolve_genotypes(T_cell, input_df)
        solutions.append(op.ul.solution(T_cell, 
                                        T_mut, 
                                        input_df,
                                        output_df, 
                                        fp,
                                        fn,
                                        output,
                                        time))
    if remove_temp_dir:
        shutil.rmtree(temp_path)
        
    if n_solutions == 1:
         return solutions[0]
    else:
        return solutions
    