# infSCITE.py
import shutil, os
import oncophylo as op 

def infSCITE(input_df, 
             repetitions=1, 
             chain_length=90000, 
             fp=0.001, 
             fn=0.2, 
             e=0,
             density=10000,
             rec="",
             infer_doublets=False,
             doublet_rate=0.0,
             n_solutions=1, 
             seed=None, 
             remove_temp_dir=True):
    """
    infSCITE Python Wrapper
    
    Parameters
    ----------
    input_df: pd.DataFrame
        A Pandas dataframe where the row are cells and the columns are mutations
    repetitions: int
        The number of restarts to run (-r for infSCITE)
    chain_length: int
        The number of trees to sample during MCMC (-l for infSCITE)
    fp: float
        The false positive rate (-fd for infSCITE). Default = 0.001
    fn: float
        The false negative rate (-ad for infSCITE). Default = 0.2
    e: float
        Error rate jump probability (-e for infSCITE). Default = 0
    density: int
        Number of samples from the posterior distribution (-p for infSCITE). Default = 10000
    rec: str
        The name of the mutation to evaluate as recurrent. This should match the name in the input_df
    infer_doublets: bool
        Flag to tell infSCITE to infer the doublet rate (-d for infSCITE). Default = False
    doublet_rate: float
        A provided doublet rate (passed to -d for infSCITE). Default = 0
    n_solutions: int
        The number of solutions to return. Default = 1
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files. Default = True

    Returns
    --------
    list or dictionary
        A dictionary or a list of dictionaries containing solutions and their corresponding meta data
    """
    
    temp_path = os.path.join(os.path.abspath(""),"infSCITE_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "input.txt")
    output_fn = os.path.join(temp_path, "output")
    gene_names_fn = os.path.join(temp_path, "gene_names.txt")
    cell_names_fn = os.path.join(temp_path, "cell_names.txt")
    
    m, n = input_df.shape
    input_df.T.to_csv(input_fn, sep=" ", index=False, header=False)
    input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)
    input_df.index.to_series().to_csv(cell_names_fn, index=False, header=False)
    cells = input_df.index.to_series().values
    
    args = [
            "-r", "%s" % repetitions,
            "-l", "%d" % chain_length,
            "-n", "%d" % n,
            "-m", "%d" % m,
            "-i", "%s" % input_fn,
            "-o", "%s" % output_fn,
            "-names", "%s" % gene_names_fn,
            "-samples", "%s" % cell_names_fn,
            "-fd", "%.6f" % fp,
            "-ad", "%.6f" % fn,
            "-seed", "%d" % seed if seed is not None else "",
            "-s",
            "-a",
            "-e", "%.4f" % e,
            "-p", "%d" % density
    ]
    
    # add recurrent mutation
    if rec != "":
        args += ["-rec", "%d" % (input_df.columns.get_loc(rec) + 1)]
        
    # add doublet detection
    if doublet_rate > 0:
        args += ["-d", "%.4f" % doublet_rate]
    elif infer_doublets:
        args += ["-d"] 

    print(args)
    # run infSCITE
    output, time = op.ul.subprocess([op.ul.binary_path("infSCITE")] + args)

    solutions = []
    for i in range(n_solutions):  
        dot_fn = os.path.join(output_fn + "_map%d.gv" % i)
        mapping = {'s_%d' % i:cells[i] for i in range(len(cells))}
        mapping["%s_copy" % rec] = "+%s (1)" % rec
        T_cell, T_mut = op.io.load_dot(dot_fn, 
                                       mutations = list(input_df.columns), 
                                       cells = list(input_df.index), 
                                       mapping=mapping,
                                       _type="cell_tree",
                                       set_id_to_label=True)
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