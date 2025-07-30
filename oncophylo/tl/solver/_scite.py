# SCITE 
import shutil, os
import oncophylo as op 

def SCITE(character_matrix, 
          repetitions=3, 
          chain_length=90000, 
          fp=0.001, 
          fn1=0.2, 
          fn2=0.0, 
          cc=0.0,
          e=0,
          n_solutions=1, 
          seed=None, 
          remove_temp_dir=True,
          destination_dir = ""):
    """
    SCITE Python Wrapper
    
    Parameters
    ----------
    character_matrix: pd.DataFrame
        A character matrix where the row are cells and the columns are mutations
    repetitions: int
        The number of restarts to run (-r for SCITE).
    chain_length: int
        The number of trees to sample during MCMC (-l for SCITE).
    fp: float
        The false positive rate (-fd for SCITE). Default = 0.001.
    fn1: float
        The false negative rate for heterozygous loci (-ad for SCITE). Default = 0.2.
    fn2: float
        The flase negative rate for homozygous loci. Default = 0.0.
    cc: float
        Estimated rate of non-mutated sites called as homozygous mutations. Default = 0.0.
    e: float
        Error rate jump probability (-e for SCITE). Default = 0.
    n_solutions: int
        The number of solutions to return. Default = 1.
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files. Default = True.

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

    fn_params = f"{fn1}"
    if fn2 > 0.0:
        fn_params = f"{fn1} {fn2}"

    args_dict = {
        "-r": str(repetitions),
        "-l": str(chain_length),
        "-n": str(n),
        "-m": str(m),
        "-i": input_fn,
        "-o": output_fn,
        "-names": gene_names_fn,
        "-fd": f"{fp:.6f}",
        "-ad": fn_params,
        "-cc": f"{cc:.6f}",
        "-e": f"{e:.4f}",
    }

    if seed is not None:
        args_dict["-seed"] = str(seed)

    args = op.ul.convert_args(args_dict) # convert to list for subprocess
    args.append("-a")

    # run SCITE
    output, time = op.ul.subprocess([op.ul.binary_path("SCITE")] + args)

    # collect solutions
    solutions = []
    for i in range(n_solutions):  
        dot_fn = os.path.join(output_fn + "_ml%d.gv" % i,)
        mapping = {'s%d' % i:cells[i] for i in range(len(cells))}
        T_cell, T_mut = op.io.load_dot(dot_fn, 
                                       mutations = list(character_matrix.columns), 
                                       cells = list(character_matrix.index), 
                                       mapping=mapping,
                                       _type="cell_tree")
        (T_cell, corrected_character_matrix) = op.ul.resolve_genotypes(T_cell, character_matrix)
        solutions.append(op.ul.solution(T_cell, 
                                        T_mut, 
                                        character_matrix,
                                        corrected_character_matrix, 
                                        output,
                                        time))

        # save output files into a directory if provided a valid directory
        op.ul.save_output_files(destination_dir, [dot_fn])

    if remove_temp_dir:
        shutil.rmtree(temp_path)
        
    if n_solutions == 1:
         return solutions[0]
    else:
        return solutions
    