# _ltorchard.py

import os, sys, shutil
import pandas as pd 
import numpy as np 

import oncophylo as op 

def ltOrchard(character_matrix, 
              k=10, 
              e=10, 
              fp=0.001, 
              fn=0.2, 
              only_add_leaves=False,
              greedy=False,
              min_cells=1,
              missing_character=-1,
              seed=None, 
              n_solutions=1, 
              remove_temp_dir=True):
    """
    scOrchard Python Wrapper
    
    Parameters
    ----------
    input_df: pd.DataFrame
        A Pandas dataframe where the row are cells and the columns are mutations
    k: int
        The beam width -- which controls the number of unique tree returned
    e: int 
        The branching factor -- which controls how many solutions are expanded from the queue at each iteration
    fp: float
        The false positive rate to use for likelihood calculations
    fn: float
        The false negative rate to use for likelihood calculations
    only_add_leaves: bool
        Flag to only consider adding mutations as leaves in the tree. Default = False
    greedy: bool
        Flag to run scOrchard in its greedy search mode. Default = False
    min_cells: int
        Remove all mutations that occur in this or fewer number of cells
    missing_character: int
        The value that represents a missing value to Cassiopeia
    seed: int
        The random seed to use
    n_solutions: int
        The number of solutions to return. Default = 1
    remove_temp_dir: bool
        Flag to remove the temporary directory used to store scOrchard's input/output files. Default = True

    Returns
    --------
    dictionary
        A dictionary of results including the mutation tree, clonal tree, run time, and likelihood
        for each returned solution.
    """
    
    # preprocessing 
    input_df, binary_character_matrix,  mapping = op.ul.prep_lineage_tracing_data(character_matrix, 
                                                                                  min_cells=min_cells, 
                                                                                  missing_character=missing_character) 
    
    muts_to_add = binary_character_matrix.loc[:,np.sum(binary_character_matrix == 1, axis=0) == 1]

    
    temp_path = os.path.join(os.path.abspath(""),"ltOrchard_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "input.txt")
    cell_names_fn = os.path.join(temp_path, "cell_names.txt")
    gene_names_fn = os.path.join(temp_path, "gene_names.txt")

    output_path = os.path.join(temp_path)
    output_prefix = "out"
    n, m = input_df.shape
        
    input_df.T.to_csv(input_fn, sep=" ", index=False, header=False)
    input_df.index.to_series().to_csv(cell_names_fn, index=False, header=False)
    input_df.columns.to_series().to_csv(gene_names_fn, index=False, header=False)
    
    args = [
            "-k", "%s" % k,
            "-e", "%d" % e,
            "-n", "%d" % n,
            "-m", "%d" % m,
            "-i", "%s" % input_fn,
            "-o", "%s" % output_path,
            "-p", "%s" % output_prefix,
            "-g", "%s" % gene_names_fn,
            "-c", "%s" % cell_names_fn,
            "-fd", "%.6f" % fp,
            "-ad", "%.6f" % fn,
            "-seed", "%d" % seed if seed is not None else "",
            "-greedy" if greedy else "",
            "-onlyAddLeaves" if only_add_leaves else "",
    ]

    # run ltOrchard
    output, time = op.ul.subprocess([op.ul.binary_path("ltOrchard")] + args)

    # process solutions
    solutions = []
    for i in range(n_solutions):  
        dot_fn = os.path.join(output_path, output_prefix + "_ml%d.gv" % i)
        T_cell, T_mut = op.io.load_dot(dot_fn, 
                     mutations = list(input_df.columns), 
                     cells = list(input_df.index), 
                     _type="cell_tree")
        (T_cell, predicted_genotypes) = op.ul.resolve_genotypes(T_cell, input_df)
        reconstructed_tree = op.ul.post_process_celltree(T_cell, character_matrix, muts_to_add, mapping)
        solutions.append(reconstructed_tree)
    
    if remove_temp_dir:
        shutil.rmtree(temp_path)
        
    if n_solutions == 1:
         return solutions[0]
    else:
        return solutions