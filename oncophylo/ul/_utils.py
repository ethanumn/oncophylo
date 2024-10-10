# utils.py
import shutil, os
import subprocess as sp
from time import perf_counter

import oncophylo as op 

def binary_path(binary_name):
    """Finds a binary executable path. First searches the local machine, then looks
    in the OncoPhylo library
    
    Input
    ------
    binary_name: str
        The name of the binary file on the system
    Returns
    --------
    str
        The path of the binary executable
    """
    ex = shutil.which(binary_name)
    if ex is None:
        ex = os.path.join(op.ul.BIN_PATH, binary_name)
        if os.path.exists(ex):
            return ex
        else:
            print(f"Cannot locate binary executable at {ex}")
    return ex

def script_path(library, script_name):
    """Finds a script within a library
    
    Input
    ------
    script_name: str
        The relative path inside the library where the script is found
    library: str
        The name of the libary. Looks in /path/to/oncophylo/lib/library/script_name

    Returns
    --------
    str
        The path to the script
    """
    ex = os.path.join(op.ul.LIB_PATH, library, script_name)
    if os.path.exists(ex):
        return ex
    else:
        print(f"Cannot locate script at {ex}")

def subprocess(call, stdout=None):
    """
    Input
    ------
    call: str
        A call to a binary executable 

    Returns
    --------
    str
        The terminal output of the program
    float
        The run time, in seconds, the executable took to complete
    """
    start = perf_counter()
    res = None
    if stdout == None:
        res = sp.run(call, capture_output=True)
    else:
        sp.run(call, stdout=stdout)
    end = perf_counter()

    return res, end-start

def solution(T_cell, 
             T_mut, 
             character_matrix,
             output_df, 
             fp,
             fn,
             output,
             time,
             T_clonal=None,
             var_reads=None,
             total_reads=None,
             ado_precision=None):
    """Returns a dictionary with all of the outputs and performance metrics for a solution. Depending on which inputs are provided,
    different measurements will available in the dictionary return by this function.
    
    Input
    ------
    T_cell: Networkx.DiGraph
        A cell tree where internal nodes are mutations and leaves are cells
    T_mut: Networkx.DiGraph
        A mutation tree
    fp: float
        The false positive rate for scoring the character matrix
    fn: float
        The false negative rate for scoring the character matrix
    character_matrix: pd.DataFrame
        The character matrx input the algorithm
    output_df: pd.DataFrame
        The predicted character matrix
    output: str
        The subprocess output (i.e., the terminal output)
    time: float
        The amount of time in seconds it took for the subprocess to complete
    T_clonal: Networkx.DiGraph, optional
        A clonal tree where nodes are groups of cells and edges are groups of mutations.
        Some methods output a clonal tree, so providing this will make it so a a clonal tree isn't computed.
    var_reads: pd.DataFrame, optional
        A cell by mutation matrix of variant read counts 
    total_reads: pd.DataFrame, optional
        A cell by mutation matrix of total read counts 
    ado_precision: float, optional
        The allelic dropout precision parameter. This (along with the var_reads and total_reads) is used to calculated the likelihood of the tree under a beta binomial model.
    """
    return {op.ul.DATA.CELL_TREE: T_cell, 
            op.ul.DATA.CLONAL_TREE: op.ul.to_clonal_tree(T_cell, output_df) if T_clonal is None else T_clonal,
            op.ul.DATA.MUTATION_TREE:T_mut,
            op.ul.DATA.PRED_DATA: output_df,
            op.ul.EVAL_KEYS.RUNTIME:time,
            op.ul.EVAL_KEYS.MATRIX_ERROR: op.tl.score.matrix_error(output_df, character_matrix),
            op.ul.EVAL_KEYS.LLH_OE: op.tl.score.score_observation_errors(output_df, character_matrix, fp, fn),
            op.ul.EVAL_KEYS.LLH_BB: op.tl.score.score_beta_binomial(output_df, var_reads, total_reads, ado_precision, fp),
            op.ul.DATA.TERMINAL_OUTPUT:output}