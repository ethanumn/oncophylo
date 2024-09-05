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

def subprocess(call):
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
    res = sp.run(call, capture_output = True)
    end = perf_counter()

    return res, end-start

def solution(T_cell, 
             T_mut, 
             input_df,
             output_df, 
             fp,
             fn,
             output,
             time):
    """Returns a dictionary with all of the outputs and performance metrics for a solve
    
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
    input_df: pd.DataFrame
        The character matrx input the algorithm
    output_df: pd.DataFrame
        The predicted character matrix
    output: str
        The subprocess output (i.e., the terminal output)
    time: float
        The amount of time in seconds it took for the subprocess to complete
    """
    return {op.ul.CONST.CELL_TREE: T_cell, 
            op.ul.CONST.CLONAL_TREE:op.ul.to_clonal_tree(T_cell, output_df),
            op.ul.CONST.MUTATION_TREE:T_mut,
            op.ul.CONST.PRED_DATA: output_df,
            op.ul.CONST.RUNTIME:time,
            op.ul.CONST.MATRIX_ERROR: op.tl.score.matrix_error(output_df, input_df),
            op.ul.CONST.LLH: op.tl.score.score_genotypes(output_df, input_df, fp, fn),
            op.ul.CONST.TERMINAL_OUTPUT:output}