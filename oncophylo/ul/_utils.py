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

def convert_args(args_dict):
    """Converts arguments stored in a dictionary into a vector that can be passed to subprocess
    
    Input
    -----
    args_dict: dict
        Dictionary where keys are string
    """
    args_list = []
    for k, v in args_dict.items():
        if k != "" and v != "":
            args_list.append(k)
            args_list.append(v)
    return args_list

def save_output_files(destination_dir, file_paths):
    """Copies files into a directory

    Input
    ------
    destination_dir: str
        A directory to copy all files in file_paths to
    file_paths: list
        A list of file paths to copy
    """
    if destination_dir != "":
        if not os.path.exists(destination_dir):
            os.mkdir(destination_dir)

        for file_path in file_paths:
            if not os.path.exists(file_path):
                raise ExceptionType(f"{file_path} does not exist, yet we are trying to copy it")
            shutil.copy(file_path, destination_dir)

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

def solution(T_cell, 
             T_mut, 
             character_matrix,
             corrected_character_matrix, 
             output,
             time,
             T_clonal=None):
    """Returns a dictionary with all of the outputs and performance metrics for a solution. Depending on which inputs are provided,
    different measurements will available in the dictionary return by this function.
    
    Input
    ------
    T_cell: Networkx.DiGraph
        A cell tree where internal nodes are mutations and leaves are cells
    T_mut: Networkx.DiGraph
        A mutation tree
    character_matrix: pd.DataFrame
        The character matrx input the algorithm
    corrected_character_matrix: pd.DataFrame
        The corrected character matrix after reconstruction
    output: str
        The subprocess output (i.e., the terminal output)
    time: float
        The amount of time in seconds it took for the subprocess to complete
    T_clonal: Networkx.DiGraph, optional
        A clonal tree where nodes are groups of cells and edges are groups of mutations.
        Some methods output a clonal tree, so providing this will make it so a a clonal tree isn't computed.
    """
    return {op.ul.DATA.CELL_TREE: T_cell, 
            op.ul.DATA.CLONAL_TREE: op.ul.to_clonal_tree(T_cell, output_df) if T_clonal is None else T_clonal,
            op.ul.DATA.MUTATION_TREE:T_mut,
            op.ul.DATA.PRED_DATA: corrected_character_matrix,
            op.ul.EVAL_KEYS.RUNTIME:time,
            op.ul.DATA.TERMINAL_OUTPUT:output}


def subprocess(call, stdout=None, stderr=None):
    """
    Input
    ------
    call: str
        A call to a binary executable 
    stdout: str
        Fil

    Returns
    --------
    str
        The terminal output of the program
    float
        The run time, in seconds, the executable took to complete
    """
    start = perf_counter()
    res = None
    if stdout == None and stderr == None:
        res = sp.run(call, capture_output=True, text=True)
    else:
        sp.run(call, stdout=stdout, stderr=stderr)
    end = perf_counter()

    return res, end-start