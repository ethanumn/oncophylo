# HUNTRESS
import shutil, os
import pandas as pd 
import oncophylo as op 

def HUNTRESS(input_df, 
             fp=0.001, 
             fn=0.2, 
             num_threads=1, 
             algorithmchoice="FPNA",
             remove_temp_dir=True):

    temp_path = os.path.join(os.path.abspath(""),"huntress_temp")
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)
    input_fn = os.path.join(temp_path, "input.SC")
    output_fn = os.path.join(temp_path, "output")
    input_df = input_df.replace(-1,3)
    input_df.to_csv(input_fn, sep="\t")
    
    args = [
            "--i=%s" % input_fn,
            "--o=%s" % output_fn,
            "--t=%d" % num_threads, 
            "--algorithmchoice=%s" % algorithmchoice,
            "--fp_coeff=%.6f" % fp,
            "--fn_coeff=%.6f" % fn
           ]
    
    # run HUNTRESS
    output, time = op.ul.subprocess(["python", op.ul.script_path("HUNTRESS", "HUNTRESS.py")] + args)

    output_df = pd.read_csv(output_fn + ".CFMatrix", sep="\t", header=0, index_col=0)
    T_clonal = op.ul.conflict_free_matrix_to_clonal_tree(output_df)
    T_cell, T_mut = op.ul.clonal_to_cell_tree(T_clonal)
    
    solution  = op.ul.solution(T_cell, 
                               T_mut, 
                               input_df,
                               output_df, 
                               fp,
                               fn,
                               output,
                               time,
                               T_clonal=T_clonal)
                               
    if remove_temp_dir:
        shutil.rmtree(temp_path)
    
    return solution