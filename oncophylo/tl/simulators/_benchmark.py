# _benchmark.py

import os, sys
import pandas as pd
import numpy as np
from itertools import product 
from dataclasses import fields

import oncophylo as op 

def load_sim_params(params, cells_geq_mutations = True, generate_seed=True):
    """Loads simulation parameters"""
    assert np.all(field.name in params for field in fields(op.ul.SIM_KEYS)), "params must contain all keys in op.ul.SIM_KEYS" 
    key_order = []
    data = []
    for k, v in params.items():
        if type(v) == list:
            data.append(v)
        else:
            data.append([v])
        key_order.append(k)
        
    res = []
    for d in product(*data):
        setup = {k:v for k,v in zip(key_order, d)}
        for i in range(1, setup[op.ul.SIM_KEYS.TRIALS]+1):
            setup_i = setup.copy()
            setup_i[op.ul.SIM_KEYS.TRIAL_NUM] = i
            if cells_geq_mutations:
                if setup_i[op.ul.SIM_KEYS.CELLS] < setup_i[op.ul.SIM_KEYS.MUTATIONS]:
                        continue
            if generate_seed:
                setup_i[op.ul.SIM_KEYS.SEED] = np.random.randint(2**32-1)
            res.append(setup_i)
            
    return res

def evaluate(adata, dataset_name, model_params, results):
    """
    Evaluates methods on a data set. 

    Parameters
    -----------
    adata: AnnData
        An annotated data object that includes all relevant data for a single-cell data set
    dataset_name: str
        The name of the data set. This is used to log the evaluation results into the results dictionary
    model_params: dict
        A dictionary where keys are the name of the model, and the values are the model parameters
    results: dict
        A dictionary to store the evaluation results. Keys must be in op.ul.EVAL_KEYS and values must lists

    Returns
    --------
    None
    """
    valid_keys = [field.default for field in fields(op.ul.EVAL_KEYS)] + [field.default for field in fields(op.ul.SIM_KEYS)]
    assert isinstance(results, dict), f"results must be a dictionary but was given {type(results)}"
    for k,v in results.items():
        assert isinstance(k, str), f"Keys of results must be str, not {type(k)}"
        assert k in valid_keys, f"{k} is not in op.ul.EVAL_KEYS or op.ul.SIM_KEYS" 
        assert isinstance(v, list), f"Keys of results must be list, not {type(v)}"

    print("Running models on %s" % dataset_name)

    character_matrix = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
    n, m = character_matrix.shape 

    # add cluster id if it's the adata
    if op.ul.DATA.CLUSTER_ID in adata.obs:
        character_matrix[op.ul.DATA.CLUSTER_ID] = adata.obs[op.ul.DATA.CLUSTER_ID]
        
    T_true = adata.uns[op.ul.DATA.CLONAL_TREE]
    
    for method, inputs in model_params.items():
        
        input_df = character_matrix.copy()
        
        method_params = inputs["fixed"]
            
        # extract parameters from layers
        for k,v in inputs["layers"].items():
            method_params[k] = pd.DataFrame(adata.layers[v], index=adata.obs.index, columns=adata.var.index)
                
        # extract parameters from uns
        for k,v in inputs["uns"].items():
            method_params[k] = adata.uns[v]
                    
        if op.ul.DATA.REQUIRES_CLUSTERS in inputs:
            if inputs[op.ul.DATA.REQUIRES_CLUSTERS]:
                pass
            else:
                input_df = input_df.drop(columns=op.ul.DATA.CLUSTER_ID)
        else:
            input_df = input_df.drop(columns=op.ul.DATA.CLUSTER_ID)
            
        output = inputs["func"](input_df, **method_params)
        output_df = output[op.ul.DATA.PRED_DATA]
        T_pred = output[op.ul.DATA.CLONAL_TREE]
        
        # collect results
        if op.ul.EVAL_KEYS.MODEL in results:
            results[op.ul.EVAL_KEYS.MODEL].append(method)
        if op.ul.EVAL_KEYS.PAIRWISE_REL_ACC in results:
            results[op.ul.EVAL_KEYS.PAIRWISE_REL_ACC].append(op.tl.score.pairwise_rel_accuracy(T_pred, T_true))
        if op.ul.EVAL_KEYS.MATRIX_ERROR in results:
            results[op.ul.EVAL_KEYS.MATRIX_ERROR].append(op.tl.score.matrix_error(output_df, input_df))
        if op.ul.EVAL_KEYS.RUNTIME in results:
            results[op.ul.EVAL_KEYS.RUNTIME].append(output[op.ul.EVAL_KEYS.RUNTIME])
        if op.ul.EVAL_KEYS.DATASET in results:
            results[op.ul.EVAL_KEYS.DATASET].append(dataset_name)
        if op.ul.SIM_KEYS.MUTATIONS in results:
            results[op.ul.SIM_KEYS.MUTATIONS].append(m)
        if op.ul.SIM_KEYS.CELLS in results:
            results[op.ul.SIM_KEYS.CELLS].append(n)
        if op.ul.SIM_KEYS.FNR in results:
            results[op.ul.SIM_KEYS.FNR].append(adata.uns[op.ul.SIM_KEYS.FNR])
        if op.ul.SIM_KEYS.FPR in results:
            results[op.ul.SIM_KEYS.FPR].append(adata.uns[op.ul.SIM_KEYS.FPR])
        if op.ul.EVAL_KEYS.LLH_BB in results:
            results[op.ul.EVAL_KEYS.LLH_BB].append(output[op.ul.EVAL_KEYS.LLH_BB])
        if op.ul.EVAL_KEYS.LLH_OE in results:
            results[op.ul.EVAL_KEYS.LLH_OE].append(output[op.ul.EVAL_KEYS.LLH_OE])
            

def benchmark(adata_list, 
              model_params):
    """Benchmarks methods on a list of adata
    
    Parameters
    ----------
    adata_list: list
        A list of AnnData objects, each describing a different scDNA-seq data set
    model_params: dictionary
        A dictionary where the keys are model names and the values are the function/parameters for the model
        
    Return
    -------
    dictionary
        A dictionary of simulation results, one for each key in op.ul.EVAL_KEYS
    """
    sim_results = {
        op.ul.EVAL_KEYS.MODEL: [],
        op.ul.EVAL_KEYS.PAIRWISE_REL_ACC: [],
        op.ul.EVAL_KEYS.MATRIX_ERROR: [],
        op.ul.EVAL_KEYS.LLH_BB: [],
        op.ul.EVAL_KEYS.LLH_OE: [],
        op.ul.SIM_KEYS.MUTATIONS: [],
        op.ul.SIM_KEYS.CELLS: [],
        op.ul.SIM_KEYS.FNR: [],
        op.ul.EVAL_KEYS.RUNTIME: [],
        op.ul.EVAL_KEYS.DATASET: [],
    }
    
    for i,adata in enumerate(adata_list):
        evaluate(adata, f"Dataset {i}", model_params, sim_results)

    return pd.DataFrame(sim_results)

def sim_benchmark(simulator,
                  sim_params, 
                  model_params, 
                  cells_geq_mutations=True, 
                  generate_seed=True, 
                  return_adata=False,
                  save_path=""):
    """Function to perform model benchmarking on simulation data
    
    Parameters
    ----------
    sim_params: dictionary
        A dictionary where keys are one of the SIM_KEYS and values are the values to use for those parameters
    model_params: dictionary
        A dictionary where the keys are model names and the values are the function/parameters for the model
    cells_geq_mutations: bool, optional
        Make sure that all simulations have at least as many cells as mutations, Default = True
    generate_seed: bool, optional
        Use a random seed when creating each data set. Default = True
    return_adata: bool, optional
        Returns a list of AnnData objects for each simulation
    save_path: str, optional
        A path to save each data sets to. The data for each simulation will be stored in its own folder.
        
    Return
    -------
    dictionary
        A dictionary of simulation results, one for each key in op.ul.EVAL_KEYS
    """
    sim_results = {
        op.ul.EVAL_KEYS.MODEL: [],
        op.ul.EVAL_KEYS.PAIRWISE_REL_ACC: [],
        op.ul.EVAL_KEYS.MATRIX_ERROR: [],
        op.ul.EVAL_KEYS.LLH_BB: [],
        op.ul.EVAL_KEYS.LLH_OE: [],
        op.ul.SIM_KEYS.MUTATIONS: [],
        op.ul.SIM_KEYS.CELLS: [],
        op.ul.SIM_KEYS.FNR: [],
        op.ul.EVAL_KEYS.RUNTIME: [],
        op.ul.EVAL_KEYS.DATASET: [],
    }
    
    res = load_sim_params(sim_params, 
                          cells_geq_mutations=cells_geq_mutations, 
                          generate_seed=generate_seed)
    
    if save_path != "" and not os.path.exists(save_path):
        os.mkdir(save_path)
        
    adata_list = []
    
    for setup in res:
        
        dataset_name = "sim%d_n%d_m%d" % (setup[op.ul.SIM_KEYS.TRIAL_NUM],
                                         setup[op.ul.SIM_KEYS.CELLS],
                                         setup[op.ul.SIM_KEYS.MUTATIONS])

        # create save path if it doesn't exist
        if save_path != "" and os.path.exists(save_path):
            os.mkdir(os.path.join(save_path, dataset))
            
        # simulate data set
        function_args = {}
        for k,v in setup.items():
            if k in simulator.__code__.co_varnames:
                function_args[k] = v

        adata = simulator(**function_args)

        if return_adata:
            adata_list.append(adata)
        
        # evaluate on data set
        evaluate(adata, dataset_name, model_params, sim_results)

    return pd.DataFrame(sim_results), adata_list
