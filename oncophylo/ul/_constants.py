# constants.py
from dataclasses import dataclass, fields 
from itertools import product 

mutation_types = {
                    2:0, # 2 is a loss which corresponds to a 0
                    4:1, # 4 is homoplasy, which corresponds to a 1
                    5:1  # 5 is recurrent, which corresponds to a 1
                 }

@dataclass
class EVAL_KEYS:
    """Dataclass that defines the keys for different evaluation metrics"""
    
    MODEL: str = "Model"
    DATASET: str = "Data set"
    AD_RECALL: str = "Ancestor-Descendant Recall"
    AD_PRECISION: str = "Ancestor-Descendant Precision"
    DL_RECALL: str = "Different-Lineage Recall"
    DL_PRECISION: str = "Different-Lineage Precision"
    CC_RECALL: str = "Co-cluster Recall"
    PAIRWISE_REL_ACC: str = "Pairwise Relationship Accuracy"
    CASET: str = "CASET"
    DISC: str = "DISC"
    RUNTIME: str = "Run time"
        
@dataclass
class CONST:
    CELL_TREE: str = "cell_tree"
    MUTATION_TREE: str = "mutation_tree"
    CLONAL_TREE: str = "clonal_tree"
    CN_TREE: str = "copy_number_tree"
    OBS_DATA: str = "observed_data"
    TRUE_DATA: str = "ground_truth_data"
    CLUSTER_ID: str = "cluster_id"
    PRED_DATA: str = "predicted_data"
    RUNTIME: str = "run_time"
    TERMINAL_OUTPUT: str = "terminal_output"
    CLUSTERS: str = "clusters"
    
        
    NOISE_FREE_COUNTS: str = "noise_free_read_counts"
    NOISY_COUNTS: str = "noisy_read_counts"
    VARIANT: str = "variant"
    TOTAL: str = "total"
    FPR: str = "false_positive_rate"
    FNR: str = "false_negative_rate"
    LLH: str = "log_likelihood"
    MATRIX_ERROR: str = "matrix_error"
        
@dataclass
class SIM_KEYS:
    """Dataclass that defines all keys for simulation parameters"""
    MUTATIONS: str = "Mutations"
    CELLS: str = "Cells"
    CLUSTERS: str = "Clusters"
    MAX_LOSSES: str = "Max Losses"
    SEED: str = "Seed"
    FPR: str = "false_positive_rate"
    FNR: str = "false_negative_rate"
    MISSING_RATE: str = "Missing Rate"
    MEAN_COVERAGE: str = "Mean Coverage"
    ADO_PRECISION: str = "ADO Precision"
    VAF_THRESHOLD: str = "VAF Threshold"
    MUTATION_RATE: str = "Mutation Rate"
    MAX_CN: str = "Max Copy Number"
    RETURN_CNTREE: str = "Return CN Tree"
    RETURN_READS: str = "Return Reads"
    TRIALS: str = "Trials"
    TRIAL_NUM: str = "Trial Number"
    PREFIX: str = "Prefix"
    SAVE_PATH: str = "Save Path"