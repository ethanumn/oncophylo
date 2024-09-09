# constants.py
from dataclasses import dataclass, fields 
from itertools import product 

mutation_types = {
                    2:0, # 2 is a loss which corresponds to a 0
                    4:1, # 4 is recurrent, which corresponds to a 1
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
    MATRIX_ERROR: str = "Matrix Error"
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
    
        
    VARIANT_READS: str = "variant_reads"
    TOTAL_READS: str = "total_reads"
    VARIANT_READS_CORRUPT: str = "variant_reads_corrupt"
    TOTAL_READS_CORRUPT: str = "total_reads_corrupt"
    TOTAL_COPY_NUMBERS: str = "total_copy_numbers"
    MUTANT_COPY_NUMBERS: str = "mutant_copy_numbers"
    COPY_STATES: str = "copy_states"
    FPR: str = "false_positive_rate"
    FNR: str = "false_negative_rate"
    MISSING_RATE: str = "missing_rate"
    LLH: str = "log_likelihood"
    MATRIX_ERROR: str = "matrix_error"
    REQUIRES_CLUSTERS: str = "requires_clusters"
        
@dataclass
class SIM_KEYS:
    """Dataclass that defines all keys for simulation parameters"""
    MUTATIONS: str = "Mutations"
    CELLS: str = "Cells"
    CLUSTERS: str = "Clusters"
    MAX_LOSSES: str = "Max Losses"
    MAX_GAINS: str = "Max Gains"
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