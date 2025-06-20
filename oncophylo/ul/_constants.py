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
    LLH_OE: str = "Log Likelihood (Observation error model)"
    LLH_BB: str = "Log Likelihood (Beta-Binomial model)"
    CASET: str = "CASET"
    DISC: str = "DISC"
    RUNTIME: str = "Run time"
        
@dataclass
class DATA:
    CELL_TREE: str = "cell_tree"
    MUTATION_TREE: str = "mutation_tree"
    CLONAL_TREE: str = "clonal_tree"
    CN_TREE: str = "copy_number_tree"
    OBS_DATA: str = "observed_data"
    TRUE_DATA: str = "ground_truth_data"
    CLUSTER_ID: str = "cluster_id"
    CELL_SAMPLE: str = "cell_sample"
    PRED_DATA: str = "predicted_data"
    RUNTIME: str = "run_time"
    TERMINAL_OUTPUT: str = "terminal_output"
    CLUSTERS: str = "clusters"
    VARIANT_READS: str = "variant_reads"
    TOTAL_READS: str = "total_reads"
    VARIANT_READS_CORRUPT: str = "variant_reads_corrupt"
    TOTAL_READS_CORRUPT: str = "total_reads_corrupt"
    REGION_READS: str = "region_reads"
    TOTAL_COPY_NUMBERS: str = "total_copy_numbers"
    MUTANT_COPY_NUMBERS: str = "mutant_copy_numbers"
    COPY_STATES: str = "copy_states"
    FPR: str = "fp_rate"
    FNR: str = "fn_rate"
    MISSING_RATE: str = "missing_rate"
    REQUIRES_CLUSTERS: str = "requires_clusters"
    REGION_PROBABILITIES: str = "region_probabilities"
    DROPOUT_RATES: str = "dropout_rates"
        
@dataclass
class SIM_KEYS:
    """Dataclass that defines all keys for simulation parameters"""
    MUTATIONS: str = "num_mutations"
    CELLS: str = "num_cells"
    CLUSTERS: str = "num_clusters"
    MAX_LOSSES: str = "max_losses"
    MAX_GAINS: str = "max_gains"
    SEED: str = "seed"
    FPR: str = "fp_rate"
    FNR: str = "fn_rate"
    MISSING_RATE: str = "missing_rate"
    MEAN_COVERAGE: str = "mean_coverage"
    ADO_PRECISION: str = "ado_precision"
    VAF_THRESHOLD: str = "vaf_threshold"
    MUTATION_RATE: str = "mutation_rate"
    MAX_CN: str = "max_cn"
    TRIALS: str = "trials"
    TRIAL_NUM: str = "trial_number"
    PREFIX: str = "prefix"
    SAVE_PATH: str = "save_path"