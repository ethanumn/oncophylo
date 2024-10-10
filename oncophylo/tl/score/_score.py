import pandas as pd 
import networkx as nx 
import numpy as np 
from itertools import combinations 
from scipy.stats import betabinom

import oncophylo as op 
from oncophylo.ul import DATA

def _score_beta_binomial(B, var_reads, total_reads, ado_precision, fp):
    """Scores a genotype matrix under a beta binomial model
    
    Parameters
    ----------
    B: numpy.ndarray
        A predicted binary matrix where rows are cells and columns are mutations.
    var_reads: pandas.DataFrame, optional
        A cell by mutation matrix of variant read counts
    total_reads: pandas.DataFrame, optional
        A cell by mutation matrix of total read counts
    ado_precision: float, optional
        The allelic dropout precision parameter used to calculate the beta binomial likelihood
    fp: float
        The estimated false positive rate

    Returns
    --------
    float
        The log likelihood of the predicted matrix B under the beta binomial model
    """
    llh = 0
    
    alpha_present = 1
    beta_present = 1
    alpha_absent = fp*ado_precision
    beta_absent = (1-fp)*ado_precision 
    n,m = B.shape
    for i in range(n):
        for j in range(m):
            if B[i,j] == 1:
                llh += betabinom.logpmf(var_reads[i,j], total_reads[i,j], alpha_present, beta_present)
            elif B[i,j] == 0:
                llh += betabinom.logpmf(var_reads[i,j], total_reads[i,j], alpha_absent, beta_absent)
    return llh

def _score_observation_errors(B, B_input, fp, fn):
    """Scores a genotype matrix based on the observed genotypes and the estimate false positive and false negative rates
    
    Parameters
    ----------
    B: numpy.ndarray
        A predicted binary matrix where rows are cells and columns are mutations.
    B_input: np.ndarray
        An observed binary matrix where rows are cells and columns are mutations. This matrix can have unknown or missing values (e.g., -1 or 3)
    fp: float
        The estimated false positive rate
    fn: float
        The estimated false negative rate

    Returns
    --------
    float
        The log likelihood of the predicted matrix B under the estimated observation error rates
    """

    llh = np.sum((B == 1) & (B_input == 1)) * np.log(1-fn) + \
          np.sum((B == 0) & (B_input == 0)) * np.log(1-fp) + \
          np.sum((B == 1) & (B_input == 0)) * np.log(fn) + \
          np.sum((B == 0) & (B_input == 1)) * np.log(fp)

    return llh


def score_observation_errors(B, B_input, fp, fn):
    """  
    Scores a genotype matrix using either observational error model from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0936-x

    Parameters
    ----------
    B: numpy.ndarray
        A predicted binary matrix where rows are cells and columns are mutations.
    B_input: np.ndarray
        An observed binary matrix where rows are cells and columns are mutations. This matrix can have unknown or missing values (e.g., -1 or 3)
    fp: float
        The estimated false positive rate
    fn: float
        The estimated false negative rate

    Returns
    --------
    float
        The log likelihood of the predicted genotypes
    """

    _B = B
    _B_input = B_input
    
    # process inputs
    if isinstance(_B, pd.DataFrame):
        if DATA.CLUSTER_ID in list(_B.columns):
            _B = _B.drop(columns=DATA.CLUSTER_ID)
        B_values = _B.values
    elif isinstance(_B, np.ndarray):
        B_values = _B
    else:
        print("Cannot compute beta binomial likelihood of genotypes. Predicted genotypes must be either a pandas.DataFrame or numpy.ndarray. Aborting!")
        return 

    if isinstance(_B_input, pd.DataFrame):
        if DATA.CLUSTER_ID in list(_B_input.columns):
            _B_input = _B_input.drop(columns=DATA.CLUSTER_ID)
        B_input_values = _B_input.values
    elif isinstance(_B_input, np.ndarray):
        B_input_values = _B_input
    else:
        print("Cannot compute beta binomial likelihood of genotypes. Observed genotypes must be either a pandas.DataFrame or numpy.ndarray. Aborting!")
        return 

    assert _B.shape == _B_input.shape, \
                "Genotype matrices must be same shape! Received (%d,%d), expected (%d,%d)." % (_B.shape[0], 
                                                                                               _B.shape[1], 
                                                                                               _B_input.shape[0], 
                                                                                               _B_input.shape[1])

    return _score_observation_errors(B_values, B_input_values, fp, fn)


def score_beta_binomial(B, var_reads, total_reads, ado_precision, fp):
    """  
    Scores a genotype matrix using the beta binomial likelihood from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03106-5

    Parameters
    ----------
    B: numpy.ndarray
        A predicted binary matrix where rows are cells and columns are mutations.
    var_reads: pandas.DataFrame, optional
        A cell by mutation matrix of variant read counts
    total_reads: pandas.DataFrame, optional
        A cell by mutation matrix of total read counts
    ado_precision: float, optional
        The allelic dropout precision parameter used to calculate the beta binomial likelihood
    fp: float

    Returns
    --------
    float
        The log likelihood of the predicted genotypes
    """

    # base case
    if (var_reads is None) or (total_reads is None) or (ado_precision is None):
        return np.nan 

    _B = B
    
    # process inputs
    if isinstance(_B, pd.DataFrame):
        if DATA.CLUSTER_ID in list(_B.columns):
            _B = _B.drop(columns=DATA.CLUSTER_ID)
        B_values = _B.values
    elif isinstance(B, np.ndarray):
        B_values = _B
    else:
        print("Cannot compute beta binomial likelihood of genotypes. Predicted genotypes must be either a pandas.DataFrame or numpy.ndarray. Aborting!")
        return 

    assert isinstance(ado_precision, float) or isinstance(ado_precision, int), "ado_precision must be a float or int, not %s" % type(ado_precision) 
    assert ado_precision > 0, "ado_precision must be non-negative"
    assert var_reads.shape == total_reads.shape, "var_reads and total_reads shapes do not match!"

    if isinstance(var_reads, pd.DataFrame) and isinstance(total_reads, pd.DataFrame):
        assert np.all(var_reads.columns == total_reads.columns) and np.all(var_reads.index == total_reads.index), "Columns and indices of var_reads do not match total_reads!"
        assert np.all(var_reads.columns == _B.columns) and np.all(var_reads.index == _B.index),  "Columns and indices of var_reads/total_reads do not match the character matrix!"
        var_reads_values = var_reads.values 
        total_reads_values = total_reads.values 
    elif isinstance(var_reads, np.ndarray) and isinstance(total_reads, pd.ndarray):
        var_reads_values = var_reads 
        total_reads_values = total_reads
    else:
        print("Types of var_reads and total_reads must both be either pandas.DataFrame or numpy.ndarray. Aborting!")
        return np.nan

    return _score_beta_binomial(B_values, var_reads_values, total_reads_values, ado_precision, fp)


def matrix_error(B, B_input):
    """Computes the matrix error metric  
    
    (Difference between observed in B and B_input) / (# observed in both B and B_input)
    """ 
    _B = B
    _B_input = B_input
    
    # process inputs
    if isinstance(_B, pd.DataFrame):
        if DATA.CLUSTER_ID in list(_B.columns):
            _B = _B.drop(columns=DATA.CLUSTER_ID)
        B_values = _B.values

    if isinstance(_B_input, pd.DataFrame):
        if DATA.CLUSTER_ID in list(_B_input.columns):
            _B_input = _B_input.drop(columns=DATA.CLUSTER_ID)
        B_input_values = _B_input.values
    
    diff = np.abs(B_values - B_input_values)
    observed = ((B_values == 0) | (B_values == 1)) & ((B_input_values == 0) | (B_input_values == 1))
    
    return np.sum(diff[observed]) / np.sum(observed)

def pairwise_rel_accuracy(T, T_true):
    """Pairwise relationship accuracy"""
    assert T.graph["type"] == DATA.CLONAL_TREE and T_true.graph["type"] == DATA.CLONAL_TREE, "Input tree must be of type %s" % DATA.CLONAL_TREE

    # get ancestral-descendant pairs
    ad_gt = get_ad(T_true)
    ad_pred = get_ad(T)
    
    # get different lineage pairs
    dl_gt = get_dl(T_true)
    dl_pred = get_dl(T)
    
    # get cocluster pairs
    cl_pred = get_coclusters(T)
    cl_true = get_coclusters(T_true)
    
    correct_pairs = len(cl_pred.intersection(cl_true)) + len(dl_pred.intersection(dl_gt)) + len(ad_pred.intersection(ad_gt))
    total_pairs = len(ad_gt) + len(dl_gt) + len(cl_true)
    return correct_pairs / total_pairs

def clonal_ad(T):
    """Computes the set of Ancestor-Descendant relationships in a clonal tree. Assumes nodes are cells and edges
    are labeled by mutations"""
    assert T.graph["type"] == DATA.CLONAL_TREE, "Input type of tree must be clonal! Use to_clonal_tree() to convert a mutation tree to a cell tree."
    root = op.ul.root_id(T)
    leaf_nodes = [node for node in T.nodes() if T.in_degree(node)!=0 and T.out_degree(node)==0]
    ad = set()
    for l in leaf_nodes:
        ancestral_mutations = []
        for path in nx.all_simple_paths(T, source=root, target=l):
            for i in range(1, len(path)):
                e = (path[i-1], path[i])
                muts_to_add = []
                for dm in T.edges[e]["label"].split(T.graph["splitter_mut"]):
                    for am in ancestral_mutations:
                        if am in T.graph["mutations"] and dm in T.graph["mutations"]:
                            ad.add((am, dm))
                    if dm in T.graph["mutations"]:
                        muts_to_add.append(dm)
                ancestral_mutations += muts_to_add
    return ad

def clonal_clusters(T):
    """Computes the set of clusters in a clonal tree"""
    assert T.graph["type"] == DATA.CLONAL_TREE, "Input type of tree must have nodes as cells and edges as mutations! Use to_clonal_tree() to convert a mutation tree to a cell tree."
    clusters = []
    root = op.ul.root_id(T)
    edge_dfs = nx.dfs_edges(T, source=root)
    ancestral_mutations = []
    for e in edge_dfs:
        c = T.edges[e]["label"].split(T.graph["splitter_mut"])
        if c in T.graph["mutations"]:
            clusters.append(c)
    return clusters

def clonal_dl(T):
    """Computes the set of Different-Lineage relationships in a clonal tree. Assumes the input is a clonal tree."""
    assert T.graph["type"] == DATA.CLONAL_TREE, "Input type of tree must have nodes as cells and edges as mutations! Use to_clonal_tree() to convert a mutation tree to a cell tree."
    
    # some mutations may be excluded and become germline, don't include those 
    mutations_in_tree = list(set(T.graph["mutations"]).difference(T.graph["become_germline"]))
    ad = set(frozenset(p) for p in clonal_ad(T))
    dl = set(frozenset(p) for p in combinations(mutations_in_tree,2))
    cl = set()
    for c in clonal_clusters(T):
        if len(c) > 1:
            for p in combinations(c, 2):
                cl.add(frozenset(p))
    return dl.difference(ad).difference(cl)

def mt_ad(T):
    """Gets all ancestor-descendant relationships in a mutation tree"""
    assert T.graph["type"] == DATA.MUTATION_TREE, "Input type of tree must be mutation!"
    ad = set()
    for ancestor in T.nodes():
        for descendant in nx.descendants(t, ancestor):
            if ancestor in T.graph["mutations"] and descendant in T.graph["mutations"]:
                ad.add((ancestor,descendant))
    return ad

def mt_dl(T):
    """Get all mutations on different lineages from a mutation tree"""
    assert T.graph["type"] == DATA.MUTATION_TREE, "Input type of tree must be mutation!"

    dl = set(frozenset((n1,n2)) for n1,n2 in combinations(T.graph["mutations"],2))
    for ancestor in T.nodes():
        for descendant in nx.descendants(T, ancestor):
            if ancestor in T.graph["mutations"] and descendant in T.graph["mutations"]:
                dl.remove(frozenset((ancestor,descendant)))
    return dl 

def get_ad(T):
    """Determines which function to call to get Ancestral-Descendant relationships"""
    if T.graph["type"] == DATA.MUTATION_TREE:
        ad = mt_ad(t)
    elif T.graph["type"] == DATA.CLONAL_TREE:
        ad = clonal_ad(T)
    else:
        raise Exception("Unable to determine type of input tree! Expects either %s or %s." % (DATA.MUTATION_TREE, DATA.CLONAL_TREE))
    return ad
        
def get_dl(T):
    """Determines which function to call to get Different Lineage relationships"""
    if T.graph["type"] == DATA.MUTATION_TREE:
        dl = mt_dl(T)
    elif T.graph["type"] == DATA.CLONAL_TREE:
        dl = clonal_dl(T)
    else:
        raise Exception("Unable to determine type of input tree! Expects either %s or %s." % (DATA.MUTATION_TREE, DATA.CLONAL_TREE))
    return dl

def ad_recall(T, T_true):
    """Ancestor-Descendant Recall"""
    # get ground truth ancestor-descendants
    ad_gt = get_ad(T_true)
    ad_pred = get_ad(T)
    if len(ad_pred) == 0 and len(ad_gt) == 0:
        return 1.0
    return len(ad_gt.intersection(ad_pred)) / np.maximum(1,len(ad_gt))
        
def dl_recall(T, T_true):
    """Different Lineage Recall"""
    # get ground truth ancestor-descendants
    dl_gt = get_dl(T_true)
    dl_pred = get_dl(T)
    if len(dl_pred) == 0 and len(dl_gt) == 0:
        return 1.0
    return len(dl_gt.intersection(dl_pred)) / np.maximum(1,len(dl_gt))

def get_coclusters(T):
    """Gets all pairs of mutations that are coclustered"""
    assert T.graph["type"] == DATA.CLONAL_TREE, "Input tree must be of type %s" % DATA.CLONAL_TREE

    cl = set()
    for c in clonal_clusters(T):
        if len(c) > 1:
            for p in combinations(c, 2):
                cl.add(frozenset(p))
                
    return cl

def cocluster_recall(T, T_true):
    """Co-cluster recall"""
    assert T.graph["type"] == DATA.CLONAL_TREE and T.graph["type"] == DATA.CLONAL_TREE, "Input tree must be of type %s" % DATA.CLONAL_TREE
    cl_pred = get_coclusters(T)
    cl_true = get_coclusters(T_true)
                
    if len(cl_true) == 0:
        return 1.0
    else:
        return len(cl_pred.intersection(cl_true)) / len(cl_true)