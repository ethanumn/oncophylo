import pandas as pd 
import networkx as nx 
import numpy as np 
from itertools import combinations 

import oncophylo as op 
from oncophylo.ul import CONST

def score_genotypes(B, B_input, fp, fn):
    """Scores a binary genotype matrix given the groun truth genotypes, 
    false positive rate, and false negative rate"""
    _B = B
    _B_input = B_input
    
    # process inputs
    if isinstance(B, pd.DataFrame):
        if CONST.CLUSTER_ID in list(_B.columns):
            _B = _B.drop(columns=CONST.CLUSTER_ID)
        _B = B.values

    if isinstance(_B_input, pd.DataFrame):
        if CONST.CLUSTER_ID in list(_B_input.columns):
            _B_input = _B_input.drop(columns=CONST.CLUSTER_ID)
        _B_input = _B_input.values
        

    assert _B.shape == _B_input.shape, \
                "Genotype matrices must be same shape! Received (%d,%d), expected (%d,%d)." % (_B.shape[0], 
                                                                                               _B.shape[1], 
                                                                                               _B_input.shape[0], 
                                                                                               _B_input.shape[1])
    llh = np.sum((_B == 1) & (_B_input == 1)) * np.log(1-fn) + \
          np.sum((_B == 0) & (_B_input == 0)) * np.log(1-fp) + \
          np.sum((_B == 1) & (_B_input == 0)) * np.log(fn) + \
          np.sum((_B == 0) & (_B_input == 1)) * np.log(fp)
    
    return llh

    from itertools import combinations

def matrix_error(B, B_input):
    """Computes the matrix error metric  
    
    (Difference between observed in B and B_input) / (# observed in both B and B_input)
    """ 
    _B = B
    _B_input = B_input
    
    # process inputs
    if isinstance(_B, pd.DataFrame):
        if CONST.CLUSTER_ID in list(_B.columns):
            _B = _B.drop(columns=CONST.CLUSTER_ID)
        _B = _B.values

    if isinstance(_B_input, pd.DataFrame):
        if CONST.CLUSTER_ID in list(_B_input.columns):
            _B_input = _B_input.drop(columns=CONST.CLUSTER_ID)
        _B_input = _B_input.values
    
    diff = np.abs(_B - _B_input)
    observed = ((_B == 0) | (_B == 1)) & ((_B_input == 0) | (_B_input == 1))
    
    return np.sum(diff[observed]) / np.sum(observed)

def pairwise_rel_accuracy(T, T_true):
    """Pairwise relationship accuracy"""
    assert T.graph["type"] == CONST.CLONAL_TREE and T.graph["type"] == CONST.CLONAL_TREE, "Input tree must be of type %s" % CONST.CLONAL_TREE

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
    assert T.graph["type"] == CONST.CLONAL_TREE, "Input type of tree must be clonal! Use to_clonal_tree() to convert a mutation tree to a cell tree."
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
    assert T.graph["type"] == CONST.CLONAL_TREE, "Input type of tree must have nodes as cells and edges as mutations! Use to_clonal_tree() to convert a mutation tree to a cell tree."
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
    assert T.graph["type"] == CONST.CLONAL_TREE, "Input type of tree must have nodes as cells and edges as mutations! Use to_clonal_tree() to convert a mutation tree to a cell tree."
    
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
    assert T.graph["type"] == CONST.MUTATION_TREE, "Input type of tree must be mutation!"
    ad = set()
    for ancestor in T.nodes():
        for descendant in nx.descendants(t, ancestor):
            if ancestor in T.graph["mutations"] and descendant in T.graph["mutations"]:
                ad.add((ancestor,descendant))
    return ad

def mt_dl(T):
    """Get all mutations on different lineages from a mutation tree"""
    assert T.graph["type"] == CONST.MUTATION_TREE, "Input type of tree must be mutation!"

    dl = set(frozenset((n1,n2)) for n1,n2 in combinations(T.graph["mutations"],2))
    for ancestor in T.nodes():
        for descendant in nx.descendants(T, ancestor):
            if ancestor in T.graph["mutations"] and descendant in T.graph["mutations"]:
                dl.remove(frozenset((ancestor,descendant)))
    return dl 

def get_ad(T):
    """Determines which function to call to get Ancestral-Descendant relationships"""
    if T.graph["type"] == CONST.MUTATION_TREE:
        ad = mt_ad(t)
    elif T.graph["type"] == CONST.CLONAL_TREE:
        ad = clonal_ad(T)
    else:
        raise Exception("Unable to determine type of input tree! Expects either %s or %s." % (CONST.MUTATION_TREE, CONST.CLONAL_TREE))
    return ad
        
def get_dl(T):
    """Determines which function to call to get Different Lineage relationships"""
    if T.graph["type"] == CONST.MUTATION_TREE:
        dl = mt_dl(T)
    elif T.graph["type"] == CONST.CLONAL_TREE:
        dl = clonal_dl(T)
    else:
        raise Exception("Unable to determine type of input tree! Expects either %s or %s." % (CONST.MUTATION_TREE, CONST.CLONAL_TREE))
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
    assert T.graph["type"] == CONST.CLONAL_TREE, "Input tree must be of type %s" % CONST.CLONAL_TREE

    cl = set()
    for c in clonal_clusters(T):
        if len(c) > 1:
            for p in combinations(c, 2):
                cl.add(frozenset(p))
                
    return cl

def cocluster_recall(T, T_true):
    """Co-cluster recall"""
    assert T.graph["type"] == CONST.CLONAL_TREE and T.graph["type"] == CONST.CLONAL_TREE, "Input tree must be of type %s" % CONST.CLONAL_TREE
    cl_pred = get_coclusters(T)
    cl_true = get_coclusters(T_true)
                
    if len(cl_true) == 0:
        return 1.0
    else:
        return len(cl_pred.intersection(cl_true)) / len(cl_true)