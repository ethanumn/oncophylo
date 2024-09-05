# CASet.py

def caset(mutations, t1_ancestors, t2_ancestors):
    """
    Compute the CASet distance between two trees stored as ancestor sets.
    :param mutations: the set of mutations over which to sum
    :param t1_ancestors: a dict storing the t1 ancestor sets of every mutation
    :param t2_ancestors: a dict storing the t2 ancestor sets of every mutation
    :return: the CASet distance between t1 and t2
    """
    m = len(mutations)
    total = 0
    for i in range(m):
        for j in range(i + 1, m):
            t1_ca_set = t1_ancestors[mutations[i]] & t1_ancestors[mutations[j]]
            t2_ca_set = t2_ancestors[mutations[i]] & t2_ancestors[mutations[j]]

            total += jaccard(t1_ca_set, t2_ca_set)

    return total / (m * (m - 1) / 2)


def caset_intersection(t1, t2):
    """
    Compute the CASet distance between two trees, only summing over pairs of mutations in both trees.
    :param t1: a Newick string representation of tree 1
    :param t2: a Newick string representation of tree 2
    :return: the CASet intersection distance between tree 1 and tree 2
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    intersection = list(set(t1_ancestors.keys()) & set(t2_ancestors.keys()))

    return caset(intersection, t1_ancestors, t2_ancestors)


def caset_union(t1, t2):
    """
    Compute the CASet distance between two trees, summing over all pairs of mutations.
    :param t1: a Newick string representation of tree 1
    :param t2: a Newick string representation of tree 2
    :return: the CASet union distance between tree 1 and tree 2
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    union = list(set(t1_ancestors.keys()) | set(t2_ancestors.keys()))

    for u in union:
        if u not in t1_ancestors:
            t1_ancestors[u] = set()
        if u not in t2_ancestors:
            t2_ancestors[u] = set()

    return caset(union, t1_ancestors, t2_ancestors)


def CASet(t1, t2, choice="intersection"):
    """CASet distance metric (https://academic.oup.com/bioinformatics/article/36/7/2090/5637226)"""
    # check and convert trees if they're networkx DiGraphs
    if isinstance(t1, networkx.classes.digraph.DiGraph):
        t1 = tree_to_newick(t1)
    if isinstance(t2, networkx.classes.digraph.DiGraph):
        t2 = tree_to_newick(t2)
        
    if choice == "intersection":
        return caset_intersection(t1, t2)
    elif choice == "union":
        return caset_union(t1, t2)
