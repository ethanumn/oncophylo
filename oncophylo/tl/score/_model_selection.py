import networkx as nx
import numpy as np 
from scipy.stats import binomtest

def log10_BF(MI, Mi, n, m):
    """
    Computes the base-10 logarithm Bayes Factor to choose between 
    two models as described here: https://genome.cshlp.org/content/27/11/1885.long
    Note that the definitions of m and n are switched compared to the referenced paper.

    Input
    ------
    MI : float
        A base-10 log score from the the baseline model
    Mi : float
        A base-10 log score from the alternative model
    n : int
        The number of cells in the input data used to generate trees from each model
    m : int 
        The number of mutations in the input data used to generate trees from each model

    """
    # CRITICAL! if we don't perform (n(n-1))/2 separately, it'll try to convert everything a float
    mKi = int(m*(m-1)/2)*(m+2)**(m-1+n) 
    KI = (m+1)**(m+n-1) # Eq 19
    KF = mKi - m*(m-1)*KI # Eq 22
    
    return np.log10((KI/KF)*(10**(Mi - MI)))

def loss_test(T, input_df, mutation, fn):
    """Test whether or not a mutation is missing in descendant cells significantly more than implied by 
    the provided false negative rate. The null hypothesis is that all relevant cells missing the mutation
    can be explained by the provided false negative rate. The alternative hypothesis is that the number of cells
    with the missing mutation is significantly higher than can be explained by the false negative rate. 
    The binomial test is used to determine significance.
    
    Parameters
    ----------
    T: Networkx.DiGraph
        A cell tree where internal nodes are mutations and leaves are cells. It's assume that the tree has
        a graph attribute with the key 'cells' 
    input_df: pd.DataFrame  
        The character matrix from which the tree was inferred
    mutation: str
        The name of the mutation to test whether or not it was lost.
    fn: float
        The false negative rate

    Returns
    --------
    float
        The p-value of whether or not to reject the null hypothesis
    int
        The number of cells that are descendants of the mutation in the tree
    int
        The number of cells that are descendants of the mutation in the tree, and the mutation is missing from input_df
    """
    # get all cells attached below mutation
    cells = list(set(nx.descendants(T, mutation)).intersection(set(T.graph["cells"])))
    
    n = len(cells)
    k = np.sum(input_df.loc[cells, mutation] == 0)
        
    return binomtest(k, n, p=fn, alternative='greater').pvalue, n, k
    
def recurrence_test(T, input_df, mutation, fp):
    """
    Test whether or cells in the character matrix contain the mutation, but their inferred genotypes do not. 
    This tests whether there are more false positives than can be explained by the false positive rate. 
    The null hypothesis is that all relevant cells with a false positive can be explained by the false positive rate. 
    The alternative hypothesis is that the number of false positives is significantly higher 
    than can be explained by the false positive rate. The binomial test is used to determine significance.
    
    Input
    ------
    T: Networkx.DiGraph
        A cell tree where internal nodes are mutations and leaves are cells. It's assume that the tree has
        a graph attribute with the key 'cells' 
    input_df: pd.DataFrame  
        The character matrix from which the tree was inferred
    mutation: str
        The name of the mutation to test whether or not it was lost.
    fn: float
        The false negative rate

    Returns
    --------
    float
        The p-value of whether or not to reject the null hypothesis
    int
        The number of cells that are descendants of the mutation in the tree
    int
        The number of cells that are descendants of the mutation in the tree, and the mutation is missing from input_df
    """
    # collect cells that do not have the mutation as an ancestor
    cells = []  
    
    # get all cells attached below mutation
    for c in T.graph["cells"]:
        if nx.lowest_common_ancestor(T, c, mutation) != mutation:
            cells.append(c)
            
    n = len(cells)
    k = np.sum(input_df.loc[cells, mutation] == 1)
        
    return binomtest(k, n, p=fp, alternative='greater').pvalue, n, k
    