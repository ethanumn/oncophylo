import numpy as np 
from itertools import combinations 

def find_loss_pairs(input_df, mutations, k):
    """Finds pairs of mutations that could have been lost together. Uses KModes clustering to cluster mutations
    and then returns all pairs of mutations that co-occur in the same cluster.
    
    Parameters
    -----------
    input_df: pd.DataFrame
        The character matrix
    mutations: list
        A list of mutations to cluster
    k: int
        The number of clusters for KModes to resolve


    Returns
    --------
    list
        A list of tuples, where each tuple is a pair of mutations that co-occur in the same cluster
    """
    
    from kmodes.kmodes import KModes

    if isinstance(mutations, list):
        mutations = np.array(mutations)
        
    # Apply K-Modes clustering
    km = KModes(n_clusters=k, init='Cao', n_init=5, verbose=0)
    labels = km.fit_predict(input_df[mutations].values.T)  # Transpose to cluster features
    clusters = [mutations[labels == num].tolist() for num in np.unique(labels)]
    pairs = []
    for c in clusters:
        for i,j in combinations(c, 2):
            pairs.append((i,j))
            
    return pairs