import numpy as np
from kmodes.kmodes import KModes as KM
from sklearn.metrics import silhouette_score

def KModes(input_df, 
           columns, 
           criterion="silhouette",
           k=None,
           kMin=2, 
           kMax=8,
           n_init=5,
           return_best=True):
    """Clusters cells in the input_df using a subset of mutation
    
    Parameters
    -----------
    input_df: pd.DataFrame
        A character matrix
    columns: list
        A list of column names to use for clustering
    criterion: str
        The criterion used to determine which number of clusters to use. Default = "silhouette"
    k: int
        The number of clusters to find in the data. If None, a clustering will be searched for between
        sizes kMin and kMax. If defined, k clusters will be returned. This overrides kMin and kMax.
    kMin: int
        The minimum number of clusters to consider. Default = 2
    kMax: int
        The maximum number of clusters to consider. Default = 8
    kColumn: int
        Cluster the columns of the input so that mutations are grouped in the output matrix.
        This is used for visualization purpose to visualize the data after clustering. If 0, no column
        clustering will be performed.
    n_init: int
        The number of different initialization to allow the clustering method to perform. Default = 5
    return_best: bool
        Flag to return only the best clustering. If False, returns data for all clusterings. Default = True
        
    Returns
    -------
    list
        A list of clusterings. If k is defined, then it will only contain a single clustering. Otherwise,
        kMax - kMin clusterings will be returned.
    list
        A list of scores for each clustering. The scores are determined using the criterion provided.
    """

    # lists to return
    clusterings = []
    scores = []
    
    # extract data for clustering
    X = input_df[columns].values
    
    # if k is defined, overwrite kMin and kMax
    if isinstance(k, int):
        kMin = k
        kMax = k+1
                
    # cluster cells
    for k in range(kMin, kMax):
        km = KM(n_clusters=k, init='Cao', n_init=n_init)
        labels = km.fit_predict(X) 
        clusters = [np.flatnonzero(labels == num) for num in np.unique(labels)]
        clusterings.append(clusters)
        scores.append(silhouette_score(X, labels))
    
    # return only the best 
    if return_best:
        clusterings = [clusterings[np.argmax(scores)]]
        scores = [scores[np.argmax(scores)]]

    return clusterings, scores
        