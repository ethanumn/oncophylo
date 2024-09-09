import numpy as np 
   
def show_matrix(input_df, 
                columns, 
                clusters, 
                show_xlabels=False,
                show_ylabels=False,
                num_column_clusters=3, 
                cluster_kwargs={"n_init":5, "init":"Cao"}):
    """
    Shows a (binary) character matrix as a heatmap. Shows only the columns provided from
    the input_df, and groups the rows according to the clusters. The columns can also 
    be clustered into a specific number of groups.
    
    Parameters
    -----------
    input_df: pd.DataFrame
        A character matrix
    columns: list
        A list of column names to visualize
    clusters: list
        A list of lists where each sublist is the index of the rows that are in the same cluster
    num_column_clusters: int
        A number of clusters to group the columns into. Default = 3
    cluster_kwargs: dict
        The keyword arguments passed to KModes for clustering the columns.
    Returns
    --------
    None
    """    
    import seaborn as sns
    
    if isinstance(columns, list):
        columns = np.array(columns)

    if num_column_clusters > 1:
        
        from kmodes.kmodes import KModes as KM

        # cluster columns together
        X = input_df[columns].values.T
        km = KM(n_clusters=num_column_clusters, **cluster_kwargs)
        labels = km.fit_predict(X) 
        col_clusters = [np.where(labels == num)[0] for num in np.unique(labels)]

        # plot the clusters with the largest amount of data first
        cluster_order = np.argsort([np.max(X[:,c].sum(axis=0)) for c in col_clusters])[::-1]
        col_ordering = []

        for idx in cluster_order:
            c = col_clusters[idx]
            c_prime = c[np.argsort(X[:,c].sum(axis=0))[::-1]]
            col_ordering += c_prime.tolist()
    else:
        col_ordering = np.arange(len(columns))
            
    # plot the clusters with the largest amount of data first
    cluster_order = np.argsort([np.max(input_df.iloc[c].sum(axis=1)) for c in clusters])[::-1]
    row_ordering = []

    for idx in cluster_order:
        c = clusters[idx]
        c_prime = c[np.argsort(input_df.iloc[c].sum(axis=1))[::-1]]
        row_ordering += c_prime.tolist()
        
    columns = columns[col_ordering]
    rows = np.array(input_df.index[row_ordering])
    
    # show xlabels
    xlabels= columns if show_xlabels else []
    ylabels= rows if show_ylabels else []
        
    sns.heatmap(input_df.loc[rows, columns].values, xticklabels=xlabels, yticklabels=ylabels, cmap="Reds")
    