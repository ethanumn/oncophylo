import numpy as np 
import seaborn as sns
import matplotlib as mp 
import matplotlib.pyplot as plt

def show_matrix(input_df, 
                columns, 
                clusters, 
                show_xlabels=False,
                show_ylabels=False,
                num_column_clusters=3, 
                cluster_kwargs={"n_init":5, "init":"Cao"},
                figsize=(17,10),
                cmap="Reds",
                cluster_cmap="Blues",
                legend_bbox=(0.95,1),
                show_clusters=True,
                show_legend=False):
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
    value_map = {0: "Absent", 1: "Present", -1: "Unknown", 3: "Unknown", 2: "Lost", 4: "Recurrent"}

    discrete_colors = [mp.colormaps[cluster_cmap](i / (len(clusters) - 1)) for i in range(len(clusters))]
    cmap = sns.color_palette(cmap, as_cmap=True)

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
    row_colors = []
    
    for idx in cluster_order:
        c = clusters[idx]
        c_prime = c[np.argsort(input_df.iloc[c].sum(axis=1))[::-1]]
        row_ordering += c_prime.tolist()
        row_colors += [discrete_colors[idx % len(discrete_colors)]]*len(c_prime)
        
    columns = columns[col_ordering]
    rows = np.array(input_df.index[row_ordering])
    
    # show xlabels
    xlabels= columns if show_xlabels else []
    ylabels= rows if show_ylabels else []
    
    # make heatmap and cluster strip
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(input_df.loc[rows, columns].values, ax=ax, xticklabels=xlabels, yticklabels=ylabels, cmap="Reds", cbar=False)
    
    if show_clusters:
        # Add a color ribbon (rectangle patches) on the side
        for i, color in enumerate(row_colors):
            rect = mp.patches.Rectangle((input_df.shape[1] + 0.1, i), input_df.shape[1] / 30, 1, linewidth=0, edgecolor=None, facecolor=color)
            ax.add_patch(rect)

        plt.xlim([0, input_df.shape[1] + 2])
        
    if show_legend:
        # Custom legend
        norm = mp.colors.Normalize(vmin=np.min(input_df.values), vmax=np.max(input_df.values))
        values = np.unique(input_df.values)  # Values for which to extract colors

        colors = [cmap(norm(value)) for value in values]

        # Display the colors in a legend
        legend_elements = [mp.patches.Patch(color=color, label=f'{value_map[value]}') 
                           for value, color in zip(values, colors)]

        plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=legend_bbox, fontsize="large")
    plt.show()
    
    return fig, ax