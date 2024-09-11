from sklearn.neighbors import NearestNeighbors
import pandas as pd
import numpy as np 


def custom_distance(cell1, cell2, missing_value):
    """Calculate the distance between two cells, ignoring positions with a missing value."""
    mask = (cell1 != missing_value) & (cell2 != missing_value)  # Ignore positions with 3's in either cell
    if np.any(mask):  # Check if there are any valid positions to compare
        return np.sum(np.abs(cell1[mask] - cell2[mask]))  # Use Manhattan distance as an example
    else:
        return np.inf  # If no valid positions, return an infinite distance

def knn_impute(input_df, k=3, missing_value=3):
    """Impute the missing values in the data using K-nearest neighbors."""
    data = input_df.values
    n_samples, n_features = data.shape
    imputed_data = data.copy()
    
    knn_model = NearestNeighbors(n_neighbors=k, metric=lambda a, b: custom_distance(a, b, missing_value=missing_value))
    knn_model.fit(data)
    
    for i in range(n_samples):
        cell = data[i]
        
        # Find the k-nearest neighbors for the current cell
        distances, indices = knn_model.kneighbors([cell], n_neighbors=k)
        neighbors = data[indices[0]]
        
        # Impute missing values  in the current cell
        for j in range(n_features):
            if cell[j] == missing_value:  # If value is missing
                # Get the corresponding values from the k-nearest neighbors
                neighbor_values = neighbors[:, j]
                valid_values = neighbor_values[neighbor_values != missing_value]  # Ignore neighbors with missing values
                
                if len(valid_values) > 0:  # If there are valid neighbors
                    # Use majority vote for imputation
                    imputed_value = np.bincount(valid_values).argmax()
                    imputed_data[i, j] = imputed_value
    
    return pd.DataFrame(imputed_data, index=input_df.index, columns=input_df.columns)