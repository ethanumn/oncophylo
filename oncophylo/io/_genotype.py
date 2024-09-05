import os
import pandas as pd

def load(data_fn, 
         cell_names_fn="", 
         gene_names_fn="", 
         sep=",",
         skiprows=0,
         header=0,
         index_col=0):
    
    """Loads a mutation matrix from a file.
    
    Parameters
    ----------
    data_fn: str
        A path to a delimited data file where the rows are cells and the columns are mutations. The alphabet should
        be valid for whatever method the data will be passed to.
    cell_names_fn: str
        A path to a file containing the names of each cell in the data_fn. Each row should be a distinct cell name.
        There should be no header or index.
    gene_names_fn: str
        A path to a file containing the names of each gene in the data_fn. Each row should be a distinct gene name.
        There should be no header or index.
    skiprows: int, optional
        The number of rows to skip when reading data_fn. Default = 0
    header: int, optional
        The row in data_fn that designates the header. Default = None
    index_col: int, optional
        The column in data_fn that designates the index. Default = None
        
    Returns
    -------
    pd.DataFrame
        A data frame where the rows are cells and the columns are mutations. The index labels will be populated
        by those names listed in cell_names_fn and the column labels will be populated by those names in 
        gene_names_fn.
    """
    if os.path.isfile(data_fn) and data_fn != "":
        data_df = pd.read_csv(data_fn, sep=sep, skiprows=skiprows, header=header, index_col=index_col)
        
    if os.path.isfile(cell_names_fn) and cell_names_fn != "":
        cells_df = pd.read_csv(cell_names_fn, index_col=False, header=None)
        assert data_df.shape[0] == cells_df.shape[0], \
                "Data is different length than row labels: %d vs. %d" % (data_df.shape[0], cells_df.shape[0])
        data_df.index = cells_df.values.reshape(-1)
        
    if os.path.isfile(gene_names_fn) and gene_names_fn != "":
        genes_df = pd.read_csv(gene_names_fn, index_col=False, header=None)
        assert data_df.shape[1] == genes_df.shape[0], \
                "Data is different length than column labels: %d vs. %d" % (data_df.shape[1], genes_df.shape[0])
        data_df.columns = genes_df.values.reshape(-1)

    return data_df