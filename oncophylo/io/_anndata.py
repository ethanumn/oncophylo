import networkx as nx 
import pydot 
from networkx.drawing.nx_pydot import to_pydot, from_pydot
from anndata import read_h5ad 

def write_adata(adata, file_name):
    """Write an AnnData to a file, while converting the Networkx.DiGraph's in uns to a 
    saveable format
    
    Parameters
    ----------
    adata: AnnData
        The AnnData object to be saved. 
    file_name: str
        The path to save the AnnData object to 
        
    Returns
    -------
    None
    """
    # convert all DiGraphs to compatible format
    keys = list(adata.uns.keys())
    for k in keys:
        v = adata.uns[k]
        if isinstance(v, nx.DiGraph):
            adata.uns[k] = {"tree": to_pydot(v).to_string(), "graph": v.graph}
            
    adata.write(file_name)

def read_adata(file_name):
    """Reads an AnnData from a file, while converting the data for trees in uns back
    into Networkx.DiGraph's
    
    Parameters
    ----------
    file_name: str
        The file to load the AnnData object from 
        
    Returns
    -------
    object
        The AnnData object 
    """
    adata = read_h5ad(file_name)
    keys = list(adata.uns.keys())
    for k in keys:
        v = adata.uns[k]
        if isinstance(v, dict):
            if "tree" in v and "graph" in v:
                pydot_graph = pydot.graph_from_dot_data(v["tree"])[0]
                G = from_pydot(pydot_graph)
                # Convert MultiDiGraph to DiGraph by keeping the first edge between each pair of nodes
                T = nx.DiGraph()
                for node1, node2, data in G.edges(data=True):
                    if not T.has_edge(node1, node2):
                        T.add_edge(node1, node2, **data)
                T.graph = v["graph"]
                adata.uns[k] = T
                
    return adata