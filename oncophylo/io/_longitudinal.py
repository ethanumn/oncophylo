import numpy as np
import pandas as pd

def convert_inputs(variants_df, regions_df, remove_snps=False):
    """
    Converts the input dataframes to the required format for scOrchard

    Parameters
    ----------  
    variants_df : DataFrame
        DataFrame containing the variant information (format COMPASS expects)
    regions_df : DataFrame
        DataFrame containing the number of reads that fell into each region for each cell (format COMPASS expects)
    remove_snps : bool, optional
        Flag to remove SNPs from the analysis. The default is False.

    Returns
    -------
    character_matrix : DataFrame
        Character matrix for the variants
    variant_reads_df : DataFrame
        DataFrame containing variant read counts
    total_reads_df : DataFrame
        DataFrame containing total read counts
    meta_df : DataFrame
        DataFrame containing metadata for the variants
    """

    # remove variants that are snps
    if remove_snps:
        variants_df = variants_df[~(variants_df["FREQ"] > 0.0)]
        
    # reformat data
    cells = variants_df.columns[6:]
    
    character_data = []
    variant_reads = []
    total_reads = []
    for col in cells:
        c = []
        v = []
        t = []
        for val in variants_df[col]:
            split = val.split(":")
            ref = int(split[0])
            alt = int(split[1])
            genotype = int(split[2])
            c.append(genotype)
            v.append(alt)
            t.append(alt+ref)
        character_data.append(c)
        variant_reads.append(v)
        total_reads.append(t)
        
    # remove strings from names that will throw errors
    mutations = variants_df["NAME"].apply(lambda x: x.replace(" ", "_").replace(".", "").replace("*", "")).values
    variants_df["NAME"] = mutations
        
    character_matrix = pd.DataFrame(character_data, index=cells, columns=mutations)
    variant_reads_df = pd.DataFrame(variant_reads, index=cells, columns=mutations)
    total_reads_df = pd.DataFrame(total_reads, index=cells, columns=mutations)
    meta_df = variants_df.iloc[:, :6].reset_index()
    
    return character_matrix, variant_reads_df, total_reads_df, meta_df

def preprocess_longitudinal(data, remove_snps=False):
    """
    Preprocesses data for longitudinal analysis 

    Parameters
    ----------
    data : list of tuples
        Each tuple contains two dataframes:
        - variants_df: DataFrame containing variant information
        - regions_df: DataFrame containing region information
    remove_snps : bool, optional
        Whether to remove SNPs from the analysis. The default is False.
        
    Returns
    -------
    merged_character_matrix : DataFrame
        Merged character matrix from all samples.
    merged_variant_reads_df : DataFrame
        Merged variant reads dataframe from all samples.
    merged_total_reads_df : DataFrame
        Merged total reads dataframe from all samples.
    merged_meta_df : DataFrame
        Merged metadata dataframe from all samples.
    merged_regions_df : DataFrame
        Merged regions dataframe from all samples.
    cell_samples : Series
        Series indicating the sample each cell belongs to.
    """
    cm_list = []
    vr_list = []
    tr_list = []
    regions_list = []
    cell_samples = []
    merged_meta_df = pd.DataFrame()
    for i,(variants_df, regions_df) in enumerate(data):
        # converting variants and regions dataframes to inputs
        character_matrix, variant_reads_df, total_reads_df, meta_df =  convert_inputs(variants_df, 
                                                                                      regions_df, 
                                                                                      remove_snps=remove_snps)
        cm_list.append(character_matrix)
        vr_list.append(variant_reads_df)
        tr_list.append(total_reads_df)
        regions_list.append(regions_df)
        
        # add cell sample label
        cell_samples += [i]*character_matrix.shape[0]
        
        # process meta data
        if merged_meta_df.shape[0] > 0:
            meta_df.loc["SAMPLE"] = np.nan
            mutations_in_sample = ~meta_df["NAME"].isin(merged_meta_df["NAME"])
            merged_meta_df = pd.merge(merged_meta_df, meta_df, how="outer").reset_index(drop=True)
            merged_meta_df["SAMPLE"] = merged_meta_df["SAMPLE"].fillna(i)
        else:
            merged_meta_df = meta_df.copy()
            merged_meta_df["SAMPLE"] = i
            
    merged_meta_df = merged_meta_df.dropna()
    merged_meta_df["SAMPLE"] = merged_meta_df["SAMPLE"].astype(int)
    merged_meta_df["CHR"] = merged_meta_df["CHR"].astype(int)
    merged_meta_df["POS"] = merged_meta_df["POS"].astype(int)

    # combine all matrices into one
    merged_character_matrix = pd.concat(cm_list).fillna(0).astype(int).reset_index(drop=True)
    merged_variant_reads_df = pd.concat(vr_list).fillna(0).astype(int).reset_index(drop=True)
    merged_total_reads_df = pd.concat(tr_list).fillna(0).astype(int).reset_index(drop=True)
    merged_regions_df = pd.concat(regions_list, axis=1).fillna(0).astype(int)
    merged_regions_df.columns = np.arange(0, len(cell_samples))
    
    # convert cell_samples to NumPy array
    cell_samples = pd.Series(cell_samples)
    
    return merged_character_matrix, merged_variant_reads_df, merged_total_reads_df, merged_meta_df, merged_regions_df, cell_samples
        