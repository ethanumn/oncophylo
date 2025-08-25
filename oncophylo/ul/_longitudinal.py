import numpy as np
import pandas as pd

def convert_inputs(variants_df, regions_df, remove_snps=False):
    """
    Converts inputs from variants and regions dataframe format into four dataframes.

    Parameters
    -----------
    variants_df: pd.DataFrame
        A dataframe describing the variant read counts from a single-cell amplicon sequencing experiment. The
        expected columns (IN ORDER!) are: "CHR", "POS", "REF", "ALT", "REGION" "NAME", "FREQ", followed by a column
        for each cell sequenced. 
    regions_df: pd.DataFrame
        A data frame describing the region coverage from a singel-cell amlpicon sequencing experiment. The columns
        are cells and the rows (indices) are regions. The region names should be in the format of "{chromosome}_{region}",
        where chromosome is a number and region is the name of the region (e.g., a gene name). The values in each entry
        should be the number of reads that mapped to that region in a cell.
    remove_snps: bool, optional
        A flag to remove single nucleotide polymorphisms (SNPs) from the input data. This should only be used if you do 
        not want to use SNPs to detect normal (i.e., non-cancerous) cells. Default is False.

    Returns
    --------
    pd.DataFrame
        A character matrix where rows are variants and columns are cells, and each entry is 
        one of the following: (1) heterozygous, (2) homozygous, (3 or -1) unknown, (0) absent.
    pd.DataFrame
        A variant read matrix where rows are variants, columns are cells, and each entry is the number of variant 
        reads mapped the variant containing locus in each cell. 
    pd.DataFrame
        A total read matrix where rows are variants, columns are cells, and each entry is the total number of
        reads mapped to the variant containing locus in each cell. 
    pd.DataFrame
        A dataframe of meta data for each variant containing the following columns: "CHR", "POS", "REF", "ALT", "REGION",
        "NAME", "FREQ", "SAMPLE".
    """
    if remove_snps:
        variants_df = variants_df[~(variants_df["FREQ"] > 0.0)]
        
    # reformat data, assumes first 6 columns are meta data
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

def preprocess_longitudinal(data, remove_snps=False, variants_to_remove=[]):
    """
    Processes a list of tuples containing data for each longitudinal sample e.g., [(sample1_variants_df, sample1_regions_df), ...]
    and returns dataframes with all of the data combined, along with annotations for which sample each cell came from.

    Parameters
    ----------
    data: list
        An ordered list of tuples, where each tuple consists of a variants_df and regions_df describing the
        cells sequenced in sample, and the samples are ordered by time, where the first tuple contains the
        data from the earliest longitudinal sample, and the last contains the data from the last longitudinal sample.
    remove_snps: bool, optional
        A flag to remove single nucleotide polymorphisms (SNPs) from the input data. This should only be used if you do 
        not want to use SNPs to detect normal (i.e., non-cancerous) cells. Default is False.
    variants_to_remove: list, optional
        A list of variants to remove from the input. Default is [].
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

    # make sure each column is in the same order as how the variants appear in the meta data
    column_order = merged_meta_df["NAME"].tolist()
    merged_character_matrix = merged_character_matrix[column_order]
    merged_variant_reads_df = merged_variant_reads_df[column_order]
    merged_total_reads_df = merged_total_reads_df[column_order]

    # remove variants
    merged_character_matrix = merged_character_matrix.drop(columns=variants_to_remove)
    merged_variant_reads_df = merged_variant_reads_df.drop(columns=variants_to_remove)
    merged_total_reads_df = merged_total_reads_df.drop(columns=variants_to_remove)
    merged_meta_df = merged_meta_df[~merged_meta_df["NAME"].isin(variants_to_remove)]
    
    # convert cell_samples to NumPy array
    cell_samples = pd.Series(cell_samples)
    
    return merged_character_matrix, merged_variant_reads_df, merged_total_reads_df, merged_meta_df, merged_regions_df, cell_samples