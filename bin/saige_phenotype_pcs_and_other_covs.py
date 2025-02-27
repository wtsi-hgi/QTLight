#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2024-09-05'
__version__ = '0.0.1'

import pandas as pd
import numpy as np
import scanpy as sc
from kneed import KneeLocator
import sys
import scipy


def PF(X):
    cd=np.asarray(X.sum(1)).ravel()
    avg_cd=cd.mean()
    
    return scipy.sparse.diags(avg_cd/cd).dot(X)

def log1p(X):
    return X.log1p()

def quantile_normalize_vector(x):
    """
    Perform quantile normalization on a vector with random tie-breaking.
    """
    # Add random noise to break ties, scaled down to be very small
    random_noise = np.random.uniform(low=-1e-10, high=1e-10, size=len(x))
    ranks = scipy.stats.rankdata(x + random_noise, method='average')  # Add noise and rank
    return scipy.stats.norm.ppf(ranks / (len(x) + 1))

def quantile_normalize_matrix(matrix):
    """
    Perform quantile normalization on each column of a pandas DataFrame or NumPy array.
    Returns a DataFrame if input is a DataFrame.
    """
    is_dataframe = isinstance(matrix, pd.DataFrame)
    
    if is_dataframe:
        matrix_values = matrix.values  # Extract underlying NumPy array
    else:
        matrix_values = matrix
    
    quantile_matrix = np.zeros_like(matrix_values)
    for i in range(matrix_values.shape[1]):  # Iterate over columns
        quantile_matrix[:, i] = quantile_normalize_vector(matrix_values[:, i])
    
    if is_dataframe:
        return pd.DataFrame(quantile_matrix, index=matrix.index, columns=matrix.columns)
    else:
        return quantile_matrix

def quantile_normalize_rows(matrix):
    """
    Perform quantile normalization across rows for a pandas DataFrame or NumPy array.
    Returns a DataFrame if input is a DataFrame.
    """
    is_dataframe = isinstance(matrix, pd.DataFrame)
    
    if is_dataframe:
        matrix_values = matrix.values  # Extract underlying NumPy array
    else:
        matrix_values = matrix
    
    quantile_matrix = quantile_normalize_matrix(matrix_values.T).T  # Transpose for row normalization
    
    if is_dataframe:
        return pd.DataFrame(quantile_matrix, index=matrix.index, columns=matrix.columns)
    else:
        return quantile_matrix



def run_pca_and_save(file_path, output_file, n_pcs, covs):
    # Load the data from the flat text file
    counts = pd.read_csv(file_path, sep="\t", header=0, index_col=0)
    with open(covs, 'r') as file:
        first_line = file.readline().strip()
    existing_covs = first_line.split(',')
    genotype_pcs = counts[existing_covs]
    counts = counts.drop(columns=existing_covs, errors='ignore')
    # 
    # Convert to AnnData object
    adata = sc.AnnData(X=counts)
    pheno_id = counts.index.name
    adata.var_names = counts.columns
    adata.obs_names = counts.index
    
    counts_orig = counts.copy()
    # counts_orig = pd.DataFrame(PF(counts), index=counts.index, columns= counts.columns)
    # counts_orig= pd.DataFrame(sc.pp.log1p(sc.pp.normalize_total(adata,
    #                                                         target_sum=1e4,
    #                                                         exclude_highly_expressed=False,
    #                                                         inplace=False)['X']), index=counts.index, columns= counts.columns)
    # import rpy2.robjects as ro
    # from rpy2.robjects import pandas2ri

    # # Activate pandas conversion for rpy2
    # pandas2ri.activate()

    # # Load R libraries
    # ro.r('library(Seurat)')
    # "done"
    # counts_r = pandas2ri.py2rpy(counts_orig.T)

    # ro.globalenv['counts'] = counts_r
    # ro.r('''
    # seurat_obj <- CreateSeuratObject(counts = counts)
    # seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
    # normalized_counts <- as.data.frame(seurat_obj[["SCT"]]@scale.data)
    # ''')
    # normalized_counts = ro.r('normalized_counts')
    # normalized_counts = pandas2ri.rpy2py(normalized_counts)
    # min_value = normalized_counts.min().min()
    # if min_value < 0:
    #     normalized_counts += abs(min_value)  # Shift all values to make them non-negative  
    
      
    # Normalize the data
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log-transform the data
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

    # Scale the data
    sc.pp.scale(adata, max_value=10)

    # Perform PCA
    try:
        sc.tl.pca(adata, svd_solver='arpack',n_comps=n_pcs)
    except:
        print(f'Can not compute {n_pcs}; most likely not enough data')
        exit()

    # Extract the PCA loadings for the specified number of PCs
    loadings = pd.DataFrame(adata.obsm['X_pca'][:, :n_pcs], index=adata.obs.index)
    
    loadings.columns = [f'xPC{i+1}' for i in range(n_pcs)]
    del adata
    # counts_orig = quantile_normalize_rows(counts_orig.T).T
    counts_orig = counts_orig.reset_index()
    loadings = loadings.reset_index(drop=True)
    genotype_pcs = genotype_pcs.reset_index(drop=True)
    # Append the PCs to the original expression data
    
    combined_data = pd.concat([counts_orig, loadings, genotype_pcs], axis=1)
    combined_data = combined_data.set_index(pheno_id)
    # Save the combined data to the specified output file
    combined_data.to_csv(output_file, sep='\t', index=True, chunksize=50000)
    s1 = ",".join(genotype_pcs.columns)
    s2 = ",".join(loadings.columns)
    with open(f"./covariates_new.txt", 'w') as file:
        file.write(f"{s1},{s2}\n")  
        file.write(f"{s1}")     

if __name__ == "__main__":
    # Expecting the number of PCs and the input/output file paths as arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    n_pcs = int(sys.argv[3])
    covs = sys.argv[4]
    run_pca_and_save(input_file, output_file, n_pcs, covs)
