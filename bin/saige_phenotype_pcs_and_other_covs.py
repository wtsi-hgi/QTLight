#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2024-09-05'
__version__ = '0.0.1'

import pandas as pd
import numpy as np
import scanpy as sc
from kneed import KneeLocator
import sys

def run_pca_and_save(file_path, output_file, n_pcs, covs):
    # Load the data from the flat text file
    counts = pd.read_csv(file_path, sep="\t", header=0, index_col=0)
    with open(covs, 'r') as file:
        first_line = file.readline().strip()
    existing_covs = first_line.split(',')
    genotype_pcs = counts[existing_covs]
    counts = counts.drop(columns=existing_covs, errors='ignore')
    counts_orig = counts.copy()
    # Convert to AnnData object
    adata = sc.AnnData(X=counts)

    adata.var_names = counts.columns
    adata.obs_names = counts.index

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
    
    counts_orig = counts_orig.reset_index()
    loadings = loadings.reset_index(drop=True)
    genotype_pcs = genotype_pcs.reset_index(drop=True)
    # Append the PCs to the original expression data
    combined_data = pd.concat([counts_orig, loadings, genotype_pcs], axis=1)
    combined_data = combined_data.set_index('pheno_id')
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
