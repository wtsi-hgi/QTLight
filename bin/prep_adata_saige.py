#!/usr/bin/env python

__author__ = 'Bradley Harris'
__date__ = '2023-11-23'
__version__ = '0.0.1'

####### Python script to prepare input files for SAIGE

# Load libraries
##################################################
#### Bradley August 2023
# Atlassing the rectum scRNAseq
# Using all of the expression data, not just limited to that which is genotyped.
##################################################

# Load packages
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib as mp
import kneed as kd
import math
import scipy.stats as st
from scipy import stats
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from scipy.optimize import fsolve
from scipy.optimize import brentq
import argparse
from numpy import asarray
from numpy import savetxt
print("Loaded libraries")

# Define covariate process function
def preprocess_covariates(df, scale_covariates):
    processed_df = df.copy()
    for column in df.columns:
        if pd.api.types.is_categorical_dtype(df[column]):  # Check if the column is categorical
            # Create dummy variables for categorical columns
            dummy_columns = pd.get_dummies(df[column], prefix=column, drop_first=True)
            # Remove the original categorical column
            processed_df.drop(column, axis=1, inplace=True)
            # Replace True with 1, False with 0
            dummy_columns.iloc[:, 0] = dummy_columns.iloc[:, 0].replace({True: 1, False: 0})
            # Replace the dummy variable column name with the original column name
            dummy_columns.columns = [column]
            # Concatenate the dummy columns to the processed DataFrame
            processed_df = pd.concat([processed_df, dummy_columns], axis=1)
        elif df[column].dtype == 'float64' and scale_covariates == "true":  # Check if the column is of type float
                # Scale and center numerical columns
                scaler = StandardScaler()
                processed_df[column] = scaler.fit_transform(df[[column]])
    return processed_df

# Inherit option
def parse_options():    
    # Inherit options
    parser = argparse.ArgumentParser(
            description="""
                Prepping files for the SAIGEQTL analysis
                """
        )

    parser.add_argument(
            '-ph', '--phenotype__file',
            action='store',
            dest='phenotype__file',
            required=True,
            help=''
        )

    parser.add_argument(
            '-a', '--aggregate_on',
            action='store',
            dest='aggregate_on',
            required=True,
            help=''
        )

    parser.add_argument(
            '-g', '--genotype_pc__file',
            action='store',
            dest='genotype_pc__file',
            required=True,
            help=''
        )

    parser.add_argument(
            '-id', '--genotype_id',
            action='store',
            dest='genotype_id',
            required=True,
            help=''
        )

    parser.add_argument(
            '-s', '--sample_id',
            action='store',
            dest='sample_id',
            required=True,
            help=''
        )


    parser.add_argument(
            '-o', '--general_file_dir',
            action='store',
            dest='general_file_dir',
            required=True,
            help=''
        )

    parser.add_argument(
            '-p', '--nperc',
            action='store',
            dest='nperc',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-m', '--min',
            action='store',
            dest='min',
            required=True,
            help=''
        )

    parser.add_argument(
            '-col', '--condition_col',
            action='store',
            dest='condition_col',
            required=False,
            default='',
            help=''
        )

    parser.add_argument(
            '-cond', '--condition',
            action='store',
            dest='condition',
            required=False,
            default='',
            help=''
        )

    parser.add_argument(
            '-covs', '--covariates',
            action='store',
            dest='covariates',
            required=True,
            help=''
        )

    parser.add_argument(
            '-xpca', '--expression_pca',
            action='store',
            dest='expression_pca',
            required=True,
            help=''
        )

    parser.add_argument(
            '-sc', '--scale_covariates',
            action='store',
            dest='scale_covariates',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-br', '--bridge',
            action='store',
            dest='bridge',
            required=False,
            default=None,
            help=''
        )
    
    parser.add_argument(
            '-l', '--level',
            action='store',
            dest='level',
            required=True,
            help=''
        )
    
    return parser.parse_args()



# Define the main script
def main():
    inherited_options = parse_options()
    phenotype__file = inherited_options.phenotype__file
    aggregate_on = inherited_options.aggregate_on
    aggregate_on='Azimuth:predicted.celltype.l1'
    genotype_pc_file = inherited_options.genotype_pc__file
    genotype_id = inherited_options.genotype_id
    sample_id = inherited_options.sample_id
    general_file_dir = inherited_options.general_file_dir
    nperc = float(inherited_options.nperc)*10
    min_cells = int(inherited_options.min)
    condition_col = inherited_options.condition_col
    condition = inherited_options.condition
    covariates = inherited_options.covariates
    expression_pca = inherited_options.expression_pca
    scale_covariates = inherited_options.scale_covariates
    bridge = inherited_options.bridge #'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/work/cc/5a5daedc5c9805a6fd355f479de666/gt_mapping_test.tsv'
    level = inherited_options.level
    print(f"~~~~~~~~~~~~Working on: {level}~~~~~~~~~~~~~~~~")

    # Load in the adata object
    print("Loading object")
    adata = ad.read_h5ad(phenotype__file)

    # Load the genotype PCs
    print("Loading genotype PCs")
    geno_pcs = pd.read_csv(genotype_pc_file, sep = "\t")
    geno_pcs = geno_pcs.set_index("#IID")
    geno_pcs = geno_pcs.iloc[:,1:]
    geno_pcs.reset_index(inplace=True)
    geno_pcs.rename(columns={"#IID": genotype_id}, inplace=True)

    # Subset for the cells we want here (based on the condition column and value specified)
    if condition_col != "NULL":
        print("Subsetting for the condition")
        adata = adata[adata.obs[condition_col] == condition]

    # Replace the ages of those missing [SPECIFIC TO OUR DATA/SAMPLES]
    # adata.obs.loc[adata.obs['sanger_sample_id'].isin(['OTARscRNA9294497', 'OTARscRNA9294498']), 'age_imputed'] = 56.5
    covs_use=covariates.split(',')

    # Define savedir
    if condition_col != "NULL":
        savedir=f"{general_file_dir}/{aggregate_on}/{level}/{condition_col}/{condition}"
    else:
        savedir=f"{general_file_dir}/{aggregate_on}/{level}"
    
    print(f"Files will be saved in {savedir}")
    
    
    if os.path.exists(savedir) == False:
        os.makedirs(savedir, exist_ok=True)
    print("Filtering anndata")
    temp = adata[adata.obs[aggregate_on] == level]
    
    # Filter for intersection with genotypes
    if bridge!=None:
        br1 = pd.read_csv(bridge,sep='\t')
        br1 = br1.set_index('RNA')
        ob1 = pd.DataFrame(temp.obs[genotype_id])
        ob1 = ob1.reset_index().set_index(genotype_id,drop=False)
        ob1[genotype_id]=br1['Genotype']
        ob1 = ob1.set_index('index')
        temp.obs[genotype_id]=ob1[genotype_id]
        
    temp = temp[temp.obs[genotype_id].isin(geno_pcs[genotype_id])]
    
    # prep counts
    print("Filtering lowly expressed genes")
    counts=temp.X
    counts = pd.DataFrame.sparse.from_spmatrix(counts)
    counts.columns = temp.var.index.values
    counts = counts.loc[:, counts.sum() > 0]
    counts[genotype_id] = temp.obs[genotype_id].values.astype(str)
    counts.index = temp.obs.index
    
    # Subset samples if doing so
    if min_cells != "NULL":
        print(f"Working on min cells/sample of {min_cells}")
        cells_per_sample = counts[genotype_id].value_counts()
        keep_samples = cells_per_sample[cells_per_sample > min_cells].index
        counts = counts[counts[genotype_id].isin(keep_samples)]
    
    # Subset the genes based on % expressing samples
    counts_per_sample = counts.groupby(genotype_id).sum()
    min_samples = math.ceil(len(np.unique(temp.obs[genotype_id]))*(int(nperc)/100))
    count_per_column = counts_per_sample.astype(bool).sum(axis=0)
    keep_genes = np.where(count_per_column > min_samples)[0]
    counts = counts.iloc[:,keep_genes]
    print(f"Final shape is:{counts.shape}") 
    
    # Preprocess the covariates (scale continuous, dummy for categorical)
    print("Extracting and sorting covariates")
    to_add = temp.obs[covs_use]
    # Preprocess the covariates (scale continuous, dummy for categorical)
    to_add = preprocess_covariates(to_add, scale_covariates)
    # Bind this onto the counts
    counts = counts.merge(to_add, left_index=True, right_index=True)
    # Add the donor ID (genotyping ID so that we match the genotypes)
    counts[genotype_id] = temp.obs[genotype_id]
    
    # Also add the genotyping PCs
    index_use = counts.index
    counts = counts.merge(geno_pcs, on=genotype_id, how='left')
    counts.set_index(index_use, inplace=True)
    
    # If a covariate has ':' in it's name, this will throw errors in SAIGE, replace this with '_'
    counts.columns=counts.columns.str.replace(':', '_')
    
    # Compute and add the expression PCs
    if expression_pca == "true":
        print("Computing expression PCs")
        sc.pp.normalize_total(temp, target_sum=1e4)
        sc.pp.log1p(temp)
        sc.pp.highly_variable_genes(temp, flavor="seurat", n_top_genes=2000)
        sc.pp.scale(temp, max_value=10)
        sc.tl.pca(temp, svd_solver='arpack')
        pca_variance=pd.DataFrame({'x':list(range(1, 51, 1)), 'y':temp.uns['pca']['variance']})
        # Identify 'knee'
        knee=kd.KneeLocator(x=list(range(1, 51, 1)), y=temp.uns['pca']['variance'], curve="convex", direction = "decreasing")
        knee_point = knee.knee
        print('Knee: ', knee_point)
        # Save knee
        np.savetxt(f"{savedir}/knee.txt", [knee_point], delimiter=',', fmt='%s')
        # Append the PC matrix onto the count data (up to 20)
        loadings = pd.DataFrame(temp.obsm['X_pca'])
        loadings = loadings.iloc[:,0:20]
        loadings.index = temp.obs.index
        loadings.rename(columns=lambda x: f'xPC{x+1}', inplace=True)
        counts = counts.merge(loadings, left_index=True, right_index=True)
        
    # Save this as a dataframe in a directory specific to this resolution - make this if not already
    print("Saving")
    # Save list of gene names to test
    with open(f"{savedir}/test_genes.txt", 'w') as file:
        for item in counts.columns.values:
            file.write(f"{item}\n")
            
            
    # Save counts + covariates
    # NOTE: This currently isn't implemented to save sparse or compressed files (as SAIGE requires dense input) - 80k cells and 10k genes = ~2.5GB file
    #counts.to_csv(f"{savedir}/saige_filt_expr_input.txt", sep = "\t", index=False)
    # Do this per gene with the maximum number of expression PCs - may end up being quicker in SAIGE
    gene_savdir=f"{savedir}/per_gene_input_files"
    if os.path.exists(gene_savdir) == False:
        os.makedirs(gene_savdir, exist_ok=True)
    counts.to_csv(f"{savedir}/saige_filt_expr_input.tsv", sep = "\t", index=True)
    # Save per gene
    # for gene_name in counts.columns.values:
    #     selected_columns = [col for col in counts.columns if col == gene_name or not col.startswith('ENSG')]
    #     new_df = counts[selected_columns]
    #     new_df.to_csv(f"{gene_savdir}/{gene_name}_saige_filt_expr_input.txt", sep = "\t", index=False)



if __name__ == '__main__':
    main()
    
