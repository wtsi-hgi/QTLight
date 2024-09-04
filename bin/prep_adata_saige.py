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
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import kneed as kd
from sklearn.preprocessing import StandardScaler
import argparse
import gc
import sys
import gc

# Define covariate process function
def preprocess_covariates(df, scale_covariates):
    processed_df = pd.get_dummies(df, drop_first=True)  # One-hot encode all categorical columns
    if scale_covariates == "true":
        scaler = StandardScaler()
        for column in processed_df.select_dtypes(include=['float64']).columns:
            processed_df[column] = scaler.fit_transform(processed_df[[column]])
    return processed_df

def parse_options():    
    parser = argparse.ArgumentParser(description="Prepping files for the SAIGEQTL analysis")
    parser.add_argument('-ph', '--phenotype__file', required=True)
    parser.add_argument('-a', '--aggregate_on', required=True)
    parser.add_argument('-g', '--genotype_pc__file', required=True)
    parser.add_argument('-id', '--genotype_id', required=True)
    parser.add_argument('-s', '--sample_id', required=True)
    parser.add_argument('-o', '--general_file_dir', required=True)
    parser.add_argument('-p', '--nperc', required=True)
    parser.add_argument('-m', '--min', required=True)
    parser.add_argument('-col', '--condition_col', default='')
    parser.add_argument('-cond', '--condition', default='')
    parser.add_argument('-covs', '--covariates', default='')
    parser.add_argument('-xpca', '--expression_pca', required=True)
    parser.add_argument('-chr', '--chr', required=False, default=None)
    parser.add_argument('-genome', '--genome', required=False, default=None)
    parser.add_argument('-sc', '--scale_covariates', required=True)
    parser.add_argument('-br', '--bridge', default=None)
    return parser.parse_args()

def main():
    inherited_options = parse_options()
    phenotype__file = inherited_options.phenotype__file
    aggregate_on = inherited_options.aggregate_on
    genotype_pc_file = inherited_options.genotype_pc__file
    genotype_id = inherited_options.genotype_id
    general_file_dir = inherited_options.general_file_dir
    nperc = float(inherited_options.nperc) * 100
    min_cells = int(inherited_options.min)
    condition_col = inherited_options.condition_col
    condition = inherited_options.condition
    covariates = inherited_options.covariates.split(',') if inherited_options.covariates else []
    expression_pca = int(inherited_options.expression_pca)
    scale_covariates = inherited_options.scale_covariates
    bridge = inherited_options.bridge

    print("Loading object")
    adata = ad.read_h5ad(phenotype__file,backed='r')
    genes=list(adata.var.index)
    if (inherited_options.chr):
        # Here we subset down to the genes available on determined chr.
        from gtfparse import read_gtf
        df = read_gtf(inherited_options.genome)
        df = df[df.feature == 'gene']
        Gene_Chr_Start_End_Data =df[['gene_id','start','end','strand','seqname']]
        Gene_Chr_Start_End_Data.rename(columns={'gene_id':'feature_id','seqname':'chromosome'},inplace=True)
        chrs = inherited_options.chr.split(',')
        all_genes = set(Gene_Chr_Start_End_Data[Gene_Chr_Start_End_Data['chromosome'].isin(chrs)]['feature_id'])
        
        genes = list(all_genes.intersection(genes))
        # adata = adata[:, all_genes]
        del df
        del all_genes
        del Gene_Chr_Start_End_Data
        gc.collect()
        
    # condition='Platelet'
    # if condition_col != "NULL":
    #     print("Subsetting for the condition")
    #     conditions=list(set(condition.split(',')))
    #     adata = adata[adata.obs[condition_col].isin(conditions),genes].copy(filename='tmp.h5ad')
    #     gc.collect()

    if len(adata.obs)==0:
        print('The final subset adata is empty, hence we do not perform analysis')
        sys.exit()
    adata.obs.index += adata.obs.index.duplicated().cumsum().astype(str)  # Resolve duplicated indices

    print("Loading genotype PCs")
    geno_pcs = pd.read_csv(genotype_pc_file, sep = "\t")
    try:
        geno_pcs = geno_pcs.set_index("#IID")
    except:
        geno_pcs = geno_pcs.set_index("IID")
    # geno_pcs = geno_pcs.iloc[:,1:]
    geno_pcs.reset_index(inplace=True)
    geno_pcs.rename(columns={"#IID": genotype_id,"IID": genotype_id}, inplace=True)
    geno_pcs.rename(columns={"#FID": genotype_id,"IID": genotype_id}, inplace=True)
    geno_pcs = geno_pcs.set_index(genotype_id)

    levels = set(adata.obs[aggregate_on].unique())
    l1 = len(levels)
    conditions=set(condition.split(','))
    levels=levels.intersection(conditions)
    
    for level in levels:
        print(f"~~~~~~~~~~~~Working on: {level}~~~~~~~~~~~~~~~~")
        savedir = f"{general_file_dir}/{aggregate_on}/{level}"
        
        os.makedirs(savedir, exist_ok=True)
        print("Filtering anndata")

        temp = adata[adata.obs[aggregate_on] == level,genes].to_memory()  # Use copy to avoid modifying original data
        if (l1==1):
            del adata
            gc.collect() 
            
        if bridge:
            br1 = pd.read_csv(bridge, sep='\t').set_index('RNA')
            br2 = temp.obs[genotype_id].map(br1['Genotype'])
            temp.obs[genotype_id] = br2
            del br1
        temp = temp[temp.obs[genotype_id].isin(geno_pcs.index)]        

            
        print("Filtering lowly expressed genes")
        counts = pd.DataFrame.sparse.from_spmatrix(temp.X, index=temp.obs.index, columns=temp.var.index)
        counts = counts.loc[:, counts.sum(axis=0) > 0]
        counts[genotype_id] = temp.obs[genotype_id].values.astype(str)
        
        if min_cells:
            print(f"Working on min cells/sample of {min_cells}")
            cells_per_sample = counts[genotype_id].value_counts()
            keep_samples = cells_per_sample[cells_per_sample > min_cells].index
            counts = counts[counts[genotype_id].isin(keep_samples)]
            del cells_per_sample
        
        temp = temp[temp.obs[genotype_id].isin(keep_samples)]
        del keep_samples
        if nperc > 0:
            # counts_per_sample = counts.groupby(genotype_id).sum()
            # min_samples = int(np.ceil(len(np.unique(genotype_ids)) * (nperc / 100)))
            # count_per_column = np.count_nonzero(counts_per_sample.values, axis=0)
            # keep_genes = np.where(count_per_column > min_samples)[0]
            # counts2 = counts.iloc[:, keep_genes]
            
            # Precompute unique genotypes and minimum samples needed
            genotype_ids = counts[genotype_id].values  # Extract genotype IDs
            unique_genotypes = np.unique(genotype_ids)  # Get unique genotypes
            min_samples = int(np.ceil(len(unique_genotypes) * (nperc / 100)))
            counts_data = counts.drop(columns=[genotype_id]).values  # Convert counts (excluding genotype_id) to NumPy array
            summed_counts = np.zeros((len(unique_genotypes), counts_data.shape[1]))
            # Accumulate sums directly on the NumPy array
            for i, genotype in enumerate(unique_genotypes):
                mask = (genotype_ids == genotype)
                summed_counts[i, :] = counts_data[mask, :].sum(axis=0)
            count_per_column = np.count_nonzero(summed_counts, axis=0)
            keep_genes = np.where(count_per_column > min_samples)[0]
            # Use NumPy slicing for efficient selection
            filtered_counts = counts.iloc[:, keep_genes].copy()
            # Reattach genotype_id
            filtered_counts[genotype_id] = genotype_ids
            counts = filtered_counts
            del counts_data, summed_counts, genotype_ids, unique_genotypes, count_per_column  # Free memory
            
        
        print(f"Final shape is: {counts.shape}") 
        if (any(np.array(counts.shape) < 2)):
            print(f"Final shape is not acceptable, skipping this phenotype: {counts.shape}")
            continue
        
        covariates_string=''
        if covariates:
            print("Extracting and sorting covariates")
            to_add = preprocess_covariates(temp.obs[covariates], scale_covariates)
            counts = counts.join(to_add)
            covariates_string = ','+inherited_options.covariates
        

        
        counts.index = counts.index + counts.index.duplicated().cumsum().astype(str)
        temp.obs.index = temp.obs.index + temp.obs.index.duplicated().cumsum().astype(str)
        
        counts[genotype_id] = temp.obs[genotype_id]
        counts = counts.merge(geno_pcs.reset_index(), on=genotype_id, how='left').set_index(counts.index)
        counts.columns = counts.columns.str.replace(':', '_')
        
        covariates_string = ','.join(geno_pcs.columns.values)+covariates_string
        sample_covariates = ','.join(geno_pcs.columns.values)
        
        with open(f"{savedir}/test_genes.txt", 'w') as file:
            file.write("\n".join(counts.columns))
            
        # if expression_pca > 0:
        #     print("Computing expression PCs")
        #     sc.pp.normalize_total(temp, target_sum=1e4)
        #     sc.pp.log1p(temp)
        #     sc.pp.highly_variable_genes(temp, flavor="seurat", n_top_genes=2000)
        #     sc.pp.scale(temp, max_value=10)
        #     sc.tl.pca(temp, svd_solver='arpack')
        #     pca_variance = pd.DataFrame({'x': range(1, 51), 'y': temp.uns['pca']['variance']})
        #     knee = kd.KneeLocator(x=pca_variance['x'], y=pca_variance['y'], curve="convex", direction="decreasing")
        #     knee_point = knee.knee
        #     np.savetxt(f"{savedir}/knee.txt", [knee_point], delimiter=',', fmt='%s')
        #     loadings = pd.DataFrame(temp.obsm['X_pca'][:, :expression_pca], index=temp.obs.index)
        #     loadings.columns = [f'xPC{i+1}' for i in range(expression_pca)]
        #     counts = counts.join(loadings)
        #     covariates_string = covariates_string+','+','.join(loadings.columns.values)
        
        print("Saving")
        with open(f"{savedir}/covariates.txt", 'w') as file:
                file.write(f"{covariates_string}\n")  
                file.write(f"{sample_covariates}")  
        # gene_savdir = f"{savedir}/per_gene_input_files"
        # os.makedirs(gene_savdir, exist_ok=True)
        counts.set_index(genotype_id).to_csv(f"{savedir}/saige_filt_expr_input.tsv", sep="\t", index=True, chunksize=50000)

        del counts
        del temp
        gc.collect()  # Clean up memory
        try:
            os.remove('tmp.h5ad')
        except:
            _=''
        try:
            os.remove('tmp2.h5ad')
        except:
            _=''
        try:
            os.remove('tmp3.h5ad')
        except:
            _=''


if __name__ == '__main__':
    main()
