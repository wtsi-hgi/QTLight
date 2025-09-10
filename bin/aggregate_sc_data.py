#!/usr/bin/env python

__author__ = 'Matiss Ozols and Tobi Alegbe'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import scanpy as sc
import argparse
import os
import re
import gc
import numpy as np
import multiprocessing as mp
from functools import partial

_ADATA = None  # per-process AnnData
_METHOD = None

def _init_pool(h5ad_path, gt_id_column, sample_column, method):
    """Runs once per worker process; opens backed AnnData locally (no pickling)."""
    global _ADATA, _METHOD
    _METHOD = method
    _ADATA = sc.read_h5ad(h5ad_path, backed='r')

    # Make obs index unique and add phenotype id in worker too (to match main)
    _ADATA.obs.index = _ADATA.obs.index + "___" + _ADATA.obs.groupby(_ADATA.obs.index).cumcount().astype(str)
    _ADATA.strings_to_categoricals()
    _ADATA.obs['adata_phenotype_id'] = _ADATA.obs[gt_id_column].astype('str') + '_' + _ADATA.obs[sample_column].astype('str')

    # If using dMean, align X to the layer in worker memory
    if _METHOD == 'dMean' and 'dMean_normalised' in _ADATA.layers:
        # materialize this layer into X within the worker
        _ADATA = _ADATA.to_memory()
        _ADATA.X = _ADATA.layers['dMean_normalised']

def aggregate_individual(
    individual_1, cell_index, indexes, n_cells, gt_id_column, method,
    agg_col="agg", type_name="RNA", sample_covariate_columns=None, log=None
):
    global _ADATA
    adata = _ADATA  # local alias

    # 1) pick cells
    donor_mask = adata.obs['adata_phenotype_id'] == individual_1
    if donor_mask.sum() == 0:
        if log is not None: log.append(f"{individual_1}: no cells in adata")
        return None

    donot_index = set(adata.obs.index[donor_mask])
    cell_donor_index = set(cell_index).intersection(donot_index)
    if not cell_donor_index:
        if log is not None: log.append(f"{individual_1}: no overlap with cell_index")
        return None

    individual_1_adata = adata[list(cell_donor_index), indexes]
    n_obs = individual_1_adata.n_obs
    if n_obs < n_cells:
        if log is not None: log.append(f"{individual_1}: have {n_obs} cells < required {n_cells}")
        return None

    # 2) genotype
    gvals = individual_1_adata.obs[gt_id_column].unique()
    if len(gvals) == 0:
        if log is not None: log.append(f"{individual_1}: no genotype values in '{gt_id_column}'")
        return None
    Genotype = gvals[0]

    # 3) aggregate
    f = individual_1_adata.to_df()
    if method == 'dSum':
        agg = f.sum(axis=0)
    elif method == 'dMean':
        agg = f.mean(axis=0)
    else:
        if log is not None: log.append("Wrong method (use dMean or dSum)")
        return None

    type2 = f"{agg_col}-{type_name}-{method}".replace(' ', '_')
    Phenotype = f"{type2}_{individual_1}".replace(' ', '_')

    df = pd.DataFrame(agg)
    df.set_index(f.columns, inplace=True)
    df.rename(columns={0: Phenotype}, inplace=True)

    mapping = {'Genotype': Genotype, 'RNA': Phenotype, 'Sample_Category': type2}

    cov_dict = None
    if sample_covariate_columns:
        cov_dict = {'RNA': Phenotype}
        for col in sample_covariate_columns:
            vals = individual_1_adata.obs[col].unique()
            if len(vals) == 1:
                cov_dict[col] = vals[0]
            else:
                if log is not None: log.append(f"{individual_1}: cov '{col}' has multiple values, skipped")

    return (df, mapping, cov_dict)

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge 10x data. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )
# --method ${params.aggregation_method}
    # parser.add_argument(
    #     '-gp', '--genotype_phenotype',
    #     action='store',
    #     dest='genotype_phenotype',
    #     required=false,
    #     help=''
    # )    
    parser.add_argument(
        '-w', '--n_workers',
        type=int,
        default=max(1, mp.cpu_count()),
        help='Number of worker processes for parallel aggregation.'
    )
    
    parser.add_argument(
        '-method', '--method',
        action='store',
        dest='method',
        required=True,
        help=''
    )

    parser.add_argument(
        '-agg', '--agg_columns',
        action='store',
        dest='agg_columns',
        required=True,
        help=''
    )
    parser.add_argument(
        '-gt', '--gt_id_column',
        action='store',
        dest='gt_id_column',
        required=True,
        help=''
    )    
    parser.add_argument(
        '-sample', '--sample_column',
        action='store',
        dest='sample_column',
        required=True,
        help=''
    )
    parser.add_argument(
        '-ncells', '--n_cells',
        action='store',
        dest='n_cells',
        required=True,
        help=''
    )
    parser.add_argument(
        '-n_individ', '--n_individ',
        action='store',
        dest='n_individ',
        required=True,
        help=''
    )
    parser.add_argument(
        '-h5ad', '--h5ad',
        action='store',
        dest='h5ad',
        required=True,
        help=''
    )

    parser.add_argument(
        '-cpt', '--cell_percentage_threshold',
        action='store',
        dest='cell_percentage_threshold',
        required=True,
        help='cell_percentage_threshold'
    )
    
    parser.add_argument(
        '-scc', '--covariates',
        action='store',
        dest='sample_covariate_column',
        required=False,
        default='',
        help='Optional obs column to extract as covariate per sample'
    )

    options = parser.parse_args()
    methods = options.method
    methods = methods.split(",")
    # h5ad = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Franke_with_genotypes_nfCore/results/celltype/adata.h5ad'
    # agg_columns = 'Azimuth:predicted.celltype.l2'
    # agg_columns='Azimuth:predicted.celltype.l2,Celltypist:Immune_All_High,Celltypist:Immune_All_Low'
    # gt_id_column = 'donor_id'
    # sample_column = 'convoluted_samplename'
    # n_individ=30
    # n_cells=10
    h5ad = options.h5ad
    if '--' in h5ad:
        prefix=h5ad.split('--')[-1]+'--'
        prefix2=h5ad.split('--')[-1]+'___'
    else:
        prefix=''
        prefix2=''
        
    agg_columns = options.agg_columns
    agg_columns = agg_columns.split(",")
    if len(agg_columns)>1:
        for f1 in agg_columns:
            if re.sub(r'\W+', '_', f1.replace(' ', '_')) in h5ad.split("__")[0]:
                agg_columns = [f1]
             
    gt_id_column =  options.gt_id_column
    sample_column = options.sample_column
    n_individ = int(options.n_individ)
    cell_percentage_threshold = float(options.cell_percentage_threshold)
    n_cells = int(options.n_cells)
    print('Reading in data...')
    adata = sc.read_h5ad(filename=h5ad, backed='r')
    adata.obs.index = adata.obs.index + "___" + adata.obs.groupby(adata.obs.index).cumcount().astype(str)
    # make sure indexes are unique
    adata.strings_to_categoricals()
    adata.obs['adata_phenotype_id'] = adata.obs[gt_id_column].astype('str')+'_'+adata.obs[sample_column].astype('str')
    


    for method in methods:

        if (method =='dMean'):
            if 'dMean_normalised' in adata.layers:
                adata = adata.to_memory()
                adata.X = adata.layers['dMean_normalised']
                print("adata.X has been updated with 'dMean_normalised'.")
            else:
                print('probably already normalised since no layer added')
                
        for agg_col in agg_columns:

            print(agg_col)
            
            print("----------")
            try:
                data_col = adata.obs[agg_col]
            except:
                print(f'Agregation column {agg_col} doesnt exist in adata')
                continue
            for type in data_col.unique():
                genotype_phenotype_mapping = []
                aggregated_data=pd.DataFrame()
            
                print(type)
                print("----------")# 
                modified_agg_col = re.sub(r'[^a-zA-Z0-9]', '_', type)
                adata.strings_to_categoricals()
                # type='CD4 CTL'
                cell_adata = adata[adata.obs[agg_col]==type]
                cell_index = set(adata[adata.obs[agg_col]==type].obs.index)
                
                # cell_percentage_threshold = 0.1  # e.g., 10% of cells must express the gene
                if cell_percentage_threshold > 0:
                    # Calculate the proportion of cells expressing each gene
                    cell_counts = (cell_adata.X > 0).sum(axis=0).A1  # .A1 converts sparse matrix to a flat array
                    total_cells = cell_adata.shape[0]
                    cell_expression_proportion = cell_counts / total_cells
                    # Apply the filter based on cell-level expression
                    keep_genes = cell_expression_proportion >= cell_percentage_threshold
                    # Get the index IDs of the retained genes
                    indexes = cell_adata.var.index[keep_genes].tolist()
                else:
                    indexes = list(cell_adata.var.index)
                if (len(cell_adata.obs['adata_phenotype_id'].unique()) >= n_individ):
                    aggregated_data_pre = pd.DataFrame()
                    sample_covariates_pre = []
                    genotype_phenotype_mapping_pre = []

                    mgr = mp.Manager()
                    shared_log = mgr.list()

                    individual_ids = list(cell_adata.obs['adata_phenotype_id'].unique())
                    sample_covariate_columns = options.sample_covariate_column.split(',') if options.sample_covariate_column else None

                    # Build starmap args (NO adata)
                    map_args = [
                        (individual_1, cell_index, indexes, n_cells, gt_id_column, method, agg_col, type, sample_covariate_columns, shared_log)
                        for individual_1 in individual_ids
                    ]

                    n_workers = max(1, int(getattr(options, 'n_workers', mp.cpu_count())))
                    ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else mp.get_context()

                    # Important: create the pool ONCE per (method, agg_col, type) batch
                    with ctx.Pool(
                        processes=n_workers,
                        initializer=_init_pool,
                        initargs=(h5ad, gt_id_column, sample_column, method)
                    ) as pool:
                        results = pool.starmap(
                            aggregate_individual,
                            map_args,
                            chunksize=max(1, len(map_args) // (n_workers * 4) or 1)
                        )

                    for result in results:
                        if result is None:
                            continue
                        df_part, mapping_dict, cov_dict = result
                        aggregated_data_pre = pd.concat([aggregated_data_pre, df_part], axis=1)
                        genotype_phenotype_mapping_pre.append(mapping_dict)
                        if cov_dict is not None:
                            sample_covariates_pre.append(cov_dict)
                    
                    total_individuals = len(individual_ids)
                    kept_individuals  = len(aggregated_data_pre.columns)
                    status = "PASSED" if kept_individuals >= n_individ else "FAILED"

                    
                    with open(f"{prefix2}{method}__{modified_agg_col}___log.txt", "w") as fh:
                        fh.write(f"[SUMMARY] method={method} agg_col={agg_col} celltype={type}\n")
                        fh.write(f"[COUNTS] total={total_individuals} kept={kept_individuals}\n")
                        fh.write(f"[THRESHOLD] required_n_individ={n_individ} -> {status}\n")
                        if len(shared_log) > 0:
                            for line in shared_log:
                                fh.write(line + "\n")
                        # --- end parallel block ---

                    # assess whether correct number of individuals ended up having right ammount of cells
                    if (len(aggregated_data_pre.columns)>=n_individ):
                        aggregated_data=pd.concat([aggregated_data,aggregated_data_pre],axis=1)
                        genotype_phenotype_mapping= genotype_phenotype_mapping+ genotype_phenotype_mapping_pre
                else:
                    with open(f"{prefix2}{method}__{modified_agg_col}___log.txt", "w") as fh:
                        fh.write(f"[SUMMARY] method={method} agg_col={agg_col} celltype={type}\n")
                        fh.write(f"[COUNTS] total={len(cell_adata.obs['adata_phenotype_id'].unique())}\n")
                        fh.write(f"[THRESHOLD] required_n_individ={n_individ} -> FAILED\n")
                    
                genotype_phenotype_mapping = pd.DataFrame(genotype_phenotype_mapping)
                if(len(genotype_phenotype_mapping)>=10):
                    genotype_phenotype_mapping.to_csv(f'{prefix}{method}__{modified_agg_col}___genotype_phenotype_mapping.tsv',sep='\t',index=False)
                    aggregated_data.to_csv(f'{prefix}{method}__{modified_agg_col}___phenotype_file.tsv',sep='\t',index=True)
                if options.sample_covariate_column and len(sample_covariates_pre) > 0:
                    pd.DataFrame(sample_covariates_pre).to_csv(
                        f'{prefix}{method}__{modified_agg_col}___sample_covariates.tsv', sep='\t', index=False
                    )
    print('Successfully Finished')

if __name__ == '__main__':
    main()
