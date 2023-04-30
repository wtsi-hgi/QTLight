#!/usr/bin/env python

__author__ = 'Matiss Ozols and Tobi Alegbe'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import scanpy as sc
import argparse
import os


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
    agg_columns = options.agg_columns
    agg_columns = agg_columns.split(",")
    gt_id_column =  options.gt_id_column
    sample_column = options.sample_column
    n_individ = int(options.n_individ)
    n_cells = int(options.n_cells)
    # if options.genotype_phenotype:
    #     genotype_phenotype = options.genotype_phenotype
    #     genotype_phenotype = pd.read_csv(genotype_phenotype)
    # else:
        # here we estimate the genotype phenotype interaction file from the genotype, since the IDs are the same. 
    print('Reading in data...')
    adata = sc.read_h5ad(filename=h5ad)
    adata.obs['adata_phenotype_id'] = adata.obs[gt_id_column].astype('str')+'_'+adata.obs[sample_column].astype('str')
    genotype_phenotype_mapping = []
    aggregated_data=pd.DataFrame()

    for method in methods:
        if (method =='dMean'):
            adata.X = adata.layers['dMean_normalised']
        for agg_col in agg_columns:
            print(agg_col)
            print("----------")
            for type in adata.obs[agg_col].unique():
                print(type)
                print("----------")# 
                # type='CD4 CTL'
                cell_adata = adata[adata.obs[agg_col]==type]
                if (len(cell_adata.obs['adata_phenotype_id'].unique())>n_individ):
                    aggregated_data_pre=pd.DataFrame()
                    genotype_phenotype_mapping_pre = []
                    for individual_1 in cell_adata.obs['adata_phenotype_id'].unique():
                        individual_1_adata = cell_adata[cell_adata.obs['adata_phenotype_id']==individual_1]
                        if(individual_1_adata.obs.shape[0]>n_cells):
                            print(individual_1)
                            Genotype = individual_1_adata.obs[gt_id_column].unique()[0]
                            # Change this to any aggregation strategy
                            #as per https://www.medrxiv.org/content/10.1101/2021.10.09.21264604v1.full.pdf 
                            # We mapped cis-eQTL within a 1 megabase (MB) window of the TSS of each gene expressed
                            # in at least 5% of the nuclei (belonging to a broad cell type)
                            if (method =='dSum'):
                                f = individual_1_adata.to_df()
                                data_aggregated_for_cell_and_individal = pd.DataFrame(f.sum(axis = 0))
                                data_aggregated_for_cell_and_individal.set_index(f.columns,inplace=True)
                                type2= f"{agg_col}-{type}-{method}"
                            elif (method =='dMean'):
                                f = individual_1_adata.to_df()
                                data_aggregated_for_cell_and_individal = pd.DataFrame(f.mean(axis = 0))
                                data_aggregated_for_cell_and_individal.set_index(f.columns,inplace=True)
                                type2= f"{agg_col}-{type}-{method}"
                            else:
                                print('Wrong method specified, please use dMean or dSum or both as a coma seperated sting dMean,dSum')
                                break
                            Phenotype = f"{type2}_{individual_1}".replace(' ','_')
                            type2=type2.replace(' ','_')
                            data_aggregated_for_cell_and_individal.rename(columns={0:Phenotype},inplace=True)
                            aggregated_data_pre=pd.concat([aggregated_data_pre,data_aggregated_for_cell_and_individal],axis=1)
                            genotype_phenotype_mapping_pre.append({'Genotype':Genotype,'RNA':Phenotype,'Sample_Category':type2})
                    # assess whether correct number of individuals ended up having right ammount of cells
                    if (len(aggregated_data_pre.columns)>=n_individ):
                        aggregated_data=pd.concat([aggregated_data,aggregated_data_pre],axis=1)
                        genotype_phenotype_mapping= genotype_phenotype_mapping+ genotype_phenotype_mapping_pre
                        # f = pd.DataFrame(individual_1_adata.X.mean(axis=0))
    genotype_phenotype_mapping = pd.DataFrame(genotype_phenotype_mapping)
    genotype_phenotype_mapping.to_csv('genotype_phenotype_mapping.tsv',sep='\t',index=False)
    aggregated_data.to_csv('phenotype_file.tsv',sep='\t',index=True)
    print('Successfully Finished')




if __name__ == '__main__':
    main()