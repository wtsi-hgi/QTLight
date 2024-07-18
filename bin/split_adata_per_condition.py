#!/usr/bin/env python

__author__ = 'Matiss Ozols and Tobi Alegbe'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import scanpy as sc
import argparse
import os
import re

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

    parser.add_argument(
        '-agg', '--agg_columns',
        action='store',
        dest='agg_columns',
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

    # if options.genotype_phenotype:
    #     genotype_phenotype = options.genotype_phenotype
    #     genotype_phenotype = pd.read_csv(genotype_phenotype)
    # else:
        # here we estimate the genotype phenotype interaction file from the genotype, since the IDs are the same. 
    print('Reading in data...')
    adata = sc.read_h5ad(filename=h5ad)


    for agg_col in agg_columns:
        print(agg_col)
        print("----------")
        try:
            data_col = adata.obs[agg_col]
        except:
            print(f'Agregation column {agg_col} doesnt exist in adata')
            continue
            
        for type in data_col.unique():
            print(type)
            print("----------")
            cell_adata = adata[adata.obs[agg_col]==type]
            agg_col_cleaned = re.sub(r'\W+', '_', agg_col.replace(' ', '_'))
            cell_adata.write(
                f'{agg_col_cleaned}__{type}__split.h5ad',
                compression='gzip'
            )     

if __name__ == '__main__':
    main()