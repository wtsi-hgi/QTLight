#!/usr/bin/env python
__author__ = 'Matiss Ozols'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import argparse
import math

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
        '-gpc', '--genotype_pcs',
        action='store',
        dest='genotype_pcs',
        required=True,
        help=''
    )
    parser.add_argument(
        '-sc', '--sample_covariates',
        action='store',
        dest='sample_covariates',
        required=False,
        default='',
        help=''
    )

    parser.add_argument(
        '-nr_gPCs', '--nr_gPCs',
        action='store',
        dest='nr_gPCs',
        required=True,
        default=0,
        help=''
    )

    options = parser.parse_args()
    genotype_pcs=options.genotype_pcs
    covariates_df = pd.read_csv(genotype_pcs, sep='\t', index_col=0)
    covariates_df.index = covariates_df.index.astype(str)
    covariates_df=covariates_df.rename(columns={'IID':'#IID'})
    
    try:
        covariates_df=covariates_df.set_index('#IID')
        covariates_df.index = covariates_df.index.astype(str)
    except:
        print('col already set')
        
    covariates_df = covariates_df.iloc[:,:int(options.nr_gPCs)]
    covariates_df=covariates_df.T

    if (options.sample_covariates):
        print('yes')
        exctra_covs = pd.read_csv(options.sample_covariates,sep='\t',index_col=0)
        idx2 = set(exctra_covs.columns).intersection(set(covariates_df.columns))
        try:
            data = pd.concat([covariates_df.loc[:,list(idx2)],exctra_covs.loc[:,list(idx2)]])
        except:
            exctra_covs=exctra_covs.T
            data = pd.concat([covariates_df.loc[:,list(idx2)],exctra_covs.loc[:,list(idx2)]])
    else:
        print('no')
        data = covariates_df  # fallback to genotype PCs only
    data2 = data.T
    res = data2.reset_index()
    res = res.rename(columns={'index': '#IID'})
    res.to_csv('G_E_Covariates.tsv',sep='\t',index=False)

if __name__ == '__main__':
    # Convert files to the BED format
    main()