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
        help=''
    )

    parser.add_argument(
        '-sm', '--sample_mapping',
        action='store',
        dest='sample_mapping',
        required=True,
        help=''
    )

    parser.add_argument(
        '-ppc', '--phenotype_pcs',
        action='store',
        dest='phenotype_pcs',
        required=True,
        help=''
    )

    parser.add_argument(
        '-pfile', '--pfile',
        action='store_true',
        dest='pfile',
        default=False,
        help=''
    )

    options = parser.parse_args()
    genotype_pcs=options.genotype_pcs
    pfile=options.pfile
    covariates_df = pd.read_csv(genotype_pcs, sep='\t', index_col=0)

    if pfile:
        covariates_df.index.names = ['Genotype']
    else:
        covariates_df=covariates_df.rename(columns={'IID':'Genotype'})
        try:
          covariates_df=covariates_df.set_index('Genotype')
        except:
          print('col already set')

    phenotype_pcs=options.phenotype_pcs
    phenotype_pcs= pd.read_csv(phenotype_pcs, sep='\t', index_col=0)

    sample_map_file=options.sample_mapping
    sample_mapping = pd.read_csv(sample_map_file,sep='\t')
    sample_mapping= sample_mapping.set_index('RNA')
    phenotype_pcs.index = sample_mapping.loc[phenotype_pcs.index]['Genotype']
    covariates_df = covariates_df.add_prefix('Genotype ')
    phenotype_pcs = phenotype_pcs.add_prefix('Phenotype ')

    idx = list(set(phenotype_pcs.index).intersection(set(covariates_df.index)))
    covariates_df = covariates_df.loc[idx]
    phenotype_pcs = phenotype_pcs.loc[list(idx)]

    all = {}
    count=0
    for index, row in phenotype_pcs.iterrows():
        d = dict(covariates_df.loc[index])|dict(row)
        d = d|{'id':index}
        all[count]=d
        count+=1
        # print()
    data = pd.DataFrame(all)
    data =data.T
    data=data.set_index('id').sort_index()
    
    data =data.T
    if (options.sample_covariates):
        print('yes')
        exctra_covs = pd.read_csv(options.sample_covariates,sep='\t',index_col=0)
        data = pd.concat([data,exctra_covs])
    else:
        print('no')
    data.to_csv('Covariates.tsv',sep='\t')

if __name__ == '__main__':
    # Convert files to the BED format
    main()