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
        '-nr_gPCs', '--nr_gPCs',
        action='store',
        dest='nr_gPCs',
        required=True,
        default=0,
        help=''
    )

    options = parser.parse_args()
    genotype_pcs=options.genotype_pcs

    
    phenotype_pcs=options.phenotype_pcs
    covariates_df = pd.read_csv(genotype_pcs, sep='\t', index_col=0)

    sample_map_file=options.sample_mapping
    sample_mapping = pd.read_csv(sample_map_file,sep='\t')
    sample_mapping= sample_mapping.set_index('RNA')
        
    if options.phenotype_pcs=='0pcs.tsv':
        phenotype_pcs=pd.DataFrame(index=covariates_df.index)
    else:
        phenotype_pcs= pd.read_csv(phenotype_pcs, sep='\t', index_col=0)
        phenotype_pcs.index = sample_mapping.loc[phenotype_pcs.index]['Genotype']
    
    idx = list(set(phenotype_pcs.index).intersection(set(covariates_df.index)))
    if len(idx)==0:
        covariates_df = covariates_df.T
        idx = list(set(phenotype_pcs.index).intersection(set(covariates_df.index)))  

    covariates_df=covariates_df.rename(columns={'IID':'Genotype'})
    covariates_df=covariates_df.rename(columns={'#IID':'Genotype'})
    
    try:
        covariates_df=covariates_df.set_index('Genotype')
    except:
        print('col already set')
        
    covariates_df = covariates_df.iloc[:,:int(options.nr_gPCs)]
    


    covariates_df = covariates_df.add_prefix('Genotype ')
    phenotype_pcs = phenotype_pcs.add_prefix('Phenotype ')


    covariates_df = covariates_df.loc[idx]
    phenotype_pcs = phenotype_pcs.loc[list(idx)]

    all = {}
    count=0
    for index, row in phenotype_pcs.iterrows():
        d = {**dict(covariates_df.loc[index]), **dict(row)}
        d = {**d, 'id': f"{index}"}
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
        idx2 = set(exctra_covs.columns).intersection(set(data.columns))
        try:
            data = pd.concat([data.loc[:,list(idx2)],exctra_covs.loc[:,list(idx2)]])
        except:
            exctra_covs=exctra_covs.T
            data = pd.concat([data.loc[:,list(idx2)],exctra_covs.loc[:,list(idx2)]])
    else:
        print('no')
    data.to_csv('Covariates.tsv',sep='\t')

if __name__ == '__main__':
    # Convert files to the BED format
    main()