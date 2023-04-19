#!/usr/bin/env python

__author__ = 'Tobi Alegbe'
__date__ = '2023-04-18'
__version__ = '0.0.1'

# https://github.com/broadinstitute/tensorqtl
# https://zenodo.org/record/4118403#.YHclzGMo9TY

import numpy as np
import scipy.stats as stats
import glob
def main():
    import torch
    import pandas as pd
    import tensorqtl
    from tensorqtl import read_phenotype_bed, genotypeio, cis, calculate_qvalues
    print('PyTorch {}'.format(torch.__version__))
    print('Pandas {}'.format(pd.__version__))
    print('Tensorqtl {}'.format(tensorqtl.__version__))
    import argparse
    import os

    os.system('python -V')
    os.system('which python')
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

    # parser.add_argument(
    #     '-gp', '--genotype_phenotype',
    #     action='store',
    #     dest='genotype_phenotype',
    #     required=false,
    #     help=''
    # )     nperm

    parser.add_argument(
        '-cov', '--covariates_file',
        action='store',
        dest='covariates_file',
        required=True,
        help=''
    )

    parser.add_argument(
        '-window', '--window',
        action='store',
        dest='window',
        required=True,
        help=''
    )

    parser.add_argument(
        '-bed', '--expression_bed',
        action='store',
        dest='expression_bed',
        required=True,
        help=''
    )
    parser.add_argument(
        '-plink', '--plink_prefix_path',
        action='store',
        dest='plink_prefix_path',
        required=True,
        help=''
    )
    parser.add_argument(
        '-inter', '--interaction_file',
        action='store',
        dest='inter',
        required=True,
        help=''
    )
    parser.add_argument(
        '-o', '--outdir',
        action='store',
        dest='outdir',
        required=False,
        default='.',
        help=''
    )

    options = parser.parse_args()
    # ValueError: The BED file must define the TSS/cis-window center, with start+1 == end.
    # --plink_prefix_path plink_genotypes/plink_genotypes --expression_bed Expression_Data.bed.gz --covariates_file gtpca_plink.eigenvec
    plink_prefix_path=options.plink_prefix_path
    expression_bed=options.expression_bed
    covariates_file=options.covariates_file
    outdir=options.outdir


    phenotype_df, phenotype_pos_df = read_phenotype_bed(expression_bed)
    
    
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
    phenotype_df = phenotype_df[covariates_df.columns]
    # have to drop dublicate rownames. and average the repeated measures.
    phenotype_df.columns = phenotype_df.columns.str.split('.').str[0]
    covariates_df.columns = covariates_df.columns.str.split('.').str[0]

    covariates_df=covariates_df.loc[:,~covariates_df.columns.duplicated()]
    # this can be adjusted to take an average. TQTL can not account for repeated measures.
    phenotype_df=phenotype_df.loc[:,~phenotype_df.columns.duplicated()]

    covariates_df=covariates_df.T
    # not a good solution but atm

    # covariates_df=covariates_df.set_index('IID')
    # to_keep = list(set(covariates_df.index).intersection(set(phenotype_df.columns)))
    # covariates_df=covariates_df.loc[to_keep]
    # covariates_df= covariates_df

    # phenotype_df = phenotype_df[to_keep]

    print('----Fine read ------')
    if torch.cuda.is_available():
        print(f'  * using GPU ({torch.cuda.get_device_name(torch.cuda.current_device())})')
    else:
        print('  * WARNING: using CPU!')
    pr = genotypeio.PlinkReader(plink_prefix_path)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    Directory = f'{outdir}/inter_output'
    os.mkdir(Directory)
    cis.map_nominal(genotype_df, variant_df, 
                    phenotype_df.loc[phenotype_pos_df['chr']!='chrY'], 
                    phenotype_pos_df.loc[phenotype_pos_df['chr']!='chrY'],
                    covariates_df=covariates_df,prefix='cis_inter1',
                    output_dir=Directory, write_top=True, write_stats=True,
                    interaction_s=interaction_s, maf_threshold_interaction=0.05,
                    run_eigenmt=True)
    
    all_files = glob.glob(f'{Directory}/cis_inter*.parquet')
    All_Data = pd.DataFrame()
    count=0
    for bf1 in all_files:
        print(bf1)
        df = pd.read_parquet(bf1)
        df.to_csv(bf1.replace('.parquet','.tsv'),sep='\t',index=False)
        os.remove(bf1) 
        count+=1    

if __name__ == '__main__':
    main()
