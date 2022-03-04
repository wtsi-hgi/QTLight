#!/usr/bin/env python

__author__ = 'Matiss Ozols and Hannes Ponstingl'
__date__ = '2021-11-25'
__version__ = '0.0.1'

# https://github.com/broadinstitute/tensorqtl
# https://zenodo.org/record/4118403#.YHclzGMo9TY

import numpy as np
import scipy.stats as stats
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
        '-nperm', '--nperm',
        action='store',
        dest='nperm',
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

    options = parser.parse_args()
    # ValueError: The BED file must define the TSS/cis-window center, with start+1 == end.
    # --plink_prefix_path plink_genotypes/plink_genotypes --expression_bed Expression_Data.bed.gz --covariates_file gtpca_plink.eigenvec
    plink_prefix_path = 'plink_genotypes/plink_genotypes'
    expression_bed = 'Expression_Data.bed.gz'
    covariates_file = 'gtpca_plink.eigenvec'
    # /lustre/scratch123/hgi/teams/hgi/mo11/eQTL_mapping/LIMIX/work/62/c8dd517ac8a526214952fb551f1c25
    # prefix = 'bin/HipSci_all'
    # extra = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/eqtl/franke_data/work/ef/7573334e4e1958378cc11ac8d33f59/normaggrsum_NK_counts_chrAll.bed.gz'
    plink_prefix_path=options.plink_prefix_path
    expression_bed=options.expression_bed
    covariates_file=options.covariates_file


    phenotype_df, phenotype_pos_df = read_phenotype_bed(expression_bed)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
    covariates_df=covariates_df.set_index('IID')
    to_keep = list(set(covariates_df.index).intersection(set(phenotype_df.columns)))
    covariates_df=covariates_df.loc[to_keep]
    covariates_df= covariates_df

    phenotype_df = phenotype_df[to_keep]

    print('----Fine read ------')
    pr = genotypeio.PlinkReader(plink_prefix_path)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

    # cis.map_nominal(genotype_df, variant_df,
    #                 phenotype_df.loc[phenotype_pos_df['chr']!='chrY'],
    #                 phenotype_pos_df.loc[phenotype_pos_df['chr']!='chrY'],
    #                 prefix, covariates_df=covariates_df,maf_threshold_interaction=0.05,
    #                 run_eigenmt=True, output_dir='.', write_top=True, write_stats=True)

    cis_df = cis.map_cis(genotype_df, variant_df, 
                        phenotype_df.loc[phenotype_pos_df['chr']!='chrY'],
                        phenotype_pos_df.loc[phenotype_pos_df['chr']!='chrY'],nperm=int(options.nperm),
                        window=int(options.window),
                        covariates_df=covariates_df, seed=123456)
    print('----cis eQTLs processed ------')
    cis_df.head()
    cis_df.to_csv("Cis_eqtls.tsv",sep="\t")
    sv = ~np.isnan(cis_df['pval_beta'])
    print(f"Dropping {sum(sv)} variants withouth Beta-approximated p-values to\n.")
    cis_df_dropped = cis_df.loc[sv]
    r = stats.pearsonr(cis_df_dropped['pval_perm'], cis_df_dropped['pval_beta'])[0]
    calculate_qvalues(cis_df_dropped, qvalue_lambda=0.85)
    cis_df_dropped.to_csv("Cis_eqtls_qval.tsv", sep='\t')

if __name__ == '__main__':
    main()

# trans_df = tensorqtl.trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=10000,
#                            return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
# trans_df = tensorqtl.trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=5000000)
# trans_df.to_csv("Trans_eqtls.tsv",sep="\t")