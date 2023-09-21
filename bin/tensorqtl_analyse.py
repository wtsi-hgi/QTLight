#!/usr/bin/env python

__author__ = 'Matiss Ozols and Hannes Ponstingl'
__date__ = '2021-11-25'
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
    from tensorqtl import read_phenotype_bed, genotypeio, cis, calculate_qvalues,pgen 
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

    parser.add_argument(
        '-o', '--outdir',
        action='store',
        dest='outdir',
        required=False,
        default='.',
        help=''
    )

    parser.add_argument(
        '-dosage', '--dosage',
        action='store_true',
        dest='dosage',
        default=False,
        help=''
    )

    parser.add_argument(
        '-maf', '--maf',
        action='store',
        dest='maf',
        required=False,
        default=0.05,
        help=''
    )

    options = parser.parse_args()
    maf=float(options.maf)
    
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
    outdir=options.outdir
    dosage=options.dosage


    phenotype_df, phenotype_pos_df = read_phenotype_bed(expression_bed)
    
    
    # phenotype_df =  pd.read_csv(expression_bed, sep='\t', index_col=0,header=None)
    # phenotype_df.columns = phenotype_df.iloc[0]
    # phenotype_df = phenotype_df.iloc[1: , :]
    # phenotype_df = phenotype_df.reindex(phenotype_df.index.drop(0)).reset_index(drop=True)
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

    # Replacing with simplier command
    # pr = genotypeio.PlinkReader(plink_prefix_path)
    # genotype_df = pr.load_genotypes()
    # variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    genotype_df, variant_df = genotypeio.load_genotypes(plink_prefix_path, dosages=dosage)
    os.makedirs(outdir)
    cis.map_nominal(genotype_df, variant_df,
                    phenotype_df.loc[phenotype_pos_df['chr']!='chrY'],
                    phenotype_pos_df.loc[phenotype_pos_df['chr']!='chrY'],
                    covariates_df=covariates_df,prefix='cis_nominal1',
                    output_dir=outdir, write_top=True, write_stats=True)

    #     try:
    #         # Here we have Plink1 bin,bed,fam
    #         pr = genotypeio.PlinkReader(plink_prefix_path)
    #         variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    #     except:
    #         # Here we have Plink2 psam,pgen,pvar
    #         pr = pgen.PgenReader(plink_prefix_path)
    #         variant_df2 = pr.variant_dfs
    #         variant_df = pd.DataFrame()
    #         for k1 in variant_df2.keys():
    #             dic1=pd.DataFrame(variant_df2[k1])
    #             dic1['chrom']=k1
    #             variant_df=pd.concat([variant_df,dic1])
    #         variant_df.index=variant_df.index.rename('snp')
    #         variant_df=variant_df[['chrom', 'pos']]
    #         del variant_df2
    #     genotype_df = pr.load_genotypes()
    #     Directory = './nom_output'
    #     os.mkdir(Directory)
    #     cis.map_nominal(genotype_df, variant_df,
    #                     phenotype_df.loc[phenotype_pos_df['chr']!='chrY'],
    #                     phenotype_pos_df.loc[phenotype_pos_df['chr']!='chrY'],maf_threshold=maf,
    #                     covariates_df=covariates_df,window=int(options.window),prefix='cis_nominal1',
    #                     output_dir=Directory, write_top=True, write_stats=True,run_eigenmt=True)

    
    all_files = glob.glob(f'{outdir}/cis_nominal*.parquet')
    All_Data = pd.DataFrame()
    count=0
    for bf1 in all_files:
        print(bf1)
        df = pd.read_parquet(bf1)
        df.to_csv(bf1.replace('.parquet','.tsv'),sep='\t',index=False)
        os.remove(bf1) 
        count+=1    


    try:
        cis_df = cis.map_cis(genotype_df, variant_df, 
                            phenotype_df.loc[phenotype_pos_df['chr']!='chrY'],
                            phenotype_pos_df.loc[phenotype_pos_df['chr']!='chrY'],nperm=int(options.nperm),
                            window=int(options.window),
                            covariates_df=covariates_df,maf_threshold=maf)
        print('----cis eQTLs processed ------')
        cis_df.head()
        cis_df.to_csv("Cis_eqtls.tsv",sep="\t")
        sv = ~np.isnan(cis_df['pval_beta'])
        print(f"Dropping {sum(sv)} variants withouth Beta-approximated p-values to\n.")
        cis_df_dropped = cis_df.loc[sv]
        # r = stats.pearsonr(cis_df_dropped['pval_perm'], cis_df_dropped['pval_beta'])[0]
        calculate_qvalues(cis_df_dropped, qvalue_lambda=0.85)
        cis_df_dropped.to_csv("Cis_eqtls_qval.tsv", sep='\t')
    except:
        # The beta aproximation sometimes doesnt work and results in a failure of the qtl mapping. 
        # This seems to be caused by failure to aproximate the betas
        # Hence the folowing part of the code if the above fails avoiding beta aproximation and 
        cis_df = cis.map_cis(genotype_df, variant_df, 
                            phenotype_df.loc[phenotype_pos_df['chr']!='chrY'],
                            phenotype_pos_df.loc[phenotype_pos_df['chr']!='chrY'],nperm=int(options.nperm),
                            window=int(options.window),
                            covariates_df=covariates_df,maf_threshold=maf,seed=7,beta_approx=False)
            
        print('----cis eQTLs processed ------')
        cis_df.head()
        cis_df.to_csv("Cis_eqtls.tsv",sep="\t")
        sv = ~np.isnan(cis_df['pval_beta'])
        print(f"Dropping {sum(sv)} variants withouth Beta-approximated p-values to\n.")
        cis_df_dropped = cis_df.loc[sv]
        # r = stats.pearsonr(cis_df_dropped['pval_perm'], cis_df_dropped['pval_beta'])[0]
        # calculate_qvalues(cis_df_dropped, qvalue_lambda=0.85)
        cis_df_dropped.to_csv("Cis_eqtls_qval.tsv", sep='\t')


if __name__ == '__main__':
    main()

# trans_df = tensorqtl.trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=10000,
#                            return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
# trans_df = tensorqtl.trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=5000000)
# trans_df.to_csv("Trans_eqtls.tsv",sep="\t")