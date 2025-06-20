#!/usr/bin/env python

__author__ = 'Bradley Harris'
__date__ = '2024-02-16'
__version__ = '0.0.1'

# https://github.com/broadinstitute/tensorqtl
# https://zenodo.org/record/4118403#.YHclzGMo9TY

print("Loading base packages")
import numpy as np
import scipy.stats as stats
import glob
def main():
    print("Starting main")
    import torch
    import pandas as pd
    import tensorqtl
    from tensorqtl import read_phenotype_bed, genotypeio, cis, trans, calculate_qvalues,pgen 
    print('PyTorch {}'.format(torch.__version__))
    print('Pandas {}'.format(pd.__version__))
    print('Tensorqtl {}'.format(tensorqtl.__version__))
    import argparse
    import os
    from statsmodels.stats.multitest import fdrcorrection
    os.system('python -V')
    os.system('which python')
    """Run CLI."""
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge 10x data. Save to AnnData object.
            """
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
        '-pf', '--phenotype_file',
        action='store',
        dest='phenotype_file',
        required=True,
        help=''
    )
    
    # parser.add_argument(
    #     '-ppf', '--phenotype_pos_file',
    #     action='store',
    #     dest='phenotype_pos_file',
    #     required=True,
    #     help=''
    # )
    
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
    
    parser.add_argument(
        '-pt', '--pval_threshold',
        action='store',
        dest='pval_threshold',
        required=False,
        default=0.01,
        help=''
    )

    parser.add_argument(
        '-cq', '--cis_qval_results',
        action='store',
        dest='cis_qval_results',
        required=False,
        help=''
    )

    parser.add_argument(
        '-a', '--alpha',
        action='store',
        dest='alpha',
        required=False,
        default=0.05,
        help=''
    )

    # Test, so define the options here [Using secretory normalised expression bed so doesn't match subset of tuft eQTLs - fine for testing though]
    # os.chdir("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/eqtl_tobi/temp_trans")
    # covariates_file="Covariates.tsv"
    # expression_bed="Expression_Data.bed.gz" # This needs to be a bed file not a tsv
    # plink="Secretory_plink_genotypes"
    # cis_qval_results="Cis_eqtls_qval.tsv"
    # alpha = 0.01
    # optim_npheno=32

    # Get the script args
    options = parser.parse_args()
    covariates_file=options.covariates_file
    print(f"Covariates: {covariates_file}")
    window=float(options.window)
    print(f"window: {str(window)}")
    phenotype_file=options.phenotype_file
    print(f"phenotype_file: {phenotype_file}")
    # phenotype_pos_file=options.phenotype_pos_file
    # print(f"phenotype_pos_file: {phenotype_pos_file}")
    plink_prefix_path=options.plink_prefix_path
    print(f"plink_prefix_path: {plink_prefix_path}")
    outdir=options.outdir
    print(f"outdir: {outdir}")
    dosage=options.dosage
    print(f"dosage: {dosage}")
    maf=float(options.maf)
    print(f"maf: {str(maf)}")
    cis_qval_results=options.cis_qval_results
    print(f"cis_qval_results: {cis_qval_results}")
    alpha = float(options.alpha)
    print(f"alpha: {str(alpha)}")
    pval_threshold=float(options.pval_threshold)
    # Read in the phenotype file (for this test)
    # phenotype_df_pre = pd.read_csv(phenotype_file, sep = "\t")
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_file)

    # phenotype_pos_df = pd.read_csv(phenotype_pos_file, sep = "\t", index_col=0)
    print("Loaded phenotypes")

    # Load in the cis-qval significant results
    cis_q = pd.read_csv(cis_qval_results, sep = "\t")
    if 'qval' in cis_q.columns:
        sig_genes = cis_q[cis_q['qval'] < alpha]['phenotype_id']
    
        if sig_genes.shape[0] > 0:
            # Load in the genotypes / dosages - Need to make sure this is adjusted to incorporate whether or not we want to limit to the variants only with a cis-effect
            genotype_df, variant_df = genotypeio.load_genotypes(plink_prefix_path, dosages = dosage)
            print(f"Genes being tested = {len(sig_genes)}")
            covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
            covariates_df = covariates_df[list(set(phenotype_df.columns).intersection(set(covariates_df.columns)))]
            #phenotype_df = phenotype_df.loc[:,phenotype_df.columns.isin(covariates_df.columns)]
            #covariates_df = covariates_df.loc[:,covariates_df.columns.isin(phenotype_df.columns)]
            phenotype_df = phenotype_df[covariates_df.columns]
            # Subset for those with cis effects
            phenotype_df = phenotype_df[phenotype_df.index.isin(sig_genes)]
            phenotype_pos_df = phenotype_pos_df[phenotype_pos_df.index.isin(sig_genes)]
            # have to drop duplicate rownames. and average the repeated measures.
            # phenotype_df.columns = phenotype_df.columns.str.split('.').str[0]
            # covariates_df.columns = covariates_df.columns.str.split('.').str[0]
            print(f"Shape of phenotype_df:{phenotype_df.shape}")
            print(f"Shape of covariates_df:{covariates_df.shape}")
            print("Loaded genotypes, filtered genotypes and loaded covariates")

            covariates_df=covariates_df.loc[:,~covariates_df.columns.duplicated()]
            # this can be adjusted to take an average. TQTL can not account for repeated measures.
            phenotype_df=phenotype_df.loc[:,~phenotype_df.columns.duplicated()]

            covariates_df=covariates_df.T
            
            
            covariates_df = covariates_df.sort_index()
            # Make sure they are always sorted the same regardless of what run it is.
            phenotype_df = phenotype_df.loc[:,sorted(phenotype_df.columns, reverse=True)]
    
            # Run test
            print("Running trans analysis")
            trans_df_all = trans.map_trans(genotype_df, phenotype_df,
                                covariates_df = covariates_df.loc[phenotype_df.columns], batch_size=10000,
                                return_sparse=True, pval_threshold=pval_threshold, maf_threshold=maf)

            # Filter the trans for distance (1Mb)
            trans_df = trans.filter_cis(trans_df_all, phenotype_pos_df, variant_df, window=window)
            print("Filtered to remove cis- effects")

            # Calculate what the p_value thresh is and add in a new column
            pval_thresh = 5e-8/len(sig_genes)
            trans_df['genome_wide_mt_p_thresh'] = pval_thresh
            trans_df['genome_wide_mt_sig'] = trans_df['pval'] < pval_thresh
            print(f"Added pvalue thresh and annotation of significance at level: {pval_thresh}")

            # Sort and save
            trans_df = trans_df.sort_values(by='pval')
            trans_df.to_csv(f"{outdir}/trans-of-cis_all.tsv", sep = "\t", index=False)
            print("Saved the corrected results")
            
            # Also save filtered
            trans_df_filt = trans_df[trans_df['pval'] < trans_df['genome_wide_mt_p_thresh']]
            trans_df_filt.to_csv(f"{outdir}/trans-of-cis_filt.tsv", sep = "\t", index=False)
            print("Saved the corrected results")
            
            # Save all results
            #print("Saving all trans results")
            #trans_df.to_csv(f"{outdir}/trans-by-cis_all.tsv.gz", compression='gzip', sep = "\t", index=False)
    
    else:
       print("No significant genes at this nPC") 

if __name__ == '__main__':
    main()