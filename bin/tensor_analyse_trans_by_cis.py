#!/usr/bin/env python

__author__ = 'Bradley Harris'
__date__ = '2024-01-2'
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

    # Read in the phenotype file (for this test)
    # phenotype_df_pre = pd.read_csv(phenotype_file, sep = "\t")
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_file)

    # phenotype_pos_df = pd.read_csv(phenotype_pos_file, sep = "\t", index_col=0)
    print("Loaded phenotypes")

    # Load in the cis-qval significant results
    cis_q = pd.read_csv(cis_qval_results, sep = "\t")
    if 'qval' in cis_q.columns:
        sig_variants = cis_q[cis_q['qval'] < alpha]['variant_id']
    
        if sig_variants.shape[0] > 0:
            # Load in the genotypes / dosages - Need to make sure this is adjusted to incorporate whether or not we want to limit to the variants only with a cis-effect
            genotype_df, variant_df = genotypeio.load_genotypes(plink_prefix_path, dosages = dosage)
            # Subset for the cis-only
            genotype_df = genotype_df[genotype_df.index.isin(sig_variants)]
            variant_df = variant_df[variant_df.index.isin(sig_variants)]
            print(f"Variants being tested = {genotype_df.shape[0]}")
            covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
            #phenotype_df = phenotype_df.loc[:,phenotype_df.columns.isin(covariates_df.columns)]
            #covariates_df = covariates_df.loc[:,covariates_df.columns.isin(phenotype_df.columns)]
            phenotype_df = phenotype_df[covariates_df.columns]
            # have to drop dublicate rownames. and average the repeated measures.
            phenotype_df.columns = phenotype_df.columns.str.split('.').str[0]
            covariates_df.columns = covariates_df.columns.str.split('.').str[0]
            print(f"Shape of phenotype_df:{phenotype_df.shape}")
            print(f"Shape of covariates_df:{covariates_df.shape}")
            print("Loaded genotypes, filtered genotypes and loaded covariates")

            covariates_df=covariates_df.loc[:,~covariates_df.columns.duplicated()]
            # this can be adjusted to take an average. TQTL can not account for repeated measures.
            phenotype_df=phenotype_df.loc[:,~phenotype_df.columns.duplicated()]

            covariates_df=covariates_df.T
            # Run test
            print("Running trans analysis")
            trans_df_all = trans.map_trans(genotype_df, phenotype_df.loc[phenotype_pos_df['chr']!='chrY'],
                                covariates_df = covariates_df, batch_size=10000,
                                return_sparse=True, pval_threshold=1, maf_threshold=0.05)

            # Filter the trans for distance (1Mb)

            trans_df = trans.filter_cis(trans_df_all, phenotype_pos_df, variant_df, window=window)
            print("Filtered to remove cis- effects")

            # Perform some correction of these results
            # Bonferoni correct, within gene. This corrected for the number of independent variants tested per gene (which differs slightly)
            trans_df['pval_bonf'] = trans_df.groupby('phenotype_id')['pval'].transform(lambda x: x * len(x))
            print("Bonferroni corrected")

            # Subset for top hit/gene 
            trans_df_bonf = trans_df.loc[trans_df.groupby('phenotype_id')['pval_bonf'].idxmin()]

            ## Subset for significant effects within gene
            #trans_df_bonf = trans_df.loc[trans_df['pval_bonf'] < alpha]

            # Ceiling the bonf before fdr
            trans_df_bonf['pval_bonf'] = trans_df_bonf['pval_bonf'].apply(lambda x: 1 if x > 1 else x)

            # FDR across genes (This ccounts for the number of genes tested and allows results to be comparable across the conditions with different nGenes)
            trans_df_bonf['pval_bonf_fdr'] = fdrcorrection(trans_df_bonf['pval_bonf'])[1]
            # Sort
            trans_df_bonf_sorted = trans_df_bonf.sort_values(by='pval_bonf_fdr')
            print("FDR corrected and sorted")

            # Save the bonferoni/FDR corrected results
            trans_df_bonf_sorted.to_csv(f"{outdir}/trans-by-cis_bonf_fdr.tsv", sep = "\t", index=False)
            print("Saved the corrected results")
            
            # Save all results
            #print("Saving all trans results")
            #trans_df.to_csv(f"{outdir}/trans-by-cis_all.tsv.gz", compression='gzip', sep = "\t", index=False)
    
    else:
       print("No significant variants at this nPC") 

if __name__ == '__main__':
    main()