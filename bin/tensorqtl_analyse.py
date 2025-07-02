#!/usr/bin/env python

__author__ = 'Matiss Ozols and Hannes Ponstingl'
__date__ = '2021-11-25'
__version__ = '0.0.1'

# https://github.com/broadinstitute/tensorqtl
# https://zenodo.org/record/4118403#.YHclzGMo9TY

import torch
import pandas as pd
import sys
import tensorqtl
from tensorqtl import read_phenotype_bed, cis, calculate_qvalues,pgen,trans
import threading
import numpy as np
import scipy.stats as stats
import glob    
import argparse
import os
import genotypeio
import gzip
class BackgroundGenerator(threading.Thread):
    # Adapted from https://github.com/justheuristic/prefetch_generator
    def __init__(self, generator, max_prefetch=10):
        threading.Thread.__init__(self)
        self.queue = queue.Queue(max_prefetch)
        self.generator = generator
        self.daemon = True
        self.start()

    def run(self):
        try:
            for item in self.generator:
                self.queue.put(item)
        except Exception as exception:
            self.queue.put(exception)
        self.queue.put(None)

    def next(self):
        next_item = self.queue.get()
        if next_item is None:
            self.join()
            raise StopIteration
        if isinstance(next_item, Exception):
            self.join()
            raise next_item
        return next_item

    def __next__(self):
        return self.next()
    def __iter__(self):
        return self
    

def gwas_trans_mapping(batch_size = 1000,chrom_to_map =2,variant_df=None,
                               phenotype_df=None,phenotype_pos_df=None,genotype_df=None,
                               phenotype_df1=None,covariates_df=None,outdir=None,map_nominal=None):
    
    chrom_to_map = str(chrom_to_map)
    # filter down to only chr20 for calibration purposes to get all the nominal values for each gene agains the chr20 variants
    sig_variants = list(variant_df[variant_df['chrom']==chrom_to_map].index)
    genotype_df_trans = genotype_df[genotype_df.index.isin(sig_variants)]
    variant_df_trans = variant_df[variant_df.index.isin(sig_variants)]
    
    chr20_genes = list(phenotype_pos_df[phenotype_pos_df['chr']==chrom_to_map].index)
    phenotype_batch_chr20 = phenotype_df.loc[chr20_genes]
    phenotype_pos_batch_chr20 = phenotype_pos_df.loc[phenotype_batch_chr20.index]
    
    chr20_genes = list(phenotype_pos_df[phenotype_pos_df['chr']!=chrom_to_map].index)  
    phenotype_batch_NOTchr20 = phenotype_df.loc[chr20_genes]
    phenotype_pos_batch_NOTchr20 = phenotype_pos_df.loc[phenotype_batch_NOTchr20.index]
    # Open the output file for writing
    with gzip.open(f"{outdir}/trans_all__{chrom_to_map}__genes_on_{chrom_to_map}.tsv.gz", "wt") as outfile:
        header_written = False
        # Divide the phenotype dataframe into batches of 10 genes
        # I think this is asking whether these 10 genes are having any trans effects, and if nothing is available it doesnt give the results.
        for i in range(0, len(phenotype_batch_chr20), batch_size):
            # Subset the phenotype data for the current batch of 10 genes
            phenotype_batch_df = phenotype_batch_chr20.iloc[i:i + batch_size]
            phenotype_pos_batch_df = phenotype_pos_batch_chr20.loc[phenotype_batch_df.index]
            # With this we are producing all SNPs against all Genes on a chromosome defined in the that the gene is on. 
            # For the purpose of producing nominal values across all chromosomes we could modify the phenotype pos_df to the chromosome of interest and do this iteratevely.
            cis.map_nominal(genotype_df_trans, variant_df_trans,
                            phenotype_batch_df,
                            phenotype_pos_batch_df,window=3000000000,
                            covariates_df=covariates_df,prefix='cis_nominal1_30',
                            output_dir=outdir, write_top=True, write_stats=True)
            all_files = glob.glob(f'{outdir}/cis_nominal1_30*.parquet')
            All_Data = pd.DataFrame()
            for bf1 in all_files:
                bf1 = all_files[0]
                df = pd.read_parquet(bf1)
                # Append the results to the file incrementally
                if not header_written:
                    df.to_csv(outfile, sep="\t", index=False, header=True)
                    header_written = True
                else:
                    df.to_csv(outfile, sep="\t", index=False, header=False)
                os.remove(bf1) 
            print(f"Processed batch {i // batch_size + 1}")       
    
    for chr_to in set(phenotype_pos_batch_NOTchr20['chr']):  
        # Open the output file for appending    
        chr20_genes__perChr = list(phenotype_pos_batch_NOTchr20[phenotype_pos_batch_NOTchr20['chr']==chr_to].index) 
        phenotype_pos_batch_NOTchr20__perChr  = phenotype_pos_batch_NOTchr20.loc[chr20_genes__perChr]
        phenotype_batch_NOTchr20__perChr = phenotype_batch_NOTchr20.loc[phenotype_pos_batch_NOTchr20__perChr.index]
        
        with gzip.open(f"{outdir}/trans_all__{chrom_to_map}__genes_on_{chr_to}.tsv.gz", "a") as outfile:
            header_written = False
            # Divide the phenotype dataframe into batches of 10 genes
            # I think this is asking whether these 10 genes are having any trans effects, and if nothing is available it doesnt give the results.
            for i in range(0, len(phenotype_batch_NOTchr20__perChr), batch_size):
                # Subset the phenotype data for the current batch of 10 genes
                phenotype_batch_df = phenotype_batch_NOTchr20__perChr.iloc[i:i + batch_size]
                phenotype_pos_batch_df = phenotype_pos_batch_NOTchr20__perChr.loc[phenotype_batch_df.index]
                phenotype_pos_batch_df['chr2']=phenotype_pos_batch_df['chr']
                phenotype_pos_batch_df['chr']=chrom_to_map
                # With this we are producing all SNPs against all Genes on a chromosome defined in the that the gene is on. 
                # For the purpose of producing nominal values across all chromosomes we could modify the phenotype pos_df to the chromosome of interest and do this iteratevely.
                cis.map_nominal(genotype_df_trans, variant_df_trans,
                                phenotype_batch_df,
                                phenotype_pos_batch_df,window=3000000000,
                                covariates_df=covariates_df,prefix='cis_nominal1_30',
                                output_dir=outdir, write_top=True, write_stats=True)
                all_files = glob.glob(f'{outdir}/cis_nominal1_30*.parquet')
                All_Data = pd.DataFrame()
                for bf1 in all_files:
                    bf1 = all_files[0]
                    df = pd.read_parquet(bf1)
                    df.index = df['phenotype_id']
                    df['chr2']= phenotype_pos_batch_df['chr2']
                    df['start_distance']= df['chr2']+'_vs_'+chrom_to_map
                    del df['chr2']
                    # Append the results to the file incrementally
                    if not header_written:
                        df.to_csv(outfile, sep="\t", index=False, header=True)
                        header_written = True
                    else:
                        df.to_csv(outfile, sep="\t", index=False, header=False)
                    os.remove(bf1) 
                print(f"Processed batch {i // batch_size + 1}")    
    

def main():

    print('PyTorch {}'.format(torch.__version__))
    print('Pandas {}'.format(pd.__version__))
    print('Tensorqtl {}'.format(tensorqtl.__version__))


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

    parser.add_argument(
        '-chmt', '--chrom_to_map_trans',
        action='store',
        dest='chrom_to_map_trans',
        required=False,
        default=False,
        help='can contain a str such as 2,3,5 - indicating chromosomes against which we want to run gwas analysis'
    )
    
    parser.add_argument(
        '-nom', '--map_nominal',
        action='store_true',
        dest='map_nominal',
        default=False,
        help=''
    )
    
    parser.add_argument(
        '-ind', '--map_independent_qtls',
        action='store_true',
        dest='map_independent_qtls',
        default=False,
        help=''
    )

    options = parser.parse_args()
    maf=float(options.maf)
    map_nominal=options.map_nominal
    
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
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
    covariates_df = covariates_df[list(set(phenotype_df.columns).intersection(set(covariates_df.columns)))]
    phenotype_df = phenotype_df[covariates_df.columns]
    print(phenotype_df.head())
    # have to drop dublicate rownames. and average the repeated measures.
    # phenotype_df.columns = phenotype_df.columns.str.split('.').str[0]
    # covariates_df.columns = covariates_df.columns.str.split('.').str[0]

    covariates_df=covariates_df.loc[:,~covariates_df.columns.duplicated()]
    # this can be adjusted to take an average. TQTL can not account for repeated measures.
    phenotype_df=phenotype_df.loc[:,~phenotype_df.columns.duplicated()]

    
    covariates_df=covariates_df.T
    phenotype_df_genes = list(set(phenotype_pos_df.index))

    print('----Fine read ------')
    if torch.cuda.is_available():
        print(f'  * using GPU ({torch.cuda.get_device_name(torch.cuda.current_device())})')
    else:
        print('  * WARNING: using CPU!')

    genotype_df, variant_df = genotypeio.load_genotypes(plink_prefix_path, dosages=dosage)
    try:
        os.makedirs(outdir)
    except:
        print('exist')

    print(f"Dropping {len(set(phenotype_df.columns)) - len(set(phenotype_df.columns).intersection(set(genotype_df.columns)))} phenotypes because no genotype is present for these")
    phenotype_df = phenotype_df[list(set(phenotype_df.columns).intersection(set(genotype_df.columns)))]
    covariates_df = covariates_df.loc[list(set(phenotype_df.columns).intersection(set(genotype_df.columns))),:]
    covariates_df = covariates_df.sort_index()
    # Make sure they are always sorted the same regardless of what run it is.
    phenotype_df = phenotype_df.loc[phenotype_df_genes,sorted(phenotype_df.columns, reverse=True)]
    covariates_df = covariates_df.loc[phenotype_df.columns]
    
    if options.chrom_to_map_trans:
        print("Running trans analysis")
        for chr1 in options.chrom_to_map_trans.split(','):
            gwas_trans_mapping(chrom_to_map = chr1,variant_df=variant_df,
                               phenotype_df=phenotype_df,phenotype_pos_df=phenotype_pos_df,genotype_df=genotype_df,
                               phenotype_df1=phenotype_df_genes,covariates_df=covariates_df,outdir=outdir,map_nominal=map_nominal) ## Map all the genes across genome against all SNPs indicated as chrom_to_map parameter. For example if 2 then wil map all the genes across chromosomes against chr2 NPs.
    else:    
        if map_nominal:
            cis.map_nominal(genotype_df, variant_df,
                            phenotype_df.loc[phenotype_df_genes],
                            phenotype_pos_df.loc[phenotype_df_genes],maf_threshold=maf,
                            covariates_df=covariates_df.loc[phenotype_df.columns],prefix='cis_nominal1',
                            output_dir=outdir, write_top=map_nominal, write_stats=map_nominal)
        
            all_files = glob.glob(f'{outdir}/cis_nominal*.parquet')
            All_Data = pd.DataFrame()
            count=0
            for bf1 in all_files:
                # print(bf1)
                df = pd.read_parquet(bf1)
                df.to_csv(bf1.replace('.parquet','.tsv'),sep='\t',index=False)
                os.remove(bf1) 
                count+=1    
  
    phenotype_df = phenotype_df[~phenotype_df.index.duplicated(keep="first")]

    phenotype_pos_df = phenotype_pos_df[~phenotype_pos_df.index.duplicated(keep="first")]
    try:

        cis_df = cis.map_cis(genotype_df, variant_df, 
                            phenotype_df.loc[phenotype_df_genes],
                            phenotype_pos_df.loc[phenotype_df_genes],nperm=int(options.nperm),
                            window=int(options.window),
                            covariates_df=covariates_df.loc[phenotype_df.columns],maf_threshold=maf,seed=7)
        

        print('----cis eQTLs processed ------')
        cis_df.head()
        cis_df.to_csv(f"{outdir}/Cis_eqtls.tsv",sep="\t")
        sv = ~np.isnan(cis_df['pval_beta'])
        print(f"Dropping {sum(sv)} variants withouth Beta-approximated p-values to\n.")
        cis_df_dropped = cis_df.loc[sv]
        # r = stats.pearsonr(cis_df_dropped['pval_perm'], cis_df_dropped['pval_beta'])[0]
        calculate_qvalues(cis_df_dropped, qvalue_lambda=0.85)
        cis_df_dropped.to_csv(f"{outdir}/Cis_eqtls_qval.tsv", sep='\t')


    except:
        # The beta aproximation sometimes doesnt work and results in a failure of the qtl mapping. 
        # This seems to be caused by failure to aproximate the betas
        # Hence the folowing part of the code if the above fails avoiding beta aproximation and 
        print('----cis eQTLs failed to aproximate betas ------')
        #import pdb; pdb.set_trace()
        cis_df = cis.map_cis(genotype_df, variant_df, 
                            phenotype_df.loc[phenotype_df_genes],
                            phenotype_pos_df.loc[phenotype_df_genes],nperm=int(options.nperm),
                            window=int(options.window),
                            covariates_df=covariates_df.loc[phenotype_df.columns],maf_threshold=maf,seed=7,beta_approx=False)
        print('----cis eQTLs processed ------')
        cis_df.head()
        cis_df.to_csv(f"{outdir}/Cis_eqtls.tsv",sep="\t")
        sv = ~np.isnan(cis_df['pval_beta'])
        print(f"Dropping {sum(sv)} variants withouth Beta-approximated p-values to\n.")
        cis_df_dropped = cis_df.loc[sv]
        # r = stats.pearsonr(cis_df_dropped['pval_perm'], cis_df_dropped['pval_beta'])[0]
        # calculate_qvalues(cis_df_dropped, qvalue_lambda=0.85)
        # Perform conditional analysis
    #######################
    if options.map_independent_qtls:        
        try:
            indep_df = cis.map_independent(genotype_df, variant_df, cis_df_dropped,
                                            phenotype_df.loc[phenotype_df_genes],       
                                            phenotype_pos_df.loc[phenotype_df_genes],
                                            nperm=int(options.nperm), window=int(options.window),
                                            covariates_df=covariates_df,maf_threshold=maf,seed=7)
            indep_df.to_csv(f"{outdir}/Cis_eqtls_independent.tsv",sep="\t",index=False)
        except:
            print("No significant phenotypes for cis.map_independent")
        
    else:
        print('--Skipping map_independent analysis')
    
    if 'qvals' not in cis_df_dropped.columns:
        # Add 'qvals' column with None values
        cis_df_dropped['qvals'] = None
    cis_df_dropped.to_csv(f"{outdir}/Cis_eqtls_qval.tsv", sep='\t')


if __name__ == '__main__':
    main()

# trans_df = tensorqtl.trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=10000,
#                            return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
# trans_df = tensorqtl.trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=5000000)
# trans_df.to_csv("Trans_eqtls.tsv",sep="\t")