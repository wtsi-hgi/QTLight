from pandas_plink import read_plink
from bgen_reader import read_bgen
import bgen_reader
import pandas as pd
import numpy as np
import os
import pgenlib

#V0.1.1

##loader functions
def ensure_dir(file_path):
    '''Check if directory exists for output, and create it if it doesn't.'''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

##Read snp filterlist
def get_snp_df(snps_filename):
    if snps_filename:
        snp_filter_df = pd.read_csv(snps_filename,sep='\t',index_col=0)
    else:
        snp_filter_df = None
    return snp_filter_df

def get_randeff_df(randomeff_filename):
    kinship_filename = False
    readdepth_filename = False
    if ',' in randomeff_filename :
        kinship_filename, readdepth_filename = randomeff_filename.split(",")
    else :
        kinship_filename = randomeff_filename
    if kinship_filename:
        kinship_df = pd.read_csv(kinship_filename,sep='\t',index_col=0)
    else:
        kinship_df = None
    if readdepth_filename:
        readdepth_df = pd.read_csv(readdepth_filename,sep='\t',index_col=0)
    else:
        readdepth_df = None
    return kinship_df,readdepth_df

def get_samplemapping_df(sample_mapping_filename,sample_labels,key_from):
    assert(key_from in ['iid','sample'])
    if sample_mapping_filename:
        mapping_df = pd.read_csv(sample_mapping_filename,sep='\t',header=None,names=['iid','sample'])
        mapping_df.set_index(key_from,inplace=True)
    else:
        #assume the mapping is the identity mapping
        identifiers = sample_labels
        mapping_df = pd.DataFrame(data=np.vstack([identifiers,identifiers]).T,index=identifiers,columns=['iid','sample'])
    return mapping_df

def get_snp_feature_df(snp_feature_filename):
    if snp_feature_filename:
        snp_feature_df = pd.read_csv(snp_feature_filename,sep='\t')
    else :
        snp_feature_df = None
    return snp_feature_df

def get_covariate_df(covariates_filename):
    if covariates_filename:
        covariate_df = pd.read_csv(covariates_filename,sep='\t',index_col=0)
    else:
        covariate_df = None
    return covariate_df

def read_plink2(geno_prefix):
    # Paths to PLINK2 files
    pgen_file = f"{geno_prefix}.pgen"
    psam_file = f"{geno_prefix}.psam"
    pvar_file = f"{geno_prefix}.pvar"

    # Read .psam (sample data, equivalent to .fam)
    fam = pd.read_csv(psam_file, sep="\t")
    try:
        fam = fam.rename(columns={'#IID':'iid'})
    except:
        _='different version'

    # Read .pvar (variant data, equivalent to .bim)
    bim = pd.read_csv(pvar_file, sep="\t", comment="#")
    bim.columns = ['chrom', 'pos', 'snp', 'ref', 'alt', 'qual',  'info']
    # Open .pgen file for genotype data
    pgen = pgenlib.PgenReader(pgen_file.encode('utf-8'))

    return bim, fam, pgen

def get_genotype_data(geno_prefix, plinkGenotype):
    if(plinkGenotype):
        try:
            bim,fam,bed = read_plink(geno_prefix,verbose=False)
            # bim2,fam2,bed2 = read_plink("/lustre/scratch127/humgen/teams/hgi/mo11/tmp_projects127/cardinal_QTLs/ELGH/shuangs_annotations_genomeWide/pakistani/results_maf0pt01_rn_transform/genotypes/plink_genotypes_bed/plink_genotypes",verbose=False)
            
            pgen = None
        except:
            bim, fam, pgen = read_plink2(geno_prefix)
            bed = None
        fam.set_index('iid',inplace=True)
        bgen=None
    else :
        bgen = read_bgen(geno_prefix+'.bgen', verbose=False)
        bed=None
        fam =bgen['samples']
        fam = fam.to_frame("iid")
        fam.set_index('iid',inplace=True)
        
        bim = bgen['variants'].compute()
        bim = bim.assign(i = range(bim.shape[0]))
        bim['id'] = bim['rsid']
        bim = bim.rename(index=str, columns={"id": "snp"})
        bim['a1'] = bim['allele_ids'].str.split(",", expand=True)[0]
        bim.index = bim["snp"].astype(str).values
        bim.index.name = "candidate"
        
        ##Fix chromosome ids
        bim['chrom'].replace('^chr','',regex=True,inplace=True)
        bim['chrom'].replace(['X', 'Y', 'XY', 'MT'], ['23', '24', '25', '26'],inplace=True)
        ##Remove non-biallelic & non-ploidy 2 (to be sure). (These can't happen in binary plink files).
        print("Warning, the current software only supports biallelic SNPs and ploidy 2")
        bim.loc[np.logical_and(bim['nalleles']<3,bim['nalleles']>0),:]
    
    return bim,fam,bed,bgen,pgen

def get_annotation_df(anno_filename):
    annotation_col_dtypes = {'feature_id':np.object,
                         'gene_id':np.object,
                         'gene_name':np.object,
                         'chromosome':np.object,
                         'start':np.int64,
                         'end':np.int64,
                         'strand':np.object}
    annotation_df = pd.read_csv(anno_filename,sep='\t',index_col=0,dtype=annotation_col_dtypes)
    return annotation_df

def get_env_df(env_filename):
    return pd.read_csv(env_filename,sep='\t',index_col=0)

def get_phenotype_df(pheno_filename):
    return pd.read_csv(pheno_filename,sep='\t',index_col=0, na_values=['.'])

def get_grs_subset_df(grs_filename, relSnps):
    iter_csv = pd.read_csv(grs_filename, chunksize=1000, sep='\t',index_col=0, na_values=['.'])
    risk_df = None
    for chunk in iter_csv:
        if any(chunk.index.isin(relSnps)):
            risk_df = pd.concat([risk_df,chunk.reindex(labels=relSnps,axis ='index')])
    return risk_df

def get_top_qtl_results(top_qtl_results_filename):
    if top_qtl_results_filename:
        top_qtl_results_df = pd.read_csv(top_qtl_results_filename,sep='\t',index_col=0)
    else:
        top_qtl_results_df = None
    return top_qtl_results_df
