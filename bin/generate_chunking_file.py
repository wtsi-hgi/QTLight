#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import math
import argparse
import os

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
        '-covar', '--covar_file',
        action='store',
        dest='covar_file',
        required=True,
        help=''
    )   

    parser.add_argument(
        '-condition', '--condition',
        action='store',
        dest='condition',
        required=True,
        help=''
    )   

    parser.add_argument(
        '-ga', '--genome_annotation',
        action='store',
        dest='genome_annotation',
        required=True,
        help='File containing genome annotation positions for each of the genes.'
    )    

    parser.add_argument(
        '-pf', '--phenotype_file',
        action='store',
        dest='phenotype_file',
        required=True,
        help='Phenotype file to subset annotation file, to limit number of jobs generated.'
    )    

    parser.add_argument(
        '-cs', '--chunk_size',
        action='store',
        dest='chunk_size',
        required=True,
        help='the size of chunks to use'
    )
    options = parser.parse_args()



    genes = int(options.chunk_size)
    Data = pd.read_csv(options.genome_annotation,sep='\t')

    Phenotype_file = pd.read_csv(options.phenotype_file,sep='\t',index_col=0)
    Phenotype_file_index = list(Phenotype_file.index)

    Phenotype_file_index
    Data2=Data.set_index('feature_id')
   
    idx_overlap = list(set(Data2.index).intersection(set(Phenotype_file_index)))
    Data =Data2.loc[idx_overlap]
    Data=Data.reset_index()
    chroms = [*range(1,23,1)]
    chroms.extend(['X','Y','MT'])
    data_to_export=[]

    for chr in chroms:    
        chr=str(chr)
        chr_genes = Data[Data.chromosome == chr]
        chr_genes = chr_genes.sort_values('start')
        chr_genes = chr_genes.reset_index()
        nr_total_genes = chr_genes.shape[0]
        nr_itterations = int(nr_total_genes/genes)
        for i in range(0,math.ceil(nr_itterations)+1):
            start_pos = i*genes
            end_pos=start_pos+genes

        
            if end_pos>nr_total_genes:
                end_pos=nr_total_genes
            if start_pos==end_pos:
                continue
            chr_genes_set = chr_genes[start_pos:end_pos+1]
            st_st=list(chr_genes_set.start)+list(chr_genes_set.end)
            min_pos = min(st_st)
            max_pos = max(st_st)
            data_to_export.append(f'{chr}:{min_pos}-{max_pos}')
            
    data_to_export = pd.DataFrame(data_to_export)
    data_to_export.to_csv('Chunging_file.tsv',header=None,index=None)
    
    data_to_export2=data_to_export.rename(columns={0:'Range'})
    data_to_export2['condition']=options.condition
    data_to_export2['phenotypeFile']=f"{os.getcwd()}/normalised_phenotype.tsv"
    data_to_export2['covars']=f"{os.getcwd()}/pcs.tsv"
    data_to_export2['anotation_file']=f"{os.getcwd()}/annotation_file_processed.tsv"
    
    data_to_export2.to_csv('limix_chunking.tsv',sep='\t',index=None)
    
    Data.to_csv('annotation_file_processed.tsv',index=None,sep='\t')
    print('Done')


if __name__ == '__main__':
    main()
