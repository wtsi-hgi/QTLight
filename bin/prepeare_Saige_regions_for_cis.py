#!/usr/bin/env python
__author__ = 'Matiss Ozols'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import argparse
import math
from gtfparse import read_gtf



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
        '-af', '--annotation_file',
        action='store',
        dest='annotation_file',
        required=True,
        help=''
    )
    parser.add_argument(
        '-genes', '--genes',
        action='store',
        dest='genes',
        required=True,
        help=''
    )

    parser.add_argument(
        '-gtf', '--gtf_type',
        action='store',
        dest='gtf_type',
        required=True,
        help='gtf_type'
    )

    parser.add_argument(
        '-wd', '--window',
        action='store',
        dest='window',
        required=True,
        help='window +- from tss region to test the QTLs'
    )



    options = parser.parse_args()
    
    annotation_file = options.annotation_file
    expression_file = options.genes
    window = int(options.window)
    gtf_type = options.gtf_type

    Expression_Data = pd.read_csv(expression_file,sep="\t",header=None)
    if(len(Expression_Data[0].iloc[0].split('.'))>1):
        prot_version=True
        Expression_Data = Expression_Data[0].str.split('.').str[0]
    else:
        prot_version=False
    
    try:
        df = read_gtf(annotation_file)
        if (gtf_type=='gene'):
            df2= df.filter(df["feature"] == "gene").to_pandas()
            Gene_Chr_Start_End_Data =df2[['gene_id','start','end','strand','seqname']]
        elif (gtf_type=='transcript'):
            df2= df.filter(df["feature"] == "transcript").to_pandas()
            Gene_Chr_Start_End_Data =df2[['transcript_id','start','end','strand','seqname']]
        else:
            _ = 'you havent specified which type of analysis you are performing'
            # 
        
    except:
        Gene_Chr_Start_End_Data = pd.read_csv(annotation_file,index_col=None,sep='\t')

    if (gtf_type=='gene'):
        Gene_Chr_Start_End_Data.rename(columns={'gene_id':'feature_id','seqname':'chromosome'},inplace=True)
    elif (gtf_type=='transcript'):
        Gene_Chr_Start_End_Data.rename(columns={'transcript_id':'feature_id','seqname':'chromosome'},inplace=True)
    else:
        _ = 'you havent specified which type of analysis you are performing'

        
    Gene_Chr_Start_End_Data.drop_duplicates(inplace=True)
    Gene_Chr_Start_End_Data=Gene_Chr_Start_End_Data.set_index('feature_id')
    f = list(Expression_Data[0])
    f2 = set(Gene_Chr_Start_End_Data.index).intersection(set(f))
    Gene_Chr_Start_End_Data = Gene_Chr_Start_End_Data.loc[list(f2)]   

    BED_Formated_Data=pd.DataFrame()
    BED_Formated_Data["#chr"]=Gene_Chr_Start_End_Data.chromosome
    
    BED_Formated_Data["start"]=0
    BED_Formated_Data["end"]=0
    
    # + strand
    idx1 = Gene_Chr_Start_End_Data[Gene_Chr_Start_End_Data.strand =='+'].index
    BED_Formated_Data.loc[idx1,"start"]=Gene_Chr_Start_End_Data.loc[idx1,"start"]-1
    BED_Formated_Data.loc[idx1,"end"]=Gene_Chr_Start_End_Data.loc[idx1,"start"]

    # - strand
    idx1 = Gene_Chr_Start_End_Data[Gene_Chr_Start_End_Data.strand =='-'].index
    BED_Formated_Data.loc[idx1,"start"]=Gene_Chr_Start_End_Data.loc[idx1,"end"]-1
    BED_Formated_Data.loc[idx1,"end"]=Gene_Chr_Start_End_Data.loc[idx1,"end"]
    
    
    BED_Formated_Data["start"]=Gene_Chr_Start_End_Data.start+((Gene_Chr_Start_End_Data.end-Gene_Chr_Start_End_Data.start)/2).apply(math.ceil)
    BED_Formated_Data["end"]=BED_Formated_Data["start"]+1

    # window=10000
    BED_Formated_Data["start"] = BED_Formated_Data["start"]-window
    BED_Formated_Data["end"] = BED_Formated_Data["end"]+window

    # BED_Formated_Data["gene_id"]=BED_Formated_Data.index
    BED_Formated_Data = BED_Formated_Data.sort_values(by=['#chr','start','end'])
    # BED_Formated_Data["idx"]=Gene_Chr_Start_End_Data.feature_id
    # BED_Formated_Data=BED_Formated_Data.set_index("idx")
    BED_Formated_Data['#chr']=BED_Formated_Data['#chr'].astype(str)    
    BED_Formated_Data.to_csv("gene_regions_to_test.tsv", sep='\t', header=False)
    print('Done')     


if __name__ == '__main__':
    main()