#!/usr/bin/env python
__author__ = 'Matiss Ozols'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import argparse

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
        '-mf', '--mapping_file',
        action='store',
        dest='mapping_file',
        required=True,
        help=''
    )
    parser.add_argument(
        '-af', '--annotation_file',
        action='store',
        dest='annotation_file',
        required=True,
        help=''
    )
    parser.add_argument(
        '-ef', '--expression_file',
        action='store',
        dest='expression_file',
        required=True,
        help=''
    )

    options = parser.parse_args()
    # print("performing data encoding in BED format")
    # Load the base that determines the gene starts and ends.
    # annotation_file = "/lustre/scratch123/hgi/teams/hgi/mo11/eQTL_mapping/LIMIX/nf_core_eqtl/assets/annotation_file.txt"
    # expression_file = '/lustre/scratch123/hgi/teams/hgi/mo11/eQTL_mapping/LIMIX/work/19/b59c2c95d89532026adfc4bdef7d63/normalised_phenotype.tsv'
    # mapping_file = '/lustre/scratch123/hgi/teams/hgi/mo11/eQTL_mapping/LIMIX/work/19/b59c2c95d89532026adfc4bdef7d63/genotype_phenotype_mapping.tsv'
    annotation_file = options.annotation_file
    expression_file = options.expression_file
    mapping_file = options.mapping_file
    
    BED_Formated_Data=pd.DataFrame()
    Gene_Chr_Start_End_Data = pd.read_csv(annotation_file,sep="\t", )
    BED_Formated_Data["#chr"]=Gene_Chr_Start_End_Data.chromosome
    BED_Formated_Data["start"]=Gene_Chr_Start_End_Data.end-1
    BED_Formated_Data["end"]=Gene_Chr_Start_End_Data.end
    BED_Formated_Data["gene_id"]=Gene_Chr_Start_End_Data.feature_id
    BED_Formated_Data["idx"]=Gene_Chr_Start_End_Data.feature_id
    BED_Formated_Data=BED_Formated_Data.set_index("idx")
    BED_Formated_Data['#chr']='chr'+BED_Formated_Data['#chr'].astype(str)


    #Load the expression data and the mapping file
    Expression_Data = pd.read_csv(expression_file,sep="\t")
        #1 Conver the RNA_seq ids to HIPSci ids.
    Mapping_File=pd.read_csv(mapping_file,sep="\t")
    try:
        Mapping_File=Mapping_File.drop('Sample_Category',axis=1)
    except:
        print('does not exist')
    Mapping_File=Mapping_File.drop_duplicates().set_index("RNA")
    Mapping_File=Mapping_File.to_dict()['Genotype']
    # Expression_Data=Expression_Data.drop("qayj_3",axis=1)
    Expression_Data=Expression_Data.rename(columns=Mapping_File)
    #2 Combine the Expression data and ID mapper to get a bed format.
    mergedDf = BED_Formated_Data.merge(Expression_Data, left_index=True, right_index=True)



    chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10'
            ,'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
    
    mergedDf2 = mergedDf.set_index('#chr')
    chrs = list(set(mergedDf2.index).intersection(set(chrs)))
    mergedDf = mergedDf2.loc[chrs]
    mergedDf = mergedDf.reset_index()
    mergedDf.dropna(axis=0,inplace=True)

    # Check for missing genes - thre are missing genes, which shouldnt be there. Will have to investigate whether this is correct.
    # Genes_reference = pd.Series(Expression_Data.index)
    # Genes_reference[Genes_reference.isin(BED_Formated_Data.index.values)==False]
    # Missing_Genes = Genes_reference[Genes_reference.isin(BED_Formated_Data.index.values)==False]
    # Missing_Genes.to_csv("Missing_Genes.csv")

    # Now just save the file to BED file.
    # mergedDf["gene_id"]=mergedDf.index
    mergedDf.to_csv("Expression_Data.bed.gz", sep='\t', compression='gzip',index=False)



if __name__ == '__main__':
    # Convert files to the BED format
    main()