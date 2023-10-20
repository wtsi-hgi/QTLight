import pandas as pd
import re

Data = pd.read_csv('/lustre/scratch123/hgi/projects/macromap/eQTL_again/premapping.csv',sep='\t')
Data_kinship = pd.read_csv('/lustre/scratch123/hgi/projects/macromap/Analysis/LIMIX/kinship_matrix.tsv',sep='\t')

Data2=pd.DataFrame()
Data2['Genotype']=Data_kinship['Unnamed: 0']
Data2['RNA']=Data2['Genotype'].str.split('-').str[1]
# Data2 = pd.read_csv('/lustre/scratch123/hgi/projects/macromap/Analysis/LIMIX/sample_mapping.tsv',sep='\t')
Data['Genotype']=''
# Data['RNA']=''

Data2['RNA'] = Data2['RNA'].str.replace('_','')

data_all = []
list_of_missing = []
for donor2 in list(set(Data['Donor line'])):
    try:
        print(donor2)
        Genotype_id = Data2[Data2['RNA'] ==donor2].Genotype.values[0]
        # Genotype_id = Data2[Data2['RNA'] ==donor2].Genotype.values[0]
        Data.loc[Data['Donor line']==donor2,'Genotype']=Genotype_id
    except:
        print(f'not available - {donor2}')
        list_of_missing.append(donor2)

Data = Data[Data['Genotype']!='']
df = pd.DataFrame({'Genotype':Data['Genotype'],'RNA':Data['Sanger Sample ID'],'Sample_Category':Data['Sample type (Ctrl or oxLDL)']})
df = df.dropna()
df =df.drop(df[df['Genotype']==''].index)
df.to_csv('/lustre/scratch123/hgi/projects/macromap/eQTL_again/sample_mapping2.tsv',index=False,sep='\t')
print('finished')