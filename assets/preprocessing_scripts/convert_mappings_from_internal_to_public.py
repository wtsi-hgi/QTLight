import os
import pandas as pd
Data1 = pd.read_csv('/scratch/cellfunc/mo246/eqtl_analysis/atac_eqtl_analysis/Mappings.csv')
D2 = Data1.set_index('supplier_name')
Data2= pd.read_csv('/scratch/cellfunc/mo246/eqtl_analysis/eqtl_analysis/sample_mappings.tsv','\t')
# Data2 = Data2.set_index('Genotype')
for i,idx1 in  Data2.iterrows():
    # print(idx1)
    try:
        sample1 = idx1['Genotype'].split('_')[1].replace('S','E')
        repl = D2.loc[sample1,'sanger_sample_id']
    except:
        repl=None
    print(repl)
    Data2.loc[i,'RNA']=repl
Data2 = Data2.dropna()
Data2.to_csv('/scratch/cellfunc/mo246/eqtl_analysis/atac_eqtl_analysis/sample_mappings2.tsv',index=False,sep='\t')
print('done')