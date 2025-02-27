#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
Data = pd.read_csv('plink2.rel',sep='\t',header=None)
Data_labels= pd.read_csv('plink2.rel.id',sep='\t')
try:
    Data_labels = Data_labels.rename(columns={'#IID':'IID'})
except:
    _='dif bed version'
Data_indexed = Data.rename(index=Data_labels.IID, columns=Data_labels.IID)
Data_indexed.to_csv('kinship_matrix.tsv',sep='\t')