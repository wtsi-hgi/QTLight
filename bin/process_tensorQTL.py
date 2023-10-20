#!/usr/bin/env python3
import os
from shutil import copyfile,copytree,copy
from os import listdir
import glob
import argparse
import pandas as pd

__date__ = '2022-21-03'
__version__ = '0.0.2'
input_dir = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/eqtl/franke_repeat/results/TensorQTL_eQTLS'
out_dir = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Franke_with_genotypes_nfCore2/results2/handover/Summary_plots/Franke_with_genotypes_nfCore2/eQTLs___801'


# List all the folders in the directory
dir_files = [x[0] for x in os.walk(input_dir)]
try:
    os.mkdir(f'{out_dir}')
except:
    print('exists')
# Loop through each
for dir1 in dir_files[1:]:
    print(dir1)
    cell_type = dir1.split('/')[-1].split('-')[0]
    aggregation_method = dir1.split('/')[-1].split('-')[1]
    Data = pd.read_csv(f'{dir1}/Cis_eqtls_qval.tsv',sep='\t')
    Data = Data[Data.qval<=0.05]
    if (len(Data)>0):
        Data = Data.sort_values(by='qval')
        try:
            os.mkdir(f'{out_dir}/{cell_type}')
        except:
            print('exists')
        Data.to_csv(f'{out_dir}/{cell_type}/{aggregation_method}_TQ__Cis_eqtls_qval.tsv',sep='\t',index=False)
    else:
        print('skipp')
    

# Open Each Cis_eqtls file

# Filter based on the qval

# If any then record in a summary folder


