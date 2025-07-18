#!/usr/bin/env python

__author__ = 'M.Ozols'
__date__ = '2023-09-18'
__version__ = '0.0.1'

import pandas as pd
import argparse
parser = argparse.ArgumentParser(
    description="""
        Replace GT filed
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-gp', '--genotype_phenotype_mapping',
    action='store',
    dest='genotype_phenotype_mapping',
    required=True,
    help=''
)

parser.add_argument(
    '-m', '--mappings',
    action='store',
    dest='mappings',
    required=True,
    help=''
)

options = parser.parse_args()

Mapping_File = pd.read_csv(options.mappings, sep='\t', dtype={'Genotype': str,'RNA':str})
Mapping_File = Mapping_File.set_index('RNA')
genotype_phenotype_mapping = pd.read_csv(options.genotype_phenotype_mapping, sep='\t', dtype={'Genotype': str,'RNA':str})
genotype_phenotype_mapping = genotype_phenotype_mapping.set_index('Genotype')
genotype_phenotype_mapping.insert(0,'Genotype','')
genotype_phenotype_mapping['Genotype'] = Mapping_File['Genotype'].astype(str)
if not list(set(genotype_phenotype_mapping['Genotype']))[0]==list(set(genotype_phenotype_mapping['Genotype']))[0]:
    # genotype_phenotype_mapping = genotype_phenotype_mapping.drop(columns='Genotype').reset_index()
    genotype_phenotype_mapping = genotype_phenotype_mapping[genotype_phenotype_mapping['Genotype']==genotype_phenotype_mapping['Genotype']]
genotype_phenotype_mapping.to_csv(f'remap_{options.genotype_phenotype_mapping}',sep='\t',index=False)
print('Done')