#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import math
import argparse

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter the sample mapping file to only contain the required columns for limix.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-gp', '--genotype_phenotype',
        action='store',
        dest='genotype_phenotype',
        required=True,
        help='Genotype Phenotype file with extra entries'
    )    
    options = parser.parse_args()
    Data = pd.read_csv(options.genotype_phenotype,sep='\t')
    
    Data2=Data[['Genotype','RNA']]
    Data2.to_csv('genotype_phenotype.tsv',index=False,sep='\t')



if __name__ == '__main__':
    main()
