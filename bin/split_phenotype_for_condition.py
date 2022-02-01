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
            Filter and merge 10x data. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-c', '--condition',
        action='store',
        dest='condition',
        required=True,
        help='File containing genome annotation positions for each of the genes.'
    )    
    parser.add_argument(
        '-gp', '--genome_phenotype',
        action='store',
        dest='genome_phenotype',
        required=True,
        help='the size of chunks to use'
    )

    parser.add_argument(
        '-p', '--phenotype',
        action='store',
        dest='phenotype',
        required=True,
        help='the size of chunks to use'
    )
    options = parser.parse_args()

    print('Done')

    Phenotype = pd.read_csv(options.phenotype,sep='\t',index_col=0)
    Condition = options.condition.replace('[','').replace(']','')
    genome_phenotype = pd.read_csv(options.genome_phenotype,sep='\t')

    Samples_For_this = list(genome_phenotype[genome_phenotype['Sample_Category'] ==Condition]['RNA'])
    All_Counts = Phenotype[Samples_For_this]
    All_Counts.to_csv(f"{Condition}_phenotype.tsv",sep='\t')

if __name__ == '__main__':
    main()


