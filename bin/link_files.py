#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import math
import argparse
import os
import subprocess

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
        '-fi', '--files_input',
        action='store',
        dest='files_input',
        required=True,
        help='the size of chunks to use'
    )
    options = parser.parse_args()


    dir1 = options.condition
    os.mkdir(dir1)
    files_input=options.files_input #'/lustre/scratch123/hgi/teams/hgi/mo11/eQTL_mapping/LIMIX/work/e9/cc9a575a6712d556e29ab034ede99e/out.txt'
    
    Files = pd.read_csv(files_input,sep='\t',header=None)
    


    for file1 in Files.iloc[:,0]:
        name = file1.split('/')[-1]
        if dir1 in name:
            print(name) 
            p = subprocess.run(f'ln -s {file1} {dir1}/{name}', shell=True, check=True)
        


if __name__ == '__main__':

    main()