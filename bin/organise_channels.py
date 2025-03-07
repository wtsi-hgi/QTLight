#!/usr/bin/env python
import pandas as pd
import argparse
import math
from gtfparse import read_gtf
import os 

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge 10x data. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-files', '--files',
        action='store',
        dest='files',
        required=True,
        help=''
    )

    parser.add_argument(
        '-files2', '--files2',
        action='store',
        dest='files2',
        required=True,
        help=''
    )
    options = parser.parse_args()
    gp_f = pd.DataFrame(options.files2.split(" "),columns=['c1'])
    table = []
    for f1 in options.files.split(" "):
        name = '__'.join(f1.split('__')[:2])
        phen_file = os.getcwd()+'/'+ f1
        gp_fi =  os.getcwd()+'/'+ gp_f[gp_f['c1'].str.contains(name)]['c1'].values[0]
        table.append({'name':name,'phen_file':phen_file,'gp_fi':gp_fi})
    t1 = pd.DataFrame(table)
    t1.to_csv('clean_table.tsv',sep='\t',index=False)
    
main()