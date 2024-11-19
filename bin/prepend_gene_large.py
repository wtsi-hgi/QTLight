#!/usr/bin/env python

import glob
import pandas as pd
import argparse

def main(pattern, column, outfile):
    
    all_files = pd.DataFrame(glob.glob("*"+pattern),columns=['col1'])
    all_files = all_files[~all_files['col1'].str.contains('.index')]
    all_files = list(all_files['col1'])
    # Initialize the output file with the header from the first file
    first_file = all_files[0]
    pattern_pre = pattern.split('*')
    d1 = pd.read_csv(first_file, sep='\t')
    try:
        name = first_file.split(pattern_pre[0])[-1].split(pattern_pre[1], '')[0]
    except:
        name = first_file.split(pattern_pre[0])[-1]
    d1.insert(loc=0, column=column, value=name)
    d1.to_csv(outfile, sep='\t', index=False, mode='w')

    # Iterate over the remaining files and append their content
    for f1 in all_files[1:]:
        d1 = pd.read_csv(f1, sep='\t')
        try:
            name = f1.split(pattern_pre[0])[-1].split(pattern_pre[1], '')[0]
        except:
            name = f1.split(pattern_pre[0])[-1]
        d1.insert(loc=0, column=column, value=name)
        d1.to_csv(outfile, sep='\t', index=False, mode='a', header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate files with an additional column.")
    parser.add_argument('--pattern', type=str, required=True, help='Glob pattern to match files.')
    parser.add_argument('--column', type=str, required=True, help='Name of the new column to be added.')
    parser.add_argument('--outfile', type=str, required=True, help='Output file name.')

    args = parser.parse_args()

    main(args.pattern, args.column, args.outfile)