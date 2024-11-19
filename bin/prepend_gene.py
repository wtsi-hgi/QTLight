#!/usr/bin/env python

import glob
import pandas as pd
import argparse

def main(pattern, column, outfile):
    pattern_pre = pattern.split('*')
    all_files = glob.glob(f"*{pattern}")
    
    combo = pd.DataFrame()
    
    for f1 in all_files:
        d1 = pd.read_csv(f1, sep='\t')
        name = f1.replace(pattern_pre[0], '').replace(pattern_pre[1], '')
        d1.insert(loc=0, column=column, value=name)
        combo = pd.concat([combo, d1])
    
    combo.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate files with an additional column.")
    parser.add_argument('--pattern', type=str, required=True, help='Glob pattern to match files.')
    parser.add_argument('--column', type=str, required=True, help='Name of the new column to be added.')
    parser.add_argument('--outfile', type=str, required=True, help='Output file name.')

    args = parser.parse_args()

    main(args.pattern, args.column, args.outfile)