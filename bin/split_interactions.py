#!/usr/bin/env python

import sys
import os
import pandas as pd
import argparse

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Split interaction file into separate files.
            """
    )
    parser.add_argument(
        '-i', '--input',
        action='store',
        dest='input',
        required=True,
        help=''
    )

    parser.add_argument(
        '-o', '--outdir',
        action='store',
        dest='outdir',
        required=True,
        help=''
    )

    options = parser.parse_args()
    input_file=str(options.input)
    output_dir = str(options.outdir)

    os.makedirs(output_dir, exist_ok=True)

    # Read the input file
    df = pd.read_csv(input_file, sep='\t', dtype=str)

    # Assuming the first column is sample IDs
    sample_col = df.columns[0]

    # For each interaction column write a file
    for col in df.columns[1:]:
        # Select rows where the value is not NA
        sub_df = df[[sample_col, col]].dropna(subset=[col])
        # Write to output file
        output_file = os.path.join(output_dir, f"{col}.tsv")
        sub_df.to_csv(output_file, sep='\t', index=False, header=False)
        
if __name__ == "__main__":
    main()
