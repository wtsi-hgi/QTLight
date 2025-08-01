#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2025-07-16'
__version__ = '0.0.2'

import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(
        description="Split BED file by 4th column (feature IDs) into equal-sized chunks."
    )

    parser.add_argument(
        '-bed', '--bed_file',
        required=True,
        help='BED file with features in 4th column.'
    )

    parser.add_argument(
        '-cs', '--chunk_size',
        required=True,
        type=int,
        help='Number of features per chunk.'
    )

    parser.add_argument(
        '-o', '--output_prefix',
        default='Chunking_file',
        help='Prefix for output chunk files.'
    )

    options = parser.parse_args()

    # Read BED file (assumes standard BED format: chrom, start, end, name, ...)
    bed_df = pd.read_csv(options.bed_file, sep='\t')
    
    if bed_df.shape[1] < 4:
        raise ValueError("BED file must have at least 4 columns.")

    feature_ids = bed_df.iloc[:,3].dropna().unique()  # 4th column = feature ID
    feature_ids = sorted(feature_ids)  # optional: sort for consistency

    # Chunk into lists
    chunks = [feature_ids[i:i + options.chunk_size] for i in range(0, len(feature_ids), options.chunk_size)]

    # Output each chunk
    for i, chunk in enumerate(chunks, start=1):
        out_df = pd.DataFrame({'feature_id': chunk})
        out_file = f'{options.output_prefix}__chunk_{i}.tsv'
        out_df.to_csv(out_file, sep='\t', index=False, header=True)
    print(f"Done. Created {len(chunks)} chunks.")

if __name__ == '__main__':
    main()
