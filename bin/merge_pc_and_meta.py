#!/usr/bin/env python3

import pandas as pd
import argparse

def main(meta_file, pc_file, output_file):
    # Load input files
    meta = pd.read_csv(meta_file, sep='\t', index_col=0)
    pcs = pd.read_csv(pc_file, sep='\t', index_col=0)

    # Clean quotes from PC matrix index if present
    pcs.index = pcs.index.str.replace('"', '')

    # Merge on index (sample ID)
    merged = pcs.join(meta, how='inner')

    # One-hot encode metadata columns
    meta_encoded = pd.get_dummies(merged[meta.columns], prefix_sep='_', dtype=int)

    # Combine PCs with one-hot metadata
    final = pd.concat([merged.drop(columns=meta.columns), meta_encoded], axis=1)

    # Save result
    final.to_csv(output_file, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge and one-hot encode metadata + PC matrix")
    parser.add_argument('--meta', required=True, help='Metadata file (TSV with sample ID as first column)')
    parser.add_argument('--pc', required=True, help='PC matrix file (TSV with sample ID as first column)')
    parser.add_argument('--output', required=True, help='Output merged + one-hot encoded TSV')
    args = parser.parse_args()

    main(args.meta, args.pc, args.output)
