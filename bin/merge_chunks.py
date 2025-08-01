#!/usr/bin/env python3

"""
Author: Matiss Ozols
Date: 2025-07-17
Description: 
    This script merges all `*cis_score.tsv.gz` files in the current directory, 
    calculates FDR-adjusted q-values using the Benjamini-Hochberg procedure, 
    and writes the output to a specified file.
"""

import pandas as pd
import glob
import os
import argparse
from statsmodels.stats.multitest import multipletests

def main():
    parser = argparse.ArgumentParser(description="Merge cis_score chunks and add FDR-adjusted q-values.")
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output filename (e.g., merged_cis_scores.tsv.gz)"
    )
    args = parser.parse_args()

    # Pattern to match chunk files in the current directory
    pattern = "*cis_score.tsv.gz"
    files = sorted(glob.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No files matching pattern: {pattern}")

    print(f"Found {len(files)} files to merge.")

    # Read and concatenate all files
    dfs = []
    for f in files:
        print(f"Reading {f}")
        df = pd.read_csv(f, sep="\t", compression="gzip")
        dfs.append(df)

    merged_df = pd.concat(dfs, ignore_index=True)

    print(f"Calculating FDR-adjusted q-values for {merged_df.shape[0]} rows...")

    # Apply Benjamini-Hochberg correction
    valid = merged_df["pval_nominal"].notnull()
    qvals = pd.Series([None] * len(merged_df), index=merged_df.index)
    qvals[valid] = multipletests(merged_df.loc[valid, "pval_nominal"], method="fdr_bh")[1]
    merged_df["qval"] = qvals

    # Write output
    merged_df.to_csv(args.output, sep="\t", index=False, compression="gzip")

    print(f"Merged file written to: {args.output}")

if __name__ == "__main__":
    main()
