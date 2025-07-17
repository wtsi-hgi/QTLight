#!/usr/bin/env python3

import pandas as pd
import glob
import os
from statsmodels.stats.multitest import multipletests

# Pattern to match chunk files in the current directory
pattern = "*cis_score.tsv.gz"

# Find all matching files
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

# Output file
out_file = "merged_cis_scores.tsv.gz"
merged_df.to_csv(out_file, sep="\t", index=False, compression="gzip")

print(f"Merged file written to: {out_file}")
