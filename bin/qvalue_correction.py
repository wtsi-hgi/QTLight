#!/usr/bin/env python3

# Matiss Aug 2024

import argparse
import pandas as pd
from statsmodels.stats.multitest import multipletests
import qtl_qvalue
# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Process some arguments.")
    parser.add_argument("-f", "--file", required=True, help="Input file path")
    parser.add_argument("-c", "--column", type=int, required=True, help="Column number to extract p-values")
    parser.add_argument("-n", "--new_q_val_column", required=True, help="New column name for q-values")
    parser.add_argument("-w", "--within_gene", required=True, help="Flag to save the minimum q-value result within the gene")
    
    return parser.parse_args()

# Parse the arguments
args = parse_args()

# Assign arguments to variables
file = args.file
column = args.column - 1  # Python uses 0-based indexing
new_column = args.new_q_val_column
within_gene = args.within_gene

# Load data
res = pd.read_csv(file, sep='\t')
if pd.isna(pd.to_numeric(res.iloc[0, column], errors='coerce')):
    res.columns = res.iloc[0]
    res = res[1:]

# Extract p-values
p_values = pd.to_numeric(res.iloc[:, column], errors='coerce')

# Generate q-values using Benjamini-Hochberg method
# _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
q_values = qtl_qvalue.QValue(p_values).qvalue()
# Add q-values to the DataFrame
res[new_column] = q_values

# Replace the original file with the current one
res.to_csv(file, sep='\t', index=False, header=True)

# Save a file with the variant with the minimum q-value for this gene
min_file = f"{file}_minimum_q.txt"
to_save = res[res[new_column] == res[new_column].min()]

# If little significance, may get lots of results with the minimum q-value
if len(to_save) > 1:
    # If have lots of variants with same corrected value, take the one with the smallest original p-value
    to_save = res[res.iloc[:, column] == res.iloc[:, column].min()]
    # If have variants in complete LD - take the top
    if len(to_save) > 1:
        to_save = to_save.iloc[0:1]

if within_gene == "TRUE":
    to_save.to_csv(min_file, sep='\t', index=False, header=True)

print("Performed q-value correction")
