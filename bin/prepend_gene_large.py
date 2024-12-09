#!/usr/bin/env python

import glob
import pandas as pd
import argparse
from scipy.stats import chi2
from math import ceil

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
    d2 = d1[['gene','p.value']]
    celltype = outfile.split('__')[1]
    cell_results = []
    d2['chisq_stat'] = chi2.isf(d2['p.value'], df=1)
    gene_lambdas = d2.groupby('gene')['chisq_stat'].median() / 0.455
    cell_results.append(pd.DataFrame({
    'gene': gene_lambdas.index,
    'lambda': gene_lambdas.values,
    'cell_type': celltype
    }))
    # d2.to_csv("p_vals_lambda2.tsv", sep='\t', index=False, mode='w')
    # Iterate over the remaining files and append their content
    all_files2 = all_files[1:]
    batch_size = 20
    # Divide files into batches
    num_batches = ceil(len(all_files2) / batch_size)
    
    for batch_idx in range(num_batches):
        # batch_idx=66
        batch_files = all_files2[batch_idx * batch_size:(batch_idx + 1) * batch_size]
        batch_data = []  # Collect all files in the batch
        
        # Process each file in the batch
        for f1 in batch_files:
            f1
            try:
                d1 = pd.read_csv(f1, sep='\t')
                try:
                    name = f1.split(pattern_pre[0])[-1].split(pattern_pre[1], '')[0]
                except:
                    name = f1.split(pattern_pre[0])[-1]
                d1.insert(loc=0, column=column, value=name)
                batch_data.append(d1)
            except Exception as e:
                print(f"Issue with file {f1}: {e}")
        
        # Concatenate the batch data and write to the output file
        if batch_data:
            batch_df = pd.concat(batch_data, ignore_index=True)
            batch_df.to_csv(outfile, sep='\t', index=False, mode='a', header=(batch_idx == 0))  # Add header only for the first batch
            
            # Calculate lambdas for the batch
            d2 = batch_df[['gene', 'p.value']]
            d2['chisq_stat'] = chi2.isf(d2['p.value'], df=1)
            gene_lambdas = d2.groupby('gene')['chisq_stat'].median() / 0.455
            cell_results.append(pd.DataFrame({
                'gene': gene_lambdas.index,
                'lambda': gene_lambdas.values,
                'cell_type': celltype
            }))

    cell_results_df = pd.concat(cell_results, ignore_index=True)
    cell_results_df.to_csv('cell_lambdas.tsv', index=False,sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate files with an additional column.")
    parser.add_argument('--pattern', type=str, required=True, help='Glob pattern to match files.')
    parser.add_argument('--column', type=str, required=True, help='Name of the new column to be added.')
    parser.add_argument('--outfile', type=str, required=True, help='Output file name.')

    args = parser.parse_args()

    main(args.pattern, args.column, args.outfile)