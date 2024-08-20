#!/usr/bin/env python

__author__ = 'Matiss Ozols and Tobi Alegbe'
__date__ = '2021-11-25'
__version__ = '0.0.1'

import pandas as pd
import scanpy as sc
import argparse
import os
import re
import gc

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge 10x data. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-agg', '--agg_columns',
        action='store',
        dest='agg_columns',
        required=True,
        help=''
    )

    parser.add_argument(
        '-h5ad', '--h5ad',
        action='store',
        dest='h5ad',
        required=True,
        help=''
    )

    options = parser.parse_args()

    h5ad = options.h5ad
    agg_columns = options.agg_columns
    agg_columns = agg_columns.split(",")

    print('Reading in data...')
    adata = sc.read_h5ad(filename=h5ad, backed='r')

    for agg_col in agg_columns:
        print(agg_col)
        print("----------")
        try:
            data_col = adata.obs[agg_col]
        except KeyError:
            print(f'Aggregation column {agg_col} doesnt exist in adata')
            continue

        for type in data_col.unique():
            print(type)
            print("----------")
            cell_adata = adata[adata.obs[agg_col] == type].copy(filename='tmp.h5ad')  # can also use .to_memory() which will be faster but more memory consuming. This is a view, not a copy
            agg_col_cleaned = re.sub(r'\W+', '_', agg_col.replace(' ', '_'))
            tp2 = re.sub(r'\W+', '_', type.replace(' ', '_'))

            output_file = f'{agg_col_cleaned}__{tp2}__split.h5ad'
            print(f'Writing to {output_file}...')

            # Write directly from the view
            print(f"Final shape is: {cell_adata.obs.shape}") 
            cell_adata.write(output_file)
            cell_adata.file.close()

            # Remove the temporary file
            os.remove('tmp.h5ad')
            # del cell_adata
            gc.collect()  # Force garbage collection to free up memory

    adata.file.close()  # Close the AnnData file to free up resources

if __name__ == '__main__':
    main()
