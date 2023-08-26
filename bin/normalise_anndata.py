#!/usr/bin/env python

__author__ = 'Tobi Alegbe'
__date__ = '2023-04-27'
__version__ = '0.0.1'


import scanpy as sc
import pandas as pd
import numpy as np
# from pysctransform import vst, get_hvg_residuals, SCTransform
np.random.seed(42)
import argparse


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Computes normalised counts from 10x AnnData. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )
    parser.add_argument(
        '-h5ad', '--h5ad',
        action='store',
        dest='h5ad',
        required=True,
        help=''
    )
    parser.add_argument(
        '-m', '--method',
        action='store',
        dest='method',
        required=False,
        default='cp10k',
        help=''
    )

    options = parser.parse_args()
    h5ad = options.h5ad
    adata = sc.read_h5ad(filename=h5ad)
    method = options.method
    

    available_methods = ['cp10k', 'scT']
    if method not in list(adata.layers.keys()) + available_methods:
        raise ValueError("Method not in adata.layers or available_methods.")

    if np.sum(adata.X.todense()[1,:]).is_integer():
        adata.layers['counts'] = adata.X.copy()
    elif np.sum(adata.layers['counts'].todense()[1,:]).is_integer():
        adata.X = adata.layers['counts'].copy()
    else:
        raise ValueError("Could not find raw counts in layers or X.")

    if method in list(adata.layers.keys()):
        adata.layers['dMean_normalised'] = adata.layers[method]

    if method == 'cp10k':
        # Total-count normalize (library-size correct) the data matrix X to
        # counts per million, so that counts become comparable among cells.
        adata.layers['dMean_normalised']= sc.pp.log1p(sc.pp.normalize_total(adata,
                                                        target_sum=1e4,
                                                        exclude_highly_expressed=False,
                                                        inplace=False)['X'])
    
    # if method == 'pearson':
    #     # Uses the Pearson residuals to normalise the data similar to scTransform
    #     # but working implementation in python
    #     adata.layers['dMean_normalised']= sc.pp.log1p(sc.pp.normalize_total(adata,
    #                                                     target_sum=1e4,
    #                                                     exclude_highly_expressed=False,
    #                                                     inplace=False)['X'])

    
    # if method == 'scT':
    #     # pySCTransform v2, not working with current version of packages plus takes a lot of memory
    #     vst_out = vst(
    #         adata.layers['counts'].T,
    #         gene_names=adata.var_names.tolist(),
    #         cell_names=adata.obs_names.tolist(),
    #         method="fix-slope",
    #         exclude_poisson=True,
    #         correct_counts=True)
        
    #     # Some gene are dropped by scTransform, so we need to subset the AnnData
    #     adata = adata[:,list(vst_out_v2['model_parameters_fit'].index)]
    #     adata.layers['dMean_normalised'] = vst_out["corrected_counts"].T

    print('Saving normalised AnnData...')
    adata.write('normAnnData.h5ad')


if __name__ == '__main__':
    main()