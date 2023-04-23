#!/usr/bin/env python

__author__ = 'Tobi Alegbe'
__date__ = '2023-04-10'
__version__ = '0.0.1'


import scanpy as sc
import pandas as pd
import numpy as np
from pysctransform import vst, get_hvg_residuals, SCTransform
np.random.seed(42)
import argparse


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Run SCTransform on 10x AnnData. Save to AnnData object. Python port found here https://github.com/saketkc/pySCTransform/tree/develop
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

    options = parser.parse_args()
    h5ad = options.h5ad
    adata = sc.read_h5ad(filename=h5ad)
    

    # scTransform v1, works but needs a lot of memory
    vst_out = vst(
        adata.layers['counts'].T,
        gene_names=adata.var_names.tolist(),
        cell_names=adata.obs_names.tolist(),
        method="theta_ml",
        n_cells=2000,
        verbosity=1)

    # scTransform v2, not working
    # vst_out = vst(
    #     adata.layers['counts'].T,
    #     gene_names=adata.var_names.tolist(),
    #     cell_names=adata.obs_names.tolist(),
    #     method="fix-slope",
    #     exclude_poisson=True,
    #     correct_counts=True)
    
    norm_x = pd.DataFrame(data=vst_out["corrected_counts"].T,
                          index = vst_out['model_parameters_fit'].index)

    adata.layers['scTransform_normalised'] = norm_x

    adata.write('normAnnData.h5ad')


if __name__ == '__main__':
    main()