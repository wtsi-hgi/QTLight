#!/usr/bin/env python

__author__ = 'Tobi Alegbe'
__date__ = '2023-04-10'
__version__ = '0.0.1'


import scanpy as sc
from scipy.sparse import issparse
import rpy2.robjects as ro
import anndata2ri


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Run scTransform on 10x data. Save to AnnData object. Based on the script from this issue https://github.com/scverse/scanpy/issues/1068#issuecomment-590009483
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

    ro.r('library(Seurat)')
    ro.r('library(scater)')
    anndata2ri.activate()

    options = parser.parse_args()
    adata = sc.read_h5ad(filename=h5ad)

    sc.pp.filter_genes(adata, min_cells=5)
    
    if issparse(adata.X):
        if not adata.X.has_sorted_indices:
            adata.X.sort_indices()

    for key in adata.layers:
        if issparse(adata.layers[key]):
            if not adata.layers[key].has_sorted_indices:
                adata.layers[key].sort_indices()

    ro.globalenv['adata'] = adata

    ro.r('seurat_obj = as.Seurat(adata, counts="X", data = NULL)')

    ro.r('res <- SCTransform(object=seurat_obj, return.only.var.genes = FALSE, do.correct.umi = FALSE)')

    norm_x = ro.r('res@assays$SCT@scale.data').T

    adata.layers['scTransform_normalized'] = norm_x

    if output_file:
        adata.write('normAnnData.h5ad')