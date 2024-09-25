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
import h5py
import numpy as np
def inc_write(adata,outname):
    with h5py.File(outname, 'w') as f:
        
        # Create datasets for the essential AnnData components
        # Initialize the datasets with the appropriate shape and dtype
        X_shape = adata.X.shape
        obs_shape = adata.obs.shape
        var_shape = adata.var.shape
        
        # Assuming `adata.X` is a sparse matrix, you may want to convert it to dense or handle it differently.
        # Here we create a dataset for `X`:
        dset_X = f.create_dataset('X', shape=X_shape, dtype=adata.X.dtype)
        
        # Create datasets for `obs` and `var`
        dset_obs = f.create_dataset('obs', shape=obs_shape, dtype=h5py.string_dtype())
        dset_var = f.create_dataset('var', shape=var_shape, dtype=h5py.string_dtype())
        
        # Optionally, other components like `uns`, `obsm`, `varm` can also be created similarly.
        
        # Write the data in chunks
        chunk_size = 1000  # Define your chunk size
        
        for i in range(0, X_shape[0], chunk_size):
            end = i + chunk_size
            # Write a chunk of X
            dset_X[i:end, :] = adata.X[i:end, :].toarray() if hasattr(adata.X, "toarray") else adata.X[i:end, :]
            
            # Write a chunk of obs
            dset_obs[i:end, :] = np.array(adata.obs.iloc[i:end, :], dtype=str)
        
        # Write var data once (if it's small enough to fit in memory)
        dset_var[:, :] = np.array(adata.var, dtype=str)
<<<<<<< HEAD
=======

# def main():
#     """Run CLI."""
#     parser = argparse.ArgumentParser(
#         description="""
#             Filter and merge 10x data. Save to AnnData object.
#             """
#     )

#     parser.add_argument(
#         '-v', '--version',
#         action='version',
#         version='%(prog)s {version}'.format(version=__version__)
#     )

#     parser.add_argument(
#         '-agg', '--agg_columns',
#         action='store',
#         dest='agg_columns',
#         required=True,
#         help=''
#     )

#     parser.add_argument(
#         '-h5ad', '--h5ad',
#         action='store',
#         dest='h5ad',
#         required=True,
#         help=''
#     )

#     options = parser.parse_args()

#     h5ad = options.h5ad
#     agg_columns = options.agg_columns
#     agg_columns = agg_columns.split(",")

#     print('Reading in data...')
#     adata = sc.read_h5ad(filename=h5ad, backed='r')

#     for agg_col in agg_columns:
#         print(agg_col)
#         print("----------")
#         try:
#             data_col = adata.obs[agg_col]
#         except KeyError:
#             print(f'Aggregation column {agg_col} doesnt exist in adata')
#             continue

#         for type in data_col.unique():
#             print(type)
#             print("----------")
#             cell_adata = adata[adata.obs[agg_col] == type].copy(filename='tmp.h5ad')  # can also use .to_memory() which will be faster but more memory consuming. This is a view, not a copy
#             agg_col_cleaned = re.sub(r'\W+', '_', agg_col.replace(' ', '_'))
#             tp2 = re.sub(r'\W+', '_', type.replace(' ', '_'))

#             output_file = f'{agg_col_cleaned}__{tp2}__split.h5ad'
#             print(f'Writing to {output_file}...')

#             # Write directly from the view
#             print(f"Final shape is: {cell_adata.obs.shape}") 
#             cell_adata.write(output_file)
#             cell_adata.file.close()

#             # Remove the temporary file
#             os.remove('tmp.h5ad')
#             del cell_adata
#             gc.collect()  # Force garbage collection to free up memory

#     adata.file.close()  # Close the AnnData file to free up resources

# if __name__ == '__main__':
#     main()


#!/usr/bin/env python

# __author__ = 'Matiss Ozols and Tobi Alegbe'
# __date__ = '2021-11-25'
# __version__ = '0.0.1'

# import pandas as pd
# import scanpy as sc
# import argparse
# import os
# import re
# import gc
>>>>>>> v1.41

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
    
    parser.add_argument(
        '-conditions', '--conditions',
        action='store',
        dest='conditions',
        required=False,
        default=None,
        help=''
    )

    options = parser.parse_args()

    h5ad = options.h5ad
    agg_columns = options.agg_columns
    agg_columns = agg_columns.split(",")

    if options.conditions:
        conditions=options.conditions.split(',')
        type_analysis = 'split'
    else:
        conditions=[]
        type_analysis = 'all'
        
    print('Reading in data...')
    adata = sc.read_h5ad(filename=h5ad, backed='r')
    c=0
    for agg_col in agg_columns:
        c+=1
        
        print(agg_col)
        print("----------")
        try:
            data_col = adata.obs[agg_col]
        except KeyError:
            print(f'Aggregation column {agg_col} doesnt exist in adata')
            continue

        for type in data_col.unique():
            if type_analysis != 'all':
                if type not in conditions:
                    continue
            # print(type)
            # print("----------")
            # cell_adata = adata[adata.obs[agg_col] == type]#.copy(filename='tmp.h5ad')  # can also use .to_memory() which will be faster but more memory consuming. This is a view, not a copy
            # agg_col_cleaned = re.sub(r'\W+', '_', agg_col.replace(' ', '_'))
            # tp2 = re.sub(r'\W+', '_', type.replace(' ', '_'))

            # output_file = f'{agg_col_cleaned}__{tp2}__split.h5ad'
            # print(f'Writing to {output_file}...')

            # # Write directly from the view
            # print(f"Final shape is: {cell_adata.obs.shape}") 
            # # inc_write(cell_adata,output_file)
            # cell_adata.write(output_file)
            # cell_adata.file.close()

            # # Remove the temporary file
            # os.remove('tmp.h5ad')
            # del cell_adata
            # gc.collect()  # Force garbage collection to free up memory
            
            
            print(type)
            print("----------")
            adata.file.close() 
            del adata
            gc.collect() 
            adata = sc.read_h5ad(filename=h5ad, backed='r')
            
            cell_adata = adata[adata.obs[agg_col] == type] #.copy(filename='tmp.h5ad')  # can also use .to_memory() which will be faster but more memory consuming. This is a view, not a copy
            agg_col_cleaned = re.sub(r'\W+', '_', agg_col.replace(' ', '_'))
            tp2 = re.sub(r'\W+', '_', type.replace(' ', '_'))

            output_file = f'{agg_col_cleaned}__{tp2}__split.h5ad'
            print(f'Writing to {output_file}...')

            # Write directly from the view
            print(f"Final shape is: {cell_adata.obs.shape}") 
            cell_adata.write(output_file)
            cell_adata.file.close()
            del cell_adata

<<<<<<< HEAD
    adata.file.close()  # Close the AnnData file to free up resources
=======
>>>>>>> v1.41

if __name__ == '__main__':
    main()


#!/usr/bin/env python

# __author__ = 'Matiss Ozols and Tobi Alegbe'
# __date__ = '2021-11-25'
# __version__ = '0.0.1'

# import pandas as pd
# import scanpy as sc
# import argparse
# import os
# import re
# import gc

# def main():
#     """Run CLI."""
#     parser = argparse.ArgumentParser(
#         description="""
#             Filter and merge 10x data. Save to AnnData object.
#             """
#     )

#     parser.add_argument(
#         '-v', '--version',
#         action='version',
#         version='%(prog)s {version}'.format(version=__version__)
#     )

#     parser.add_argument(
#         '-agg', '--agg_columns',
#         action='store',
#         dest='agg_columns',
#         required=True,
#         help=''
#     )

#     parser.add_argument(
#         '-h5ad', '--h5ad',
#         action='store',
#         dest='h5ad',
#         required=True,
#         help=''
#     )

#     options = parser.parse_args()

#     h5ad = options.h5ad
#     agg_columns = options.agg_columns
#     agg_columns = agg_columns.split(",")

#     print('Reading in data...')
#     adata = sc.read_h5ad(filename=h5ad, backed='r')
    
#     for agg_col in agg_columns:
#         print(agg_col)
#         print("----------")
#         try:
#             data_col = adata.obs[agg_col]
#         except KeyError:
#             print(f'Aggregation column {agg_col} doesnt exist in adata')
#             continue

#         for type in data_col.unique():
#             print(type)
#             print("----------")
#             adata.file.close() 
#             del adata
#             gc.collect() 
#             adata = sc.read_h5ad(filename=h5ad, backed='r')
            
#             cell_adata = adata[adata.obs[agg_col] == type] #.copy(filename='tmp.h5ad')  # can also use .to_memory() which will be faster but more memory consuming. This is a view, not a copy
#             agg_col_cleaned = re.sub(r'\W+', '_', agg_col.replace(' ', '_'))
#             tp2 = re.sub(r'\W+', '_', type.replace(' ', '_'))

#             output_file = f'{agg_col_cleaned}__{tp2}__split.h5ad'
#             print(f'Writing to {output_file}...')

#             # Write directly from the view
#             print(f"Final shape is: {cell_adata.obs.shape}") 
#             cell_adata.write(output_file)
#             cell_adata.file.close()
#             del cell_adata

     # Close the AnnData file to free up resources

# if __name__ == '__main__':
<<<<<<< HEAD
#     main()
=======
#     main()
>>>>>>> v1.41
