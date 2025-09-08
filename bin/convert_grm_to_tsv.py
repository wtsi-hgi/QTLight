#!/usr/bin/env python

#take a GRM in mtx format and save as tsv with sample_id as a header for use in Quasar
import argparse
import pandas as pd
from scipy.io import mmread

def get_options():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="""
            Convert GRM in mtx format to tsv format with sample IDs as header.
            """
    )

    parser.add_argument(
        '-g', '--grm',
        action='store',
        dest='grm',
        required=True,
        help='Path to GRM in mtx format.'
    )

    parser.add_argument(
        '-s', '--samples',
        action='store',
        dest='samples',
        required=True,
        help='Path to sample IDs file.'
    )

    parser.add_argument(
        '-o', '--output',
        action='store',
        dest='output',
        required=True,
        help='Path to output tsv file.'
    )

    args = parser.parse_args()
    return args


def main():
    args = get_options()
    # Load GRM matrix
    matrix = mmread(args.grm)
    outfile = args.output
    dense_matrix = matrix.toarray()
    df = pd.DataFrame(dense_matrix)
    sample_ids = pd.read_csv(args.samples, header=None)[0].tolist()
    df = pd.DataFrame(dense_matrix, index=sample_ids, columns=sample_ids)
    df.to_csv(outfile, sep="\t", index=True, header=True, index_label="sample_id")


if __name__ == "__main__":
    main()