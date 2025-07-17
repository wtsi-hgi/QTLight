#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Transpose covariates file for jaxQTL.")
parser.add_argument("--infile", required=True, help="Path to input covariates file (rows=covariates, cols=samples)")
parser.add_argument("--outfile", required=True, help="Path to write transposed covariates file (rows=samples, cols=covariates)")
args = parser.parse_args()

df = pd.read_csv(args.infile, sep="\t", index_col=0)
df = df.transpose()
df.index.name = "sample_id"
df.to_csv(args.outfile, sep="\t")