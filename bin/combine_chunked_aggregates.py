#!/usr/bin/env python

import pandas as pd, glob
import argparse

"""Run CLI."""
parser = argparse.ArgumentParser(
    description="""
        Filter and merge 10x data. Save to AnnData object.
        """
)


parser.add_argument(
    '-n_individ', '--n_min_individ',
    action='store',
    dest='n_individ',
    required=True,
    help=''
)
options = parser.parse_args()
# n_min_individ
n_individ = int(options.n_individ)

files = sorted(glob.glob("*___phenotype_file.tsv"))
assert files, "No phenotype files matched"

dfs = []
for f in files:
    df = pd.read_csv(f, sep="\t")
    if "ensembl_id" not in df.columns:
        df.rename(columns={df.columns[0]: "ensembl_id"}, inplace=True)
    df = df.set_index("ensembl_id")
    # drop any duplicate columns within a chunk (just in case)
    df = df.loc[:, ~df.columns.duplicated()]
    dfs.append(df)

# inner join keeps genes common to all chunks; use 'outer' if you prefer a union
merged = pd.concat(dfs, axis=1, join="inner")

# if any duplicate column names occur across chunks, drop true duplicates
dup_mask = merged.columns.duplicated()
if dup_mask.any():
    merged = merged.loc[:, ~dup_mask]
files[0].split('--')[-1]

if n_individ < 1:
    raise ValueError(f"--n_min_individ must be >= 1, got {n_individ}")

if merged.shape[1] >= n_individ:
    merged.reset_index().to_csv(files[0].split('--')[-1],
                                sep="\t", index=False)
    print("Merged phenotype shape:", merged.shape)


    # ---------- MAPPING MERGE ----------
    map_files = sorted(glob.glob("*___genotype_phenotype_mapping.tsv"))
    assert map_files, "No genotype_phenotype_mapping files matched"

    maps = []
    for f in map_files:
        df = pd.read_csv(f, sep="\t", dtype=str)
        # normalize common header variants
        df.columns = [c.strip() for c in df.columns]
        expected = {"Genotype", "RNA", "Sample_Category"}
        if not expected.issubset(set(df.columns)):
            raise ValueError(f"{f}: expected columns {expected}, got {set(df.columns)}")
        maps.append(df)

    merged_map = pd.concat(maps, axis=0, ignore_index=True).drop_duplicates()

    merged_map.to_csv(map_files[0].split('--')[-1], sep="\t", index=False)
    print("Merged mapping rows:", merged_map.shape[0])