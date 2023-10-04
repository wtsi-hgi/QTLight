#!/usr/bin/env Rscript
library(edgeR)
library(DESeq2)
library(ggplot2)
# library(dplyr)
library(ggplot2)
library(PCAtools)
set.seed(2023)

args = commandArgs(trailingOnly=TRUE)

# This subsets the pcs to a defined dimensions.
# We do this independently of the normalisation to ensure the pipelines portability and dymamic properties.
print(args[1])
print(args[2])
pcs_file = args[1]
npcs = as.integer(args[2])
pcs = read.table(file = pcs_file, sep = '\t',check.names=FALSE, row.names = 1,header = TRUE)
pcs_sliced  = pcs[,1:npcs]
write.table(pcs_sliced,file=paste0(npcs,'pcs.tsv'),sep='\t')
