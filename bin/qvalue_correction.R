#!/usr/bin/env Rscript

# Bradley Dec 2023

library(dplyr)
library(qvalue)
# Function to parse command-line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  args_list <- list()
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("-f", "--file")) {
      args_list$file <- args[i + 1]
      i <- i + 1
    } else if (args[i] %in% c("-c", "--column")) {
      args_list$column <- as.numeric(args[i + 1])
      i <- i + 1
    } else if (args[i] %in% c("-n", "--new_q_val_column")) {
      args_list$new_q_val_column <- args[i + 1]
      i <- i + 1
    } else if (args[i] %in% c("-w", "--within_gene")) {
      args_list$within_gene <- args[i + 1]
      i <- i + 1
    }
    i <- i + 1
  }
  return(args_list)
}

# Parse the arguments
args <- parse_args()

# Check if all required arguments are provided
if (is.null(args$file) || is.null(args$column) || is.null(args$new_q_val_column) || is.null(args$within_gene)) {
  stop("All arguments -f, -c, -n, and -w are required")
}

# Assign arguments to variables
# file='/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/saige/v1/work/66/6e28d628324ab5ddb2378689d8b838/output_Azimuth_predicted_celltype_l1__CD8_T__split_h5ad___1/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_cis_ENSG00000117862'
# column='13'
# new_column='qval'
# within=TRUE
file <- args$file
column <- args$column
new_column <- args$new_q_val_column
within <- args$within_gene

# Print the parsed arguments (optional, for debugging)
cat("File:", file, "\n")
cat("Column:", column, "\n")
cat("New Q-value column:", new_column, "\n")
cat("Within gene:", within, "\n")

# Load data
res = read.delim(file, header=F)
if(is.na(as.numeric(res[,as.numeric(column)][1]))){
    colnames(res) = res[1,]
    res = res[-1,]
}

# Extract p_values
p_values = as.numeric(res[,as.numeric(column)])

# Generate q-value object
qobj1 = qvalue(p = p_values); 
qobj = qobj1$qvalues

# qobj =p.adjust(p_values, method = c("BH"),n = length(p_values))
# qobj <- p.adjust(p_values, method = "BY", n = length(p_values))
# Add q-values
res[,new_column] = qobj

# Add local fdr
# res$lfdr = qobj$lfdr

# Replace the original file with the current one
write.table(res, file, col.names=T, quote=F, sep = "\t",row.names=F)

# Save a file with the variant with the minimum qvalue for this gene
min_file = gsub(".txt", "", file)
min_file = paste0(min_file, "_minimum_q.txt")
to_save = res[res[,new_column] == min(res[,new_column]),]
# If little significance, may get lots of results with the minimum qvalue
if(nrow(to_save > 1)){
    # If have lots of variants with same corrected value, take the one with the smallest original pvalue
    to_save=res[res[,column] == min(res[,column]),]
    # If have variants in complete LD - take the top
    if(nrow(to_save) > 1){
        to_save = to_save[1,]
    }
}

if(within == "TRUE"){
    write.table(to_save, min_file, col.names=T, quote=F, sep = "\t", row.names=F)
}
print("Performed q-value correction")