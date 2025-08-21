
library(data.table)
library(doMC)
library(tidyverse)

folder <- paste0("TenQTL_HUVEC_Finemap2", "/") # dir for output files
phreds <- 6 # no. of threads
n_samples <- 96 # Provide number in arguments

qtl <- data.frame(fread("tenQTLs/HUVEC_RNA-Seq_10pcs_Cis_eqtls_pval005.tsv.gz"))
gc()

# Subset genes that have genome-wide significant snp associations
sigQTL <- qtl[which(qtl$pvalue < 5e-08), ]
sigQTL <- sigQTL[order(sigQTL$pvalue), ]

# How many genes have genome wide significant associations
topSigQTL <- sigQTL[!duplicated(sigQTL$gene), ]

# Prepare genes with genomewide significant QTL for eCAVIAR 
foreach(rr = 1:nrow(topSigQTL), .combine = rbind) %dopar% {
  eGene <- topSigQTL[rr, "gene"]
  topSNP <- topSigQTL[rr, "snps"]
  topSNP_bp <- as.numeric(topSigQTL[rr, "position"])
  topSNP_chr <- topSigQTL[rr, "chr"]
  
  # Get Summary Stats for gene
  eGene_SS <- subset(qtl, qtl$gene == eGene)
  eGene_SS$snps <- paste0("chr", eGene_SS$snps)
  
  # Prep for Finemap
  qtlSNPList <- data.frame(eGene_SS[, "snps"])
  qtlSNPListFile <- paste0("qtlSNPList.txt")
  write.table(qtlSNPList, file = qtlSNPListFile, quote = F, row.names = F, col.names = F)
  axe1 <- system(paste0("grep -f ", qtlSNPListFile, " /scratch/cellfunc/shared/HUVEC_genotype/Imputed_Genotypes_Data_All_HUVEC_Samples.bim  | awk \'{print $2}\'"), intern = T)
  eGene_SS <- eGene_SS[eGene_SS$snps %in% axe1, ]
  eGene_SS2 <- eGene_SS[, c("snps", "chr", "position", "allele1", "allele2", "maf", "beta", "se")]
  colnames(eGene_SS2)[1:5] <- c("rsid", "chromosome", "position", "allele1", "allele2")
  plinkInput <- paste0(folder, eGene, "_snpList.z")
  write.table(eGene_SS2, file = plinkInput, quote = F, col.names = T, row.names = F)
  
  system(paste0("plink --bfile /scratch/cellfunc/shared/HUVEC_genotype/Imputed_Genotypes_Data_All_HUVEC_Samples --allow-extra-chr --extract ", plinkInput, " --r --matrix --threads ", phreds, " --out ", plinkInput))
  z <- paste0(plinkInput)
  ld <- paste0(plinkInput, ".ld")
  snp <- paste0(folder, eGene, ".snp")
  config <- paste0(folder, eGene, ".config")
  cred <- paste0(folder, eGene, ".cred")
  log <- paste0(folder, eGene, ".log")
  k <- paste0(folder, eGene, ".k")
  master <- data.frame(z, ld, snp, config, cred, log, k, n_samples)
  masterFile <- paste0(eGene, "_finemapMaster")
  write.table(master, file = masterFile, sep = ";", quote = F, row.names = F)
  
  # Perform Finemap and import and extract the credible set from gwas and qtl data
  system(paste0("/home/c/cs806/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --sss --in-files ", masterFile, " --n-causal-snps 1 --n-iter 1000000 --n-conv-sss 500000 --n-threads ", phreds))
  system(paste0("rm ", eGene, "_finemapMaster"))
  
  credSet <- read.table(paste0(folder, eGene, ".cred1"), header = T)
  eGene_credSet <- merge(eGene_SS, credSet[, -1], by.x = "snps", by.y = "cred1")
  eGene_credSet$snps <- sub("chr", "", eGene_credSet$snps)
  
  credFile <- paste0(folder, "credset/", eGene, "_Credible_Set.txt")
  write.table(eGene_credSet, file = credFile, row.names = F)
  system(paste0("rm ", folder, eGene, "*"))
  
  
}

# combine files and export

credDF <- list.files(path=paste0(folder, "credset"), full.names = TRUE) %>% 
  lapply(read.table, stringsAsFactors = FALSE, header = TRUE) 
credDF2 <- do.call("rbind", credDF)

fwrite(credDF2, file = paste0(folder, "credset/", "tensorQTL_HUVEC_FineMap_Credible_Set.csv"))
