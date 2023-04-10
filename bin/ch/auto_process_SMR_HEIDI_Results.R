# This script takes one argument
# - the path to the output file from SMR HEIDI analyses
# This script outputs  a clean smr result in csv format

# Process command line arguments
aha <- commandArgs(trailingOnly = TRUE)
aha2 <- basename(aha)
aha2 <- sub(".smr", "", aha2)

#aha2 <- paste0("zhu_SMR/smrResults/", aha)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qqman))


# Process imputeData_Nelson2017_Result ------------------------------------
smr <- read.table(aha, header = TRUE)

# Fetch rsIDs
tt <- data.frame(smr$topSNP)
snpListFile <- paste0("getRS_heidi.txt")
write.table(tt, file = snpListFile, quote = F, row.names = F, col.names = F)
system(paste0("grep -F -f ", snpListFile, " /scratch/vasccell/cs806/colocalization/dbSNP/hg38.snp151_All.bed | cut -f4,3,10 > rsLookup_heidi.txt"))
rs <- read.table("rsLookup_heidi.txt")
colnames(rs) <- c("topSNP_bp", "topSNP_rsID", "topSNP")
# Remove duplicate rsIDs
rsdedup <- rs[!(duplicated(rs$topSNP_rsID) | duplicated(rs$topSNP_rsID, fromLast = TRUE)), ]
# Add rsID to result
smr2 <- merge(smr, rsdedup[, -1], by = "topSNP", all.x = TRUE)

# Add gene names
annot <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Annot.txt"))
smr2 <- merge(smr2, annot, by.x = "Gene", by.y = "GeneID", all.x = TRUE)

smr2 <- smr2[, c(3,23,4,5,2,22,6:21)]

# Subset for Manhattan plot
smr3 <- smr2[!is.na(smr2$p_HEIDI), ]

# Plot Manhattan and export results
if (dim(smr3)[1] != 0) {
  #manhattan plot
  png(filename = paste0("/scratch/cellfunc/cs806/huvecColoc/smrResults/manhattanPlots/", aha2, "_SMRHEIDImanhattan.png"))
  
  manhattan(smr3, chr = "topSNP_chr", bp = "topSNP_bp", p = "p_HEIDI", snp = "topSNP_rsID",
            col = c("gray10", "gray60"), chrlabs = NULL, logp = TRUE,
            suggestiveline =  1.30103, genomewideline = FALSE,
            annotateTop = FALSE, annotatePval = FALSE, highlight = NULL)
  dev.off()
  
  
  smr4 <- subset(smr2, smr2$p_SMR < 0.05 & smr2$p_HEIDI > 0.05)
  smr5 <- subset(smr2, smr2$p_SMR < 0.05)
  
  write.csv(smr2, paste0("/scratch/cellfunc/cs806/huvecColoc/smrResults/resultTables/", aha2, ".csv"), row.names = F)
  write.csv(smr4, paste0("/scratch/cellfunc/cs806/huvecColoc/smrResults/resultTables/", aha2, "_sig_p_HEIDI.csv"), row.names = F)
  write.csv(smr5, paste0("/scratch/cellfunc/cs806/huvecColoc/smrResults/resultTables/", aha2, "_sig_p_SMR.csv"), row.names = F)
  
} else if (dim(smr3)[1] == 0) {
  cat(" \n")
  cat(aha2, " does not have significant SMR HEIDI tests")
  cat(" \n")
}
