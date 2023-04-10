
# Prepare ATAC-SEQ QTL Summ Stats
# Charles Solomon
# 27/03/2022

# libraries
library(data.table)

# Assign file names
aha1 <- "tenQTLs/HUVEC_ChipSeq_Cis_eqtlsAllNominal.tsv.gz"

# Import data
genePos <- data.frame(fread("/scratch/cellfunc/shared/HUVEC_ChipSeq/eQTLs/prepeared_counts_and_annotation/chip_gtf_input.tsv"))
colnames(genePos) <- c("gene_id", "Chr", "Start", "Stop", "Strand", "Length" )

frq <- data.frame(fread("../huvecEQTL/huvecImputeGenoEQTL/imputeMAF.frq"))

geQTL <- data.frame(fread(paste0(aha1)))
geQTL <- geQTL[geQTL$pval_nominal < 0.05, ]
gc()

geQTL <- merge(geQTL, frq, by.x = "variant_id", by.y = "SNP")

geQTL <- setnames(geQTL, c("phenotype_id", "variant_id", "slope", "slope_se",
                           "pval_nominal", "tss_distance", "CHR", "A1", "A2", "MAF"),
                  c("gene", "snps", "beta", "se", "pvalue", "dist",
                    "chr", "allele1", "allele2", "maf"),
                  skip_absent=TRUE)


#geQTL <- geQTL[, c("gene_id", "variant", "rsid", "chromosome", "position", "ref",
#                   "alt", "maf", "beta", "se", "pvalue", "type", "ac",
#                   "ma_samples", "an")]
gc()

# Add essential columns 
geQTL$snps <- gsub("chr", "", geQTL$snps)
geQTL$position <- sapply(geQTL$snps, function(x){unlist(strsplit(x, split = ":"))[2]})
geQTL$statistic <- geQTL$beta / geQTL$se

gc()

# Prep export variables

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))
aha2 <- sub("AllNominal.tsv.gz", "", aha1)
aha2 <- split_path(aha2)
aha2 <- paste0("tenQTLs/", aha2[1])


aha3 <- paste0(aha2, "_pval005.tsv.gz")
aha4 <- paste0(aha2, "_pval005_SMR.txt")

fwrite(geQTL, paste0(aha3))

# Prep for SMR
geQTL2 <- geQTL[!duplicated(geQTL[c("gene", "snps")]), ]
geQTL3 <- geQTL2[, c("gene", "snps", "dist", "pvalue", "beta")]

fwrite(geQTL3, paste0(aha4), sep = "\t", col.names = F)


system(paste0("/home/c/cs806/SMR/smr_v1.3.1_linux_x86_64_static --eqtl-summary ", aha4,  " --fastqtl-nominal-format --make-besd --out ", aha2))


# Restart R to clear memory --
esi <- fread(paste0(aha2, ".esi"))
esi$id  <- 1:nrow(esi)
geQTL4 <- geQTL2[!duplicated(geQTL2[c("snps")]), ]
esi <- merge(esi, geQTL4[, c("chr", "snps", "allele1", "allele2", "maf", "position")], by.x = "V2", by.y = "snps")
esi <- esi[order(esi$id), ]
esi <- esi[, c("chr", "V2", "V3", "position", "allele1", "allele2", "maf")]
fwrite(esi, paste0(aha2, ".esi"), sep = "\t", col.names = FALSE)

epi <- fread(paste0(aha2, ".epi"))
epi$id  <- 1:nrow(epi)
epi <- merge(epi, genePos, by.x = "V2", by.y = "gene_id")
epi <- epi[order(epi$id), ]
epi <- epi[, c("Chr", "V2", "V3", "Start", "V2", "Strand")]
fwrite(epi, paste0(aha2, ".epi"), sep = "\t", col.names = FALSE)

