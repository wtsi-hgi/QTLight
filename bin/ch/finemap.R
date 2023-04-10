library(data.table)
library(qqman)

gwas <- data.frame(fread("/scratch/vasccell/cs806/colocalization/cleanGWAS_Summary_Stats/GWAS_Systolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt"))
gc()

# download cytobands from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

# Extract coordinate of fes locus from cytobands
hg38cytoBand <- read.delim("cytoRegions_hg38.txt.gz", header = F)
colnames(hg38cytoBand) <- c("Chromosome", "Start", "End", "Region", "GieStain")
# remove non contigs without regions
hg38cytoBand <- hg38cytoBand[hg38cytoBand$GieStain != "", ]
# change to 1 based system
hg38cytoBand$Start <- hg38cytoBand$Start + 1
hg38cytoBand$Chromosome <- sub("chr", "", hg38cytoBand$Chromosome)
hg38cytoBand$Region <- paste0(hg38cytoBand$Chromosome, hg38cytoBand$Region)
fesCoord <- hg38cytoBand[grep("15q26.1", hg38cytoBand$Region), ]



#manhattan(gwas, snp="hg38_markername", chr="chr", bp="hg38_bp",  p="pvalue" )

# Create variables that will be used later
phreds <- 6 # no. of threads
pthresh <- 0.00000005 # pvalue threshold

auth <- paste0("Evangelou_SBP") # author of study
folder <- paste0("testFiles", "/") # dir for output files
chromosome <- as.numeric(fesCoord$Chromosome)
begin <- as.numeric(fesCoord$Start)
finish <- as.numeric(fesCoord$End)
n_samples <- 743708 # Provide number in arguments

# extract SNPs that belong to FES loci
gwasReg <- subset(gwas, gwas$chr == chromosome & gwas$hg38_bp > begin & gwas$hg38_bp < finish)

#manhattan(gwasReg, snp="hg38_markername", chr="chr", bp="hg38_bp",  p="pvalue" )

gwasReg$axe <- ifelse(gwasReg$pvalue > pthresh, gwasReg$pvalue, 0)
gwasReg <- subset(gwasReg, gwasReg$axe == 0)
gwasReg <- gwasReg[ , c("hg38_markername", "chr", "hg38_bp", "REF", "ALT", "maf", "beta", "se")]
colnames(gwasReg) <- c("hg38_markername", "chr", "hg38_bp", "allele1", "allele2", "maf", "beta", "se")

gwasRegList <- data.frame(gwasReg[, "hg38_markername"])
snpListFile <- paste0(auth, "_15q26.1_snpList.txt")
write.table(gwasRegList, file = snpListFile, quote = F, row.names = F, col.names = F)
axe1 <- system(paste0("grep -f ", snpListFile, " /scratch/vasccell/cs806/colocalization/1000Genome/1kGMerge.bim  | awk \'{print $2}\'"), intern = T)
gwasReg <- gwasReg[gwasReg$hg38_markername %in% axe1, ]
colnames(gwasReg)[1:3] <- c("rsid", "chromosome", "position")
plinkInput <- paste0(folder, auth, "_15q26.1_Region.z")
write.table(gwasReg, file = plinkInput, quote = F, col.names = T, row.names = F)

system(paste0("plink --bfile /scratch/vasccell/cs806/colocalization/1000Genome/1kGMerge --extract ", plinkInput, " --r --matrix --threads ", phreds, " --out ", plinkInput))
z <- paste0(plinkInput)
ld <- paste0(plinkInput, ".ld")
snp <- paste0(folder, auth, "_15q26.1.snp")
config <- paste0(folder, auth, "_15q26.1.config")
cred <- paste0(folder, auth, "_15q26.1.cred")
log <- paste0(folder, auth, "_15q26.1.log")
k <- paste0(folder, auth, "_15q26.1.k")
master <- data.frame(z, ld, snp, config, cred, log, k, n_samples)
masterFile <- paste0(auth, "_15q26.1_finemapMaster")
write.table(master, file = masterFile, sep = ";", quote = F, row.names = F)

# Perform Finemap and import and extract the credible set from gwas and qtl data
system(paste0("/home/c/cs806/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --sss --in-files ", masterFile, " --n-causal-snps 1 --n-iter 1000000 --n-conv-sss 500000 --n-threads ", phreds))
credSet <- read.table(paste0(folder, auth, "_15q26.1.cred1"), header = T)
credSet1 <- credSet[ , 2]
gwasReg <- subset(gwas, gwas$hg38_markername %in% credSet1)
gwasReg <- gwasReg[ , c("hg38_markername", "permID")]

credSet2 <- merge(gwasReg, credSet[, -c(1)], by.x = "hg38_markername", "cred1")
credFile <- paste0("FES_FineMap_Results/", auth, "_15q26.1_FES_Credible_Set.csv")
write.csv(credSet2, file = credFile, row.names = F)
