#!/usr/bin/env Rscript
# .libPaths('/lustre/scratch118/humgen/hgi/users/mercury/scratch_mo11/r_server/r_libs_mo11')
# library('tximport')
library('edgeR')
library('DESeq2')
library(ggplot2)
library(dplyr)
library(devtools)
# install_github('sinhrks/ggfortify')
library(ggfortify)
library(ggplot2)

limix_analysis_folder = '/lustre/scratch123/hgi/projects/macromap/Analysis/LIMIX/Sanger_IDs_2/interaction_vst'
# t= read.table(file = '/lustre/scratch123/hgi/projects/macromap/rna_seq_5044_nextflow/results/tximport/txi_lengthScaledTPM_gene_counts.csv', sep = ',', row.names = 1,header = TRUE)
# read the Star counts
Star_counts_load = t(read.table(file = '/lustre/scratch123/hgi/projects/macromap/rna_seq_5044_nextflow/results/combined/star-fc-genecounts.txt', sep = '\t', row.names = 1,header = TRUE))
# Experimental_grops = read.table('/lustre/scratch123/hgi/projects/macromap/Analysis/LIMIX/Maria_oxLDL project groups2.csv',row.names = 1,header = TRUE,sep = ',')
Experimental_grops_load = read.table('/lustre/scratch123/hgi/projects/macromap/Analysis/QC/Maria_oxLDL project metadata2.csv',, fill = TRUE,row.names = 1,header = TRUE,sep = ',')
Experimental_grops = Experimental_grops_load
Star_counts_pre = Star_counts_load
# We merge the Expreimental groups and the counts to make sure that they are in the same order, which is cruical for analysis.
Star_counts_pre <- merge(Star_counts_pre, Experimental_grops, by=0)
# Experimental_grops$Label = rownames(Experimental_grops)
# some of the samples may not be annotated or may be part of other grop, hence here we exclude these.
Star_counts_pre = Star_counts_pre[complete.cases(Star_counts_pre$Sample.type..Ctrl.or.oxLDL.),]
row.names(Star_counts_pre) <- Star_counts_pre$Row.names
n=ncol(Experimental_grops)
Experimental_grops = (Star_counts_pre[,(ncol(Star_counts_pre)-n+1):ncol(Star_counts_pre)])
Star_counts = Star_counts_pre[,0:(ncol(Star_counts_pre)-n)]
Star_counts = subset(Star_counts, select = -c(Row.names))
# drops <- c("Sample.type..Ctrl.or.oxLDL.","Sample_Category",'Row.names','Donor.line')
# Star_counts = t(Star_counts_pre[ , !(names(Star_counts_pre) %in% drops)])
Star_counts=t(Star_counts)

# Filter genes with low counts
Star_counts = Star_counts[rowSums(Star_counts > 1) >= 10,]
removed_gene_nr = ncol(Star_counts_pre)-nrow(Star_counts)
paste('Removed',removed_gene_nr,'genes')
# For DE Might want to filter down to only Coding genes

# Check if all the data is in the correct order
if (all(rownames(Experimental_grops) == colnames(Star_counts)))
  print('The order is correct')

# To double check select also 
# heatmap(head(Star_counts))

# ##########################
#edgeR normalisation methods
# ##########################

# TMM normalise the dataset - this is using edgeR
y <- DGEList(counts=Star_counts, group=Star_counts_pre$Sample.type..Ctrl.or.oxLDL.,samples=Experimental_grops)
design <- model.matrix(~Experimental_grops$Sample.type..Ctrl.or.oxLDL.)
design<-formula(~Experimental_grops$Sample.type..Ctrl.or.oxLDL.)
# design <- model.matrix(~ Experimental_grops$oxLDL.set + Experimental_grops$Donor.line + Experimental_grops$Sample.type..Ctrl.or.oxLDL.)

keep <- filterByExpr(y,group=Star_counts_pre$Sample.type..Ctrl.or.oxLDL.,design=design)
y <- y[keep, keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
TMM_normalised_counts = cpm(y, normalized.lib.size=T,log=T, prior.count=1)
normalised_counts=TMM_normalised_counts
# ##########################
# Deseq2 normalisation
# ##########################

# https://support.bioconductor.org/p/92455/
# dds <- DESeqDataSetFromMatrix(Star_counts, Experimental_grops, design= ~ oxLDL.set + Donor.line + Sample.type..Ctrl.or.oxLDL.)
# Star_counts2 <- Star_counts[rowSums(Star_counts) >= 50,]
# dds <- DESeqDataSetFromMatrix(Star_counts, Experimental_grops, design= ~ Donor.line +Sample.type..Ctrl.or.oxLDL.)
design<-formula(~Sample.type..Ctrl.or.oxLDL.)
dds <-DESeqDataSetFromMatrix(Star_counts, Experimental_grops, design=design, ignoreRank = FALSE)
# dds <- DESeq(dds,parallel=TRUE,BPPARAM=12)
dds <- DESeq(dds,betaPrior=FALSE)
resultsNames(dds) 
vst_normalised_counts <- assay(vst(dds))
normCounts<-as.data.frame(counts(dds,normalized=TRUE))
colnames(normCounts) <- paste0(colnames(normCounts),'_normalised')
write.csv(normCounts,"DESeq_outputNormalised.csv")

save(dds, Experimental_grops, design,Star_counts, file = "Important_stuff.RData")
load("Important_stuff.RData")

normalised_counts=vst_normalised_counts
normalised_counts=Star_counts
normalised_counts= subset(normalised_counts, select = -c(MM_oxLDL8032568))
Experimental_grops= t(subset(t(Experimental_grops), select = -c(MM_oxLDL8032568)))
Experimental_grops =data.frame(Experimental_grops)
library(PCAtools)
p <- pca(normalised_counts, metadata = Experimental_grops, removeVar = 0.1)
screeplot(p, components = 1:15)

biplot(p, x='PC1',y='PC2',showLoadings = FALSE,
        labSize = 3, pointSize = 0.4, colby = 'Sample.type..Ctrl.or.oxLDL.', sizeLoadingsNames = 5,ellipse=F)

pairsplot(p,
          components = getComponents(p, c(1:10)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.4,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Sample.type..Ctrl.or.oxLDL.',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

eigencorplot(p,
             components = getComponents(p, 1:27),
             metavars = c('Donor.line','Age','Sex','Ethnicity','Pool.name','oxLDL.set','Sample.type..Ctrl.or.oxLDL.'),
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             cexCorval = 1.2,
             fontCorval = 2,
             posLab = 'all',
             rotLabX = 45,
             scale = TRUE,
             main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))

eigencorplot(p,
             components = getComponents(p, 1:27),
             metavars = c('Donor.line','Age','Sex','Ethnicity','Pool.name','oxLDL.set','Sample.type..Ctrl.or.oxLDL.'),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'PC1-27 clinical correlations',
             colFrame = 'white',
             plotRsquared = FALSE)

plotloadings(p,
             components = getComponents(p, c(1)),
             rangeRetain = 0.12, absolute = TRUE,
             col = c('black', 'pink', 'red4'),
             drawConnectors = TRUE, labSize = 4) + coord_flip()


# Loop through each of the groups split the PCs and Phenotypes accordingly and write to the file
for (group_category in unique(unlist(Star_counts_pre$Sample_Category, use.names = FALSE))) {
  # split the PCs and Phenotypes accordingly and write to the file
  
  print(group_category)
  idx <- Experimental_grops$Sample_Category == group_category
  experiment_ids = Experimental_grops[idx, 'Label']
  
  Group_countd = normalised_counts[,match(experiment_ids, colnames(normalised_counts))]
  pcs = pcs20[match(experiment_ids, rownames(pcs20)),]
  
  # replace the Sanger ids with the individual ids
  # if any of the individuals are repeated then take the first one (this may need to be changed).
  
  # for the expression data
  idl = data.frame(match(colnames(Group_countd),Experimental_grops$Label))
  idl$Cols = colnames(Group_countd)
  idl$Donor = Experimental_grops[idl$match.colnames.Group_countd...Experimental_grops.Label.,'Donor.line']
  idl = distinct(idl,Donor, .keep_all= TRUE)
  Group_countd = Group_countd[,idl$Cols]
  # colnames(Group_countd) = idl$Donor
  
  # for the PCs
  idl = data.frame(match(rownames(pcs),Experimental_grops$Label))
  idl$Cols = rownames(pcs)
  idl$Donor = Experimental_grops[idl$match.rownames.pcs...Experimental_grops.Label.,'Donor.line']
  idl = distinct(idl,Donor, .keep_all= TRUE)
  pcs = pcs[idl$Cols,]
  # rownames(pcs) = idl$Donor
  # link = paste(limix_analysis_folder,group_category,sep='/')
  # dir.create(link)
  # write.table(pcs,file=paste(link,'/',group_category,'_pcs.tsv',sep=''),sep='\t')
  # write.table(Group_countd,file=paste(link,'/',group_category,'_phenotype.tsv',sep=''),sep='\t')
}

# Plot the PCs - I used this to see if repeated individuals cluster together, but it did nont seem so.

library(PCAtools)
p <- pca(vst, metadata = Experimental_grops, removeVar = 0.1)
screeplot(p, components = 1:15)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

# horn method to estimate the optimal number of PCS.
horn <-parallelPCA(vst)
horn$n
pairsplot(p,
          components = getComponents(p, c(1:10)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.4,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Donor.line',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

biplot(p,
       colby = 'Sample.type..Ctrl.or.oxLDL.', colkey = c('M0_oxLDL' = 'forestgreen', 'M0_Ctrl' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 0,
       ellipseFillKey = c('M0_oxLDL' = 'yellow', 'M0_Ctrl' = 'pink'),
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)

plotloadings(p,
             components = getComponents(p, c(26)),
             rangeRetain = 0.12, absolute = TRUE,
             col = c('black', 'pink', 'red4'),
             drawConnectors = TRUE, labSize = 4) + coord_flip()

plotloadings(p,
             components = getComponents(p, c(2)),
             rangeRetain = 0.12, absolute = TRUE,
             col = c('black', 'pink', 'red4'),
             drawConnectors = TRUE, labSize = 4) + coord_flip()


eigencorplot(p,
             components = getComponents(p, 1:27),
             metavars = c('Donor.line','Age','Sex','Ethnicity','Pool.name','oxLDL.set'),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'PC1-27 clinical correlations',
             colFrame = 'white',
             plotRsquared = FALSE)

eigencorplot(p,
             components = getComponents(p, 1:27),
             metavars = c('Donor.line','Age','Sex','Ethnicity','Pool.name','oxLDL.set'),
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             cexCorval = 1.2,
             fontCorval = 2,
             posLab = 'all',
             rotLabX = 45,
             scale = TRUE,
             main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))

res <- results(dds, name="Sample.type..Ctrl.or.oxLDL._M0_oxLDL_vs_M0_Ctrl")
resLFC <- lfcShrink(dds, coef="Sample.type..Ctrl.or.oxLDL._M0_oxLDL_vs_M0_Ctrl", type="apeglm")
plotMA(resAsh, ylim=c(-2,2))

idx <- identify(resLFC$baseMean, resLFC$log2FoldChange)
rownames(resLFC)[idx]

resultsNames(dds)
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

plotCounts(dds, gene=4, intgroup='Sample.type..Ctrl.or.oxLDL.')
mcols(res)$description


ddsMF <- dds


library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sample.type..Ctrl.or.oxLDL.
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

