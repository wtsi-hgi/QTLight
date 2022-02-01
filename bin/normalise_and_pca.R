#!/usr/bin/env Rscript
# .libPaths('/lustre/scratch118/humgen/hgi/users/mercury/scratch_mo11/r_server/r_libs_mo11')

library('edgeR')
library('DESeq2')
library(ggplot2)
# library(dplyr)
library(ggfortify)
library(ggplot2)
library(PCAtools)


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

Star_path = args[1]
Mapping_Path = args[2]

Star_counts_pre = t(read.table(file = Star_path, sep = '\t',check.names=FALSE, row.names = 1,header = TRUE))
Experimental_grops = read.table(Mapping_Path, fill = TRUE,row.names = 2,check.names=FALSE,header = TRUE,sep = '\t')

nonzero_genes = colSums(Star_counts_pre) != 0
Star_counts_pre <- Star_counts_pre[,nonzero_genes]

# We merge the Expreimental groups and the counts to make sure that they are in the same order, which is cruical for analysis.
Star_counts_pre <- transform(merge(Star_counts_pre, Experimental_grops, by=0,all.x = TRUE,), row.names=Row.names, Row.names=NULL)
Star_counts_pre = Star_counts_pre[complete.cases(Star_counts_pre$Sample_Category),]


n=ncol(Experimental_grops)
Experimental_grops = (Star_counts_pre[,(ncol(Star_counts_pre)-n+1):ncol(Star_counts_pre)])
Star_counts = Star_counts_pre[,0:(ncol(Star_counts_pre)-n)]
# Star_counts = subset(Star_counts, select = -c(Row.names))
Star_counts=t(Star_counts)
all(rownames(Experimental_grops) == colnames(Star_counts))

y <- DGEList(counts=Star_counts, group=Star_counts_pre$Sample_Category,samples=Experimental_grops)
# design <- model.matrix(~ Experimental_grops$oxLDL.set + Experimental_grops$Donor.line + Experimental_grops$Sample.type..Ctrl.or.oxLDL.)


keep <- filterByExpr(y,group=Star_counts_pre$Sample_Category)
y <- y[keep, keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

if (lengths(unique(Experimental_grops$Sample_Category)) ==1){
  y <- estimateDisp(y)
}else{
  design <- model.matrix(~ Experimental_grops$Sample_Category)
  y <- estimateDisp(y,design)
}

# dge <- edgeR::DGEList(counts=dm)
# dge <- edgeR::calcNormFactors(dge)

# defaults:
# calcNormFactors(
#  dm,
#  lib.size = None, method = "TMM", refColumn = NULL,
#  logratioTrim = 0.3, sumTrim = 0.5, doWeighting = TRUE,
#  Acutoff = -1e10
#  )

# log2(CPM) as output
TMM_normalised_counts_log <- cpm(y, log=TRUE)

# TMM_normalised_counts = t(t(y$counts)*y$samples$norm.factors)
# norms = y$samples$norm.factors
# TMM_normalised_counts_log = log(TMM_normalised_counts+1, 2) # Apply log2 transform on the TMM normalised counts.

pcs = prcomp(TMM_normalised_counts_log, scale = TRUE)
if(ncol(TMM_normalised_counts_log)<20){
  len1=ncol(TMM_normalised_counts_log)
}else{
  len1=20
}


pcs20  = pcs$rotation[,1:len1]

write.table(pcs20,file=paste('pcs.tsv',sep=''),sep='\t')
write.table(TMM_normalised_counts_log,file=paste('normalised_phenotype.tsv',sep=''),sep='\t')

# plots
p <- pca(Star_counts, metadata = Experimental_grops, removeVar = 0.1)
if(ncol(TMM_normalised_counts_log)<15){
  len1=ncol(TMM_normalised_counts_log)
}else{
  len1=15
}
pdf("screenplot.pdf") 
screeplot(p, components = 1:len1)
dev.off()

pdf("biplot.pdf") 
biplot(p, showLoadings = TRUE,
       labSize = 2, pointSize = 4, sizeLoadingsNames = 2)
dev.off()

pdf("loadings.pdf") 
plotloadings(p,
             components = getComponents(p, c(1)),
             rangeRetain = 0.12, absolute = TRUE,
             col = c('black', 'pink', 'red4'),
             drawConnectors = TRUE, labSize = 4) + coord_flip()
dev.off()

# eigencorplot(p,
            #  components = getComponents(p, 1:27),
            #  metavars = c('Donor.line','Age','Sex','Ethnicity','Pool.name','oxLDL.set','Sample.type..Ctrl.or.oxLDL.'),
            #  col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
            #  cexCorval = 0.7,
            #  colCorval = 'white',
            #  fontCorval = 2,
            #  posLab = 'bottomleft',
            #  rotLabX = 45,
            #  posColKey = 'top',
            #  cexLabColKey = 1.5,
            #  scale = TRUE,
            #  main = 'PC1-27 clinical correlations',
            #  colFrame = 'white',
            #  plotRsquared = FALSE)