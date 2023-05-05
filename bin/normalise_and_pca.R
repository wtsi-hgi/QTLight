#!/usr/bin/env Rscript
# .libPaths('/lustre/scratch118/humgen/hgi/users/mercury/scratch_mo11/r_server/r_libs_mo11')
library(edgeR)
library(DESeq2)
library(ggplot2)
# library(dplyr)
library(ggfortify)
library(ggplot2)
library(PCAtools)
# B_naive-dMean_phenotype.tsv

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
Star_path = 'M0_100_phenotype.tsv'
Mapping_Path = 'sample_mappings2.tsv'

Star_path = args[1]
Mapping_Path = args[2]
filter_type = args[3]
number_phenotype_pcs = args[4]
sc_or_bulk = args[5]
inverse_normal = as.logical(args[6])
stopifnot(inverse_normal %in% c(TRUE, FALSE))


number_phenotype_pcs = as.numeric(unlist(strsplit(number_phenotype_pcs, ',')))
max_number_phenotype_pcs =  max(number_phenotype_pcs)

# Functions taken from https://github.com/kauralasoo/eQTLUtils/blob/master/R/matrix_operations.R
quantileNormaliseVector = function(x){
  qnorm(rank(x,ties.method = "random")/(length(x)+1))
}

quantileNormaliseMatrix <- function(matrix){
  quantile_matrix = matrix(0, nrow(matrix), ncol(matrix))
  for (i in seq_along(matrix[1,])){
    quantile_matrix[,i] = quantileNormaliseVector(matrix[,i])
  }
  #Add names
  rownames(quantile_matrix) = rownames(matrix)
  colnames(quantile_matrix) = colnames(matrix)
  return(quantile_matrix)
}

quantileNormaliseRows <- function(matrix,...){
  t(quantileNormaliseMatrix(t(matrix), ...))
}

Star_counts_pre = read.table(file = Star_path, sep = '\t',check.names=FALSE, row.names = 1,header = TRUE)
colnames(Star_counts_pre)[duplicated(colnames(Star_counts_pre))]=paste0('rep_',colnames(Star_counts_pre)[duplicated(colnames(Star_counts_pre))])
percent_of_population_expressed = 0.2 #We want to only map the values that are expressed in at least 20% of donors. 
#as per https://www.medrxiv.org/content/10.1101/2021.10.09.21264604v1.full.pdf 
# We mapped cis-eQTL within a 1 megabase (MB) window of the TSS of each gene expressed
# in at least 5% of the nuclei (belonging to a broad cell type)

keep=c()
for (row in 1:nrow(Star_counts_pre)) {
  r1 = Star_counts_pre[row,]
  number_elems = length(r1)
  number_elems_greater_than_0 = length(r1[r1>0])
  if ((number_elems_greater_than_0/number_elems>percent_of_population_expressed) & (number_elems_greater_than_0>1)){
    keep <- append (keep,row)
  }
  # here also check how many elements are with equal values.
}

Star_counts_pre = Star_counts_pre[keep, ]
Star_counts_pre = t(Star_counts_pre)


Experimental_grops = read.table(Mapping_Path, fill = TRUE,check.names=FALSE,header = TRUE,sep = '\t')
Experimental_grops[duplicated(Experimental_grops[2]),2]=paste0('rep_',Experimental_grops[duplicated(Experimental_grops[2]),2])
row.names(Experimental_grops) <- Experimental_grops[,2]
Experimental_grops = Experimental_grops[,-2]
Experimental_grops2=Experimental_grops
Experimental_grops2$RNA = rownames(Experimental_grops)
Experimental_grops2 <- Experimental_grops2[, c(1,3,2)]
write.table(Experimental_grops2, file='mappings_handeling_repeats.tsv', quote=FALSE, row.names = FALSE,sep='\t')

nonzero_genes = colSums(Star_counts_pre) != 0
Star_counts_pre <- Star_counts_pre[,nonzero_genes]

Star_counts_pre_t=t(Star_counts_pre)


# keep <- apply(Star_counts_pre[0:], 1, function(x) length(unique(x[!is.na(x)])) != 1)
# f = Star_counts_pre[keep, ]


# We merge the Expreimental groups and the counts to make sure that they are in the same order, which is cruical for analysis.
Star_counts_pre <- transform(merge(Star_counts_pre, Experimental_grops, by=0,all.x = TRUE,), row.names=Row.names, Row.names=NULL)
Star_counts_pre = Star_counts_pre[complete.cases(Star_counts_pre$Sample_Category),]
#Remove the counts where there is no genotype
Star_counts_pre = Star_counts_pre[!(is.na(Star_counts_pre$Genotype) | Star_counts_pre$Genotype==""), ]


n=ncol(Experimental_grops)
Experimental_grops = (Star_counts_pre[,(ncol(Star_counts_pre)-n+1):ncol(Star_counts_pre)])

Star_counts = Star_counts_pre[,0:(ncol(Star_counts_pre)-n)]
# Star_counts = subset(Star_counts, select = -c(Row.names))
Star_counts=t(Star_counts)
all(rownames(Experimental_grops) == colnames(Star_counts))

y <- DGEList(counts=Star_counts, group=Star_counts_pre$Sample_Category,samples=Experimental_grops)
# design <- model.matrix(~ Experimental_grops$oxLDL.set + Experimental_grops$Donor.line + Experimental_grops$Sample.type..Ctrl.or.oxLDL.)

if (filter_type=='filterByExpr'){
  # this approach is not very suitable for some scRNA datasets since they are quite sarse
  keep <- filterByExpr(y,group=Star_counts_pre$Sample_Category)
  y <- y[keep, keep.lib.sizes=TRUE]
  y <- calcNormFactors(y, method = "TMM")
  
  if (lengths(unique(Experimental_grops$Sample_Category)) ==1){
    y <- estimateDisp(y)
  }else{
    design <- model.matrix(~ Experimental_grops$Sample_Category)
    y <- estimateDisp(y,design)
  }
}else if(filter_type=='HVG'){
  # Highly variable genes (HVGs) were defined as the genes in the top two quartiles based on their squared coefficient of variation (CV^2 = variance / mean^2) calculated across all cells of each different cell-type. In this manner, we identified 21,592 HVGs for the iPSC Smart-Seq2 dataset, and 16,369 for the FPP 10X dataset.
  # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02407-x#Sec12
  rows=c()
  cvs=c()
  transposed_rows=y$counts
  row_names = row.names(transposed_rows)
  keep=c()
  for (row in 1:nrow(transposed_rows)) {
    r1 = transposed_rows[row,]
    cv2 <- sd(r1)^2 / mean(r1) ^2
    rn1 = row_names[row]
    rows = append (rows, c(rn1))
    cvs = append (cvs, c(cv2))
  }
  df <- data.frame(rows, cvs)
  row.names(df)=df$rows
  cvs= df[c('cvs')]
  median_of_hvgs = median(cvs$cvs, na.rm = TRUE)
  q3 = quantile(cvs$cvs, 0.75)
  q1 = quantile(cvs$cvs, 0.25)
  keep = cvs < q3 
  y <- y[keep, keep.lib.sizes=TRUE]
  y <- calcNormFactors(y, method = "TMM")
  
}else if(filter_type=='None'){
  y=y
  y <- calcNormFactors(y, method = "TMM")
  
  
  if (lengths(unique(Experimental_grops$Sample_Category)) ==1){
    y <- estimateDisp(y)
  }else{
    design <- model.matrix(~ Experimental_grops$Sample_Category)
    y <- estimateDisp(y,design)
  }
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

# I start to doubth about this apporach - permutations sometimes fail like this.
# log2(CPM) as output
if ((sc_or_bulk == 'bulk') || grepl('-dSum', Star_path, fixed = TRUE)){
  normalised_counts <- cpm(y, log=TRUE) #https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/cpm 
  normalised_counts = normalised_counts[complete.cases(normalised_counts), ]
} else {
  normalised_counts <- y$counts
}
# Apply inverse normal transformation to each row so traits are normally distributed
if (inverse_normal == TRUE){
  print('Applying inverse normal transformation')
  normalised_counts = quantileNormaliseRows(normalised_counts)
}

# TMM_normalised_counts = t(t(y$counts)*y$samples$norm.factors)
# norms = y$samples$norm.factors
# TMM_normalised_counts_log = log(TMM_normalised_counts+1, 2) # Apply log2 transform on the TMM normalised counts.

pcs = prcomp(normalised_counts, scale = TRUE)
if(ncol(normalised_counts) < max_number_phenotype_pcs){
  max_number_phenotype_pcs=ncol(normalised_counts)
}

for (npcs in number_phenotype_pcs){
  if (max_number_phenotype_pcs < npcs){
    npcs = max_number_phenotype_pcs
  }
  pcs_sliced  = pcs$rotation[,1:npcs]
  write.table(pcs_sliced,file=paste0(npcs,'pcs.tsv'),sep='\t')
}

<<<<<<< HEAD

=======
<<<<<<< HEAD
=======

>>>>>>> 05298b280f63cbe6582fef55c60a42859d79bbba
>>>>>>> 050309403243f31ade448f98eeb701dbd2b75c95
write.table(normalised_counts,file=paste('normalised_phenotype.tsv',sep=''),sep='\t')

# plots
p <- pca(Star_counts, metadata = Experimental_grops, removeVar = 0.1)
if(ncol(normalised_counts)<15){
  len1=ncol(normalised_counts)
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