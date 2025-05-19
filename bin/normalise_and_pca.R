#!/usr/bin/env Rscript

library(edgeR)
library(DESeq2)
library(ggplot2)
library(PCAtools)
set.seed(2023)
library(ggfortify)
library(dplyr)
library(ggrepel)
library(reshape2)
library(rlang)
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# dMean__plasma_cells___phenotype_file.tsv remap_dMean__plasma_cells___genotype_phenotype_mapping.tsv None single_cell TRUE NONE 0.1 true
Star_path = 'splicing_phenotype_input.tsv'
Mapping_Path = 'fake_file.fq'
filter_type = 'None'

# number_phenotype_pcs = args[4]
sc_or_bulk = 'single_cell'
inverse_normal = as.logical('FALSE')
stopifnot(inverse_normal %in% c(TRUE, FALSE))
norm_method = 'NONE'
percent_of_population_expressed = '0.2'
pc_strat='true'


Star_path = args[1]
Mapping_Path = args[2]
filter_type = args[3]
# number_phenotype_pcs = args[4]
sc_or_bulk = args[4]
inverse_normal = as.logical(args[5])
stopifnot(inverse_normal %in% c(TRUE, FALSE))


norm_method = args[6]
percent_of_population_expressed = as.numeric(args[7])
if (args[8]=='false'){
  pc_strat = 'FALSE'
}

if (args[8]=='true'){
  pc_strat = 'TRUE'
}

use_sample_pca = as.logical(pc_strat)

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

original_cols <- colnames(Star_counts_pre)
dup_tracker <- table(original_cols)
new_cols <- character(length(original_cols))
name_counter <- list()

for (i in seq_along(original_cols)) {
  name <- original_cols[i]
  if (dup_tracker[name] > 1) {
    count <- name_counter[[name]] %||% 1  # use %||% from rlang or just `if (is.null())`
    new_cols[i] <- paste0("rep", count, "_", name)
    name_counter[[name]] <- count + 1
  } else {
    new_cols[i] <- name
  }
}

colnames(Star_counts_pre) <- new_cols
# percent_of_population_expressed = 0.2 #We want to only map the values that are expressed in at least 20% of donors. 
#as per https://www.medrxiv.org/content/10.1101/2021.10.09.21264604v1.full.pdf 
# We mapped cis-eQTL within a 1 megabase (MB) window of the TSS of each gene expressed
# in at least 5% of the nuclei (belonging to a broad cell type)

if (percent_of_population_expressed>0){
  keep=c()
  for (row in 1:nrow(Star_counts_pre)) {
    r1 = Star_counts_pre[row,]
    number_elems = length(r1)
    number_elems_greater_than_0 = length(r1[r1>0])
    if ((number_elems_greater_than_0/number_elems>=percent_of_population_expressed) & (number_elems_greater_than_0>1)){
      keep <- append (keep,row)
    }
    # here also check how many elements are with equal values.
  }
  Star_counts_pre = Star_counts_pre[keep, ]
  
}
if ("gene_name" %in% colnames(Star_counts_pre)) {
  Star_counts_pre$gene_name <- NULL
}
Star_counts_pre = t(Star_counts_pre)


if (grepl("fake_file", Mapping_Path)) {
  # Assume Star_counts_pre is already loaded in the environment
  row_data <- rownames(Star_counts_pre)
  Sample_Category <- tools::file_path_sans_ext(basename(Star_path))
  Experimental_grops2 <- data.frame(
    Genotype = row_data,
    RNA = row_data,
    Sample_Category = Sample_Category
  )
  Experimental_grops2 <- Experimental_grops2[Experimental_grops2$Genotype != "gene_name", ]
  # For fake mapping path:
  Experimental_grops2$Genotype <- sub("^rep\\d+_", "", Experimental_grops2$Genotype)

  rownames(Experimental_grops2) = Experimental_grops2$RNA
} else {
  Experimental_grops = read.table(Mapping_Path, fill = TRUE,check.names=FALSE,header = TRUE,sep = '\t')
  original_ids <- Experimental_grops[[2]]
  id_counts <- table(original_ids)
  new_ids <- character(length(original_ids))
  id_tracker <- list()

  for (i in seq_along(original_ids)) {
    id <- original_ids[i]
    if (id_counts[id] > 1) {
      count <- if (is.null(id_tracker[[id]])) 1 else id_tracker[[id]]
      new_ids[i] <- paste0("rep", count, "_", id)
      id_tracker[[id]] <- count + 1
    } else {
      new_ids[i] <- id
    }
  }

  # Update column 2 with new unique IDs
  Experimental_grops[[2]] <- new_ids
  row.names(Experimental_grops) <- Experimental_grops[,2]
  Experimental_grops = Experimental_grops[,-2]
  Experimental_grops2=Experimental_grops
  Experimental_grops2$RNA = rownames(Experimental_grops)
  Experimental_grops2 <- Experimental_grops2[, c(1,3,2)]
}

write.table(Experimental_grops2, file='mappings_handeling_repeats.tsv', quote=FALSE, row.names = FALSE,sep='\t')
Experimental_grops <- Experimental_grops2

nonzero_genes = colSums(Star_counts_pre) != 0
Star_counts_pre <- Star_counts_pre[,nonzero_genes]

Star_counts_pre_t=t(Star_counts_pre)


# keep <- apply(Star_counts_pre[0:], 1, function(x) length(unique(x[!is.na(x)])) != 1)
# f = Star_counts_pre[keep, ]


# We merge the Expreimental groups and the counts to make sure that they are in the same order, which is cruical for analysis.
Star_counts_pre = merge(Star_counts_pre, Experimental_grops, by=0,all.x = TRUE,)
rownames(Star_counts_pre) <- Star_counts_pre[,'Row.names']
# Star_counts_pre = Star_counts_pre[,!(names(Star_counts_pre) %in% c("Row.names"))]


Star_counts_pre = Star_counts_pre[complete.cases(Star_counts_pre$Sample_Category),]
#Remove the counts where there is no genotype
Star_counts_pre = Star_counts_pre[!(is.na(Star_counts_pre$Genotype) | Star_counts_pre$Genotype==""), ]

n=ncol(Experimental_grops)
Experimental_grops = (Star_counts_pre[,(ncol(Star_counts_pre)-n+1):ncol(Star_counts_pre)])

Star_counts = Star_counts_pre[, seq_len(ncol(Star_counts_pre) - n)]

Star_counts$Row.names <- NULL
# Star_counts = subset(Star_counts, select = -c(Row.names))
Star_counts=t(Star_counts)
all(rownames(Experimental_grops) == colnames(Star_counts))

if (norm_method=='NONE' && filter_type=='None'){
    normalised_counts=Star_counts
}else{
    y <- DGEList(counts=Star_counts, group=Star_counts_pre$Sample_Category,samples=Experimental_grops)

    # design <- model.matrix(~ Experimental_grops$oxLDL.set + Experimental_grops$Donor.line + Experimental_grops$Sample.type..Ctrl.or.oxLDL.)
    if (filter_type=='filterByExpr'){
      # this approach is not very suitable for some scRNA datasets since they are quite sarse
      keep <- filterByExpr(y,group=Star_counts_pre$Sample_Category)
      y <- y[keep, keep.lib.sizes=TRUE]
      y <- calcNormFactors(y, method = "TMM")
      
      if (length(unique(Experimental_grops$Sample_Category)) == 1){
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
      y <- calcNormFactors(y, method = "TMM")
      if (length(unique(Experimental_grops$Sample_Category)) ==1){
        y <- estimateDisp(y)
      }else{
        design <- model.matrix(~ Experimental_grops$Sample_Category)
        y <- estimateDisp(y,design)
      }
    }

    if ((sc_or_bulk == 'bulk') || grepl('-dSum', Star_path, fixed = TRUE)){
      if (norm_method=='TMM'){
        normalised_counts <- cpm(y, log=TRUE) #https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/cpm 
        normalised_counts = normalised_counts[complete.cases(normalised_counts), ]
      }else if(norm_method=='DESEQ'){
        counts = y$counts
        all(colnames(counts) %in% rownames(Experimental_grops))
        all(colnames(counts) == rownames(Experimental_grops))
        Experimental_grops$Sample_Category=as.numeric(factor(Experimental_grops$Sample_Category))
        mode(counts) <- "integer"
        dds <- DESeqDataSetFromMatrix(countData = counts, colData = Experimental_grops, design = ~ 1)
        ddsF <- dds[ rowSums(counts(dds)) > ncol(dds), ]
        vst=varianceStabilizingTransformation(ddsF)
        normalised_counts <- as.data.frame(assay(vst))
      }else if(norm_method=='NONE'){
        normalised_counts <- y$counts
      }
    } else {
      normalised_counts <- y$counts
    }
}

# Apply inverse normal transformation to each row so traits are normally distributed
if (inverse_normal == TRUE){
  print('Applying inverse normal transformation')
  normalised_counts = quantileNormaliseRows(normalised_counts)
}

if (use_sample_pca) {
  pcs = prcomp(t(normalised_counts), scale = TRUE)  # PCA on samples
  write.table(pcs$x, file = "all__pcs.tsv", sep = "\t")
} else {
  pcs = prcomp(normalised_counts, scale = TRUE)  # PCA on genes
  write.table(pcs$rotation, file = "all__pcs.tsv", sep = "\t")
}

write.table(normalised_counts,file=paste('normalised_phenotype.tsv',sep=''),sep='\t')


pc_var <- pcs$sdev^2
explained_variance <- pc_var / sum(pc_var) * 100
cumulative_variance <- cumsum(explained_variance)

# Limit to first 10 PCs
num_pcs <- 10
explained_variance <- explained_variance[1:num_pcs]
cumulative_variance <- cumulative_variance[1:num_pcs]

# Set up the plotting area
par(mar = c(5, 5, 2, 2))  # Adjust margins for better spacing

pdf("screeplot.pdf")
par(bg = "white")
# First, draw an empty plot to define the plot area and add grid
bar_positions <- barplot(explained_variance, names.arg = 1:num_pcs,
                         col = NA, border = NA, ylim = c(0, 100),
                         xlab = "Principal component", ylab = "Explained variation (%)")
# Add grid lines (horizontal every 10%, vertical at each bar position)
grid(nx = length(bar_positions), ny = 10, col = "lightgray", lty = "dotted")
# Redraw the bar plot on top of the grid
bar_positions <- barplot(explained_variance, names.arg = 1:num_pcs,
                         col = "dodgerblue", border = NA, ylim = c(0, 100),
                         xlab = "Principal component", ylab = "Explained variation (%)", add = TRUE)
# Overlay cumulative variance as red line
lines(bar_positions, cumulative_variance, col = "red", type = "b", pch = 16, lwd = 2)
dev.off()

pdf("biplot.pdf", width = 8, height = 6) 

# Extract variance explained
pc_variance <- summary(pcs)$importance[2,] * 100  # Convert proportion to percentage

# Create PCA plot without metadata
p <- autoplot(pcs, scale = 0, loadings = FALSE, alpha = 0.9, size = 1.2) +
  labs(
    x = paste0("PC1 (", round(pc_variance[1], 2), "% variation)"),
    y = paste0("PC2 (", round(pc_variance[2], 2), "% variation)"),
    title = "PCA Biplot"
  ) +
  theme_minimal() +  # Enables grid lines
  theme(
    text = element_text(size = 5),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey80", size = 0.4),
    panel.grid.minor = element_line(color = "grey90", size = 0.2),
    legend.position = "none"
  )

# Label only the top N most extreme points
num_labels <- 20
scores <- as.data.frame(pcs$x[, 1:2])
top_samples <- rownames(scores)[order(abs(scores$PC1) + abs(scores$PC2), decreasing = TRUE)[1:num_labels]]
scores_subset <- scores[top_samples, ]

# Add labels with smaller font
p <- p + 
  geom_text_repel(data = scores_subset, aes(x = PC1, y = PC2, label = rownames(scores_subset)), 
                  size = 1, max.overlaps = 25)

print(p)
dev.off()


pdf("loadings.pdf") 
pca=pcs
loadings <- as.data.frame(pca$rotation)
loadings$Variable <- rownames(loadings)
loadings_melted <- reshape2::melt(loadings, id.vars = "Variable", variable.name = "PrincipalComponent", value.name = "LoadingScore")
loadings_melted <- loadings_melted %>% group_by(PrincipalComponent) %>% arrange(desc(abs(LoadingScore))) %>% slice(1:5)

ggplot(loadings_melted[loadings_melted$PrincipalComponent %in% c(paste0("PC",1:9)),], aes(x = Variable, y = LoadingScore)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~PrincipalComponent, scales = "free", ncol = 3) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title = "PCA Loading Scores - top 5 genes",
    x = "Gene",
    y = "Loading Score"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()