#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

# Adapted from script in https://github.com/andersonlab/sc_nf_diffexpression/

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option(c('-i', '--iegenes'),
                        type = 'character',
                        default = '',
                        help = 'TSV containing summary statistics. The
                          following columns are required:'
  ),

  optparse::make_option(c('-g', '--group_var'),
                        type = 'character',
                        default = 'coef_value',
                        help = 'Variable to group DE results by.'
  ),

  optparse::make_option(c('-r', '--ranking_var'),
                        type = 'character',
                        default = 'test_statistic',
                        help = 'Column of de_result to rank the genes by.'
  ),

  optparse::make_option(c('-q', '--sample_size'),
                        type = 'integer',
                        default = 101,
                        help = 'Sample size to use during error calculation.'
  ),

  optparse::make_option(c('-l', '--min_set_size'),
                        type = 'integer',
                        default = 1,
                        help = 'Minimum number of genes for a gene set.'
  ),

  optparse::make_option(c('-m', '--max_set_size'),
                        type = 'double',
                        default = Inf,
                        help = 'Maximum number of genes for a gene set.'
  ),

  optparse::make_option(c('-e', '--eps'),
                        type = 'double',
                        default = 1e-10,
                        help = 'Boundary for calculating the p-value.'
  ),

  optparse::make_option(c('-s', '--score_type'),
                        type = 'character',
                        default = 'std',
                        help = 'GSEA score type. Possible options: ("std",
                        "pos", "neg")'
  ),

  optparse::make_option(c('-t', '--unsigned_ranking'),
                        action = 'store_true',
                        default = FALSE,
                        help = 'If `ranking_var` should be signed.'
  ),

  optparse::make_option(c('-z', '--gsets_gene_matrix'),
                        type = 'character',
                        default = '',
                        help = 'TSV file containing gene sets with genes.'
  ),

  optparse::make_option(c('-y', '--gsets_info_file'),
                        type = 'character',
                        default = '',
                        help = 'TSV file containing gene set info.'
  ),

  optparse::make_option(c('-d', '--database'),
                        type = 'character',
                        default = 'all',
                        help = 'Comma-separated list of databases to use from
                        the iDEA humanGeneSets compilation. Descriptions for
                        databas can be found:
                        https://www.gsea-msigdb.org/gsea/msigdb/index.jsp

                        Options:
                        Key Name        Description
                        c2.cgp          Chemical and genetic perturbations
                        c2.cp.biocarta  BioCarta
                        c2.cp.kegg      KEGG
                        c2.cp.reactome  Reactome
                        c2.cp           PID
                        c5.bp           GO biological process
                        c5.cc           GO cellular component
                        c5.mf           GO molecular function
                        c6.all          Oncogenic signatures
                        c7.all          Immunologic signatures
                        all             All gene sets'
  ),

  optparse::make_option(c('-o', '--output_file'),
                        type = 'character',
                        default = '',
                        help = 'Base output name.'
  ),

  optparse::make_option(c('-x', '--exclude_genes'),
                        type = 'character',
                        default = '',
                        help = 'Genes to filter out.'
  ),

  optparse::make_option(c('-n', '--n_cores'),
                        type = 'integer',
                        default = 1,
                        help = 'Number of cores to use.'
  ),

  optparse::make_option(c('-v', '--verbose'),
                        action = 'store_true',
                        default = TRUE,
                        help = ''
  )
)

parser <- optparse::OptionParser(
  usage = '%prog',
  option_list = optionList,
  description = paste0(
    'Perform gene set enrichment using fGSEA.'
  )
)

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
  optionStrings <- character()
  for (item in parserObj@options) {
    optionStrings <- append(optionStrings,
                            c(item@short_flag, item@long_flag))
  }
  optionStrings
}

optStrings <- getOptionStrings(parser)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)

##############
# Variables
##############

# ranking_var <-  'b_gi' 
# sample_size <-  101 
# score_type <-  'std'
# min_set_size <-  1 
# max_set_size <-  Inf 
# eps <-  0 
# verbose <- TRUE
# gsets_gene_matrix <-'assets/data/gene_set_gene_matrix.tsv.gz' 
# gsets_info_file <-  'assets/data/gene_set_info.tsv.gz' 
# database <-  'c2.cp.reactome'
# unsigned_ranking <- TRUE
# n_cores <-  8 
# exclude_genes <- ''


################################################################################

######################## Required Packages #####################################
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(Matrix))
################################################################################

################################ Functions #####################################

get_unique_symbols <- function(df) {
  return(df %>%
           group_by(gene_symbol) %>%
           group_split() %>%
           lapply( . , function(x) {
             ## If more than one unique value, get random and return
             if (nrow(x) > 1) {
               print(sprintf(paste('Gene symbol %s has more than 1 Ensembl ID.',
                                   'Selecting one at random to carry through.'),
                             x$gene_symbol[1]))
               index <- sample(1:nrow(x), 1)
               return(x[index,])
             }
             return(x)
           }) %>%
           do.call(rbind, .)
  )
}

################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file_base <- arguments$options$output_file
n_cores <- arguments$options$n_cores
unsigned_ranking <- arguments$options$unsigned_ranking
score_type <- arguments$options$score_type
eps <- arguments$options$eps
max_set_size <- arguments$options$max_set_size
min_set_size <- arguments$options$min_set_size
sample_size <- arguments$options$sample_size
ranking_var <- arguments$options$ranking_var
group_var <- arguments$options$group_var
gsets_gene_matrix <- arguments$options$gsets_gene_matrix
gsets_info_file <- arguments$options$gsets_info_file
database <- arguments$options$database
exclude_genes <- arguments$options$exclude_genes
iegenes <- arguments$options$iegenes




if (n_cores != 1) {
  n_cores <- n_cores - 1
}

# Read in data
if (verbose) {
  print('Reading in the data...')
}
# Read in exclusion genes
if (exclude_genes != '') {
  exclude.genes <- read_tsv(exclude_genes)
} else {
  exclude.genes <- data.frame(ensembl_gene_id = character(0))
}

# Read in summary statistics
sumstat.df <- read_tsv(iegenes) %>%
  filter(!phenotype_id %in% exclude.genes$ensembl_gene_id) %>%
  mutate(gene_symbol = mapIds(org.Hs.eg.db, keys = phenotype_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")) %>%
  drop_na(gene_symbol)


# Get annotation data
gene_set_genes <- read.csv(gsets_gene_matrix,
                           sep='\t',
                           header=T,
                           row.names='gene')

gene_set_info <- read.csv(gsets_info_file,
                          sep='\t',
                          header=T)

rownames(gene_set_info) <- gene_set_info$gene_set ## Index by gset for later
if (database == 'all') {
  annot_data <- gene_set_genes
} else {
  if (verbose) {
    print(sprintf('Subsetting annotations to those apart of %s pathways...',
                  database))
  }
  databases <- strsplit(x = database,
                        split = ',',
                        fixed = T)[[1]]
  gs_ids <- gene_set_info$gene_set[
    which(gene_set_info$category %in% databases)
  ]
  annot_data <- gene_set_genes[ , as.character(gs_ids)]
  if (verbose) {
    print(sprintf('After subsetting, there are %s gene sets remaining.',
                  ncol(annot_data)))
  }
}

## Format annotation data to be compatible with fGSEA
annot_data <- sapply(colnames(annot_data), function(x) {
  return(list(
    rownames(annot_data)[which(annot_data[[x]] == 1)]
  ))
})

x <- get_unique_symbols(sumstat.df)

## Rank genes by ranking_var
if (unsigned_ranking == FALSE){
  summary_data <- x[[ranking_var]]
} else {
  summary_data <- abs(x[[ranking_var]])
}
names(summary_data) <- x$gene_symbol
summary_data <- na.omit(summary_data)


if (verbose) {
  print('Calculating gene set enrichments for coefficient')
  cat(sprintf(paste('Performing fGSEA using the arguments:\n',
                    'Sample size parameter: %s\n',
                    'Maximum gene set size: %s\n',
                    'Minimum gene set size: %s\n',
                    'EPS value: %s\n',
                    'Score type: %s\n'),
              sample_size,
              max_set_size,
              min_set_size,
              eps,
              score_type))
}

## This will run fgsea.Multilevel, which gives more precise p-vals at a
## consequence for longer compute time. Could adjust to using fgsea.Simple
## See: https://github.com/ctlab/fgsea/blob/master/R/fgsea.R
fgseaRes <- fgsea::fgseaMultilevel(pathways = annot_data,
                                    stats = summary_data,
                                    sampleSize = sample_size,
                                    eps = eps,
                                    minSize = min_set_size,
                                    maxSize = max_set_size,
                                    scoreType = score_type,
                                    nproc = n_cores,
                                    gseaParam = 1,
                                    BPPARAM = NULL)

## Rename fgsea results ot match iDEA
fgseaRes <- fgseaRes %>%
  dplyr::rename(
    annot_id = pathway,
    annot_coef = NES,
    pvalue = pval,
    qvalue_bh = padj, ## fgsea automatically does BH correction
    n_gene_contained = size
  )

## Add shared run-specific information back to result
fgseaRes$gsea_method <- 'fgsea-multi' ## in case we add beta effect model
fgseaRes$sample_size <- sample_size
fgseaRes$min_set_size <- min_set_size
fgseaRes$max_set_size <- max_set_size
fgseaRes$eps <- eps
fgseaRes$score_type <- score_type
fgseaRes$signed_ranking <- !unsigned_ranking

## Add categorical information for each id
fgseaRes$gset_database <- gene_set_info[fgseaRes$annot_id, 'category']
fgseaRes$gset_size <- Matrix::colSums(gene_set_genes[ , fgseaRes$annot_id])

if (verbose) {
  print('Writing GSEA results...')
}

saveRDS(fgseaRes, file=sprintf('%s-fgsea.rds', output_file_base))
write_tsv(x=fgseaRes,
          file = sprintf('%s-gsea_results.tsv.gz', output_file_base),
          col_names=TRUE)
