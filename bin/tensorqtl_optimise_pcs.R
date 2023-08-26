library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

sumstats.dir = args[1]
alpha = as.numeric(args[2])

##################
# Read in files
##################

sumstats.dir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/results/TensorQTL_eQTLS/category__machine-B_Cell-dMean/'
alpha <- 0.05
annot.name <- basename(sumstats.dir)

sumstat.files <- list.files(sumstats.dir, pattern="*Cis_eqtls_qval.tsv", full.names=FALSE, recursive=TRUE)

list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
  # Get info from filename
  split.file.name <- str_split(file.name,'/')[[1]]
  n.pcs <- split.file.name[1]
  
  
  # Read file
  df <- read_tsv(paste0(sumstats.dir, file.name),show_col_types = FALSE)
  
  # Add new columns
  mutate(df, num_PCs = as.integer(gsub('pcs', '', n.pcs)))
})

sumstat.df <- bind_rows(list.of.sumstat.dfs)

# Remove all non-significant hits
sig.sumstat.df <- sumstat.df %>% 
  filter(qval < alpha)

###############
# Make PC figure
###############

# Wrangle
pc.sig.sumstat.df <- sig.sumstat.df  %>% 
  group_by(num_PCs) %>% 
  add_count(name = 'num_eGenes') %>%
  ungroup() %>% 
  mutate(is_max_egenes = case_when(num_eGenes == max(num_eGenes) ~ TRUE,
                                   TRUE ~ FALSE)) %>% 
  group_by(is_max_egenes) %>% 
  mutate(is_min_PCs = case_when(num_PCs == min(num_PCs) ~ TRUE,
                                TRUE ~ FALSE)) %>% 
  ungroup() %>% 
  mutate(optimal_PC_eGenes = is_min_PCs & is_max_egenes) %>% 
  select(-is_min_PCs, -is_max_egenes)

# Only data needed for plotting PC plot
plot.pc.sig.sum.df <- pc.sig.sumstat.df %>%
  group_by(num_PCs) %>%
  slice(1) %>% 
  ungroup() %>%
  select(num_PCs, num_eGenes, optimal_PC_eGenes)

plot.optimal.pc.sig.sum.df <- plot.pc.sig.sum.df %>% 
  mutate(mid_PCs = mean(c(min(num_PCs),max(num_PCs)) ),
         mid_eGenes = mean(c(min(num_eGenes), max(num_eGenes)) )) %>%
  filter(optimal_PC_eGenes == TRUE) %>%
  slice(1)

###############
# Make PC figure and save
###############

ggplot(plot.pc.sig.sum.df, aes(x=num_PCs, y=num_eGenes)) +
  geom_text(data = plot.optimal.pc.sig.sum.df,
            aes(label=num_PCs, x=mid_PCs, y=mid_eGenes),
            alpha = 0.3, size = 30, colour = 'red')+ # Add text to show which PC is the optimum
  geom_point() +
  geom_point(data=plot.optimal.pc.sig.sum.df,
             colour='red') +
  labs(x = "Number of PCs", y = "Number of eGenes", title = annot.name) +
  theme_bw()

out.dir <- paste0(sumstats.dir, '/optim/')

ggsave(paste0(out.dir,"-optimise_nPCs-FDR",gsub('[.]', 'pt', alpha),".pdf"),
         width = 3, height = 3)

# Save a file containing name and optimal PCs
write_tsv(plot.pc.sig.sum.df, paste0(out.dir,"-optimise_nPCs-FDR",gsub('[.]', 'pt', alpha)".txt"))

