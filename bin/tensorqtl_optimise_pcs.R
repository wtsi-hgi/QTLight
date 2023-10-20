#!/usr/bin/env Rscript

library(ggplot2)
library(stringr)
library(readr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if(!interactive()){
   pdf(NULL)
}


sumstats.dir = './'
alpha = as.numeric(0.05)
interaction.name = 'NA'
annot.name = 'Celltypist:Immune_All_High:predicted_labels-Monocytes-dMean'
out.dir = './OPTIM_pcs/base_output__base'

sumstats.dir = args[1]
alpha = as.numeric(args[2])
interaction.name = args[3]
annot.name = args[4]
out.dir = args[5]

##################
# Read in files
##################

if (interaction.name == "NA") {
  interaction.name = NA
}

if (is.na(interaction.name)){
  sumstat.files <- list.files(sumstats.dir, pattern="*Cis_eqtls_qval.tsv", full.names=FALSE, recursive=TRUE)
} else {
  sumstat.files <- list.files(sumstats.dir, pattern="*cis_inter1.cis_qtl_top_assoc.txt.gz", full.names=FALSE, recursive=TRUE)
}
list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
  # Get info from filename
  split.file.name <- str_split(file.name,'/')[[1]]
  symlink.dir <- split.file.name[2]
  split.symlink.dir <- str_split(symlink.dir,'__')[[1]]
  n.pcs <- split.symlink.dir[1]
  
  
  # Read file
  df <- read_tsv(paste0(sumstats.dir, file.name),show_col_types = FALSE)
  
  # Add new columns
  mutate(df, num_PCs = as.integer(gsub('pcs', '', n.pcs)))
})

sumstat.df <- bind_rows(list.of.sumstat.dfs)

# Remove all non-significant hits

if (is.na(interaction.name)){
  sig.sumstat.df <- sumstat.df %>% 
    filter(qval < alpha)
} else {
  sig.sumstat.df <- sumstat.df %>% 
    filter(pval_adj_bh < alpha)
}

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
  theme_bw() +
  theme(title=element_text(size=6,face="bold"))


out.path <- paste0(out.dir,"/optimise_nPCs-FDR",gsub('[.]', 'pt', alpha))


# Save plot
ggsave(paste0(out.path,".pdf"),width = 3, height = 3)

# Save a file containing name and optimal PCs
write_tsv(plot.pc.sig.sum.df, paste0(out.path,".txt"))
