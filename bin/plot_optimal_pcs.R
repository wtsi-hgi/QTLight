#!/usr/bin/env Rscript

library(ggplot2)
library(stringr)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
working_dir <- ifelse(length(args) >= 1, args[1], ".")
alpha <- ifelse(length(args) >= 2, as.numeric(args[2]), 0.05)
annot_name <- ifelse(length(args) >= 3, args[3], "auto_detect")
out_dir <- ifelse(length(args) >= 4, args[4], "./output")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Get files matching __merged_cis_scores.tsv.gz
files <- list.files(path = working_dir,
                    pattern = "__merged_cis_scores.tsv.gz$",
                    full.names = TRUE)

if (length(files) == 0) {
  stop("No __merged_cis_scores.tsv.gz files found in the directory.")
}

cat("Found", length(files), "files\n")

# Load and tag each file with its PC count
dfs <- lapply(files, function(f) {
  pcs <- str_extract(f, "\\d+pcs") |> str_remove("pcs") |> as.integer()
  df <- read_tsv(f, show_col_types = FALSE)
  df$num_PCs <- pcs
  return(df)
})

df_all <- bind_rows(dfs)

# Filter by q-value threshold
sig_df <- df_all %>% filter(qval < alpha)

# Count all unique PCs seen in input
all_pcs <- sort(unique(df_all$num_PCs))

if (nrow(sig_df) == 0) {
  message("⚠️ No significant hits found — continuing with 0 eGenes per PC.")

  pc_summary <- data.frame(
    num_PCs = all_pcs,
    num_eGenes = 0,
    optimal_PC_eGenes = FALSE
  )

  plot.optimal.pc.sig.sum.df <- NULL
} else {
  pc_summary <- sig_df %>%
    group_by(num_PCs) %>%
    summarise(num_eGenes = n(), .groups = "drop") %>%
    mutate(is_max_egenes = num_eGenes == max(num_eGenes)) %>%
    mutate(optimal_PC_eGenes = is_max_egenes & num_PCs == min(num_PCs[is_max_egenes]))

  plot.optimal.pc.sig.sum.df <- pc_summary %>% filter(optimal_PC_eGenes)
}


# Get annotation name from first filename if not provided
if (annot_name == "auto_detect") {
  annot_name <- str_split(basename(files[[1]]), "__")[[1]][2]
}

# Compute midpoint for label placement
mid_x <- mean(range(pc_summary$num_PCs))
mid_y <- mean(range(pc_summary$num_eGenes))

# Plot
p <- ggplot(pc_summary, aes(x = num_PCs, y = num_eGenes)) +
  geom_point() +
  labs(x = "Number of PCs", y = "Number of eGenes", title = annot_name) +
  theme_bw()

# Only add red highlight and label if we have optimal points
if (!is.null(plot.optimal.pc.sig.sum.df) && nrow(plot.optimal.pc.sig.sum.df) > 0) {
  mid_x <- mean(range(pc_summary$num_PCs))
  mid_y <- mean(range(pc_summary$num_eGenes))

  p <- p +
    geom_text(data = plot.optimal.pc.sig.sum.df,
              aes(x = mid_x, y = mid_y, label = num_PCs),
              size = 10, alpha = 0.3, colour = 'red') +
    geom_point(data = plot.optimal.pc.sig.sum.df, colour = "red")
}

# Output paths
out_base <- file.path(out_dir, paste0("optimise_nPCs-FDR", gsub("\\.", "pt", alpha)))
ggsave(paste0(out_base, ".pdf"), plot = p, width = 4, height = 4)
write_tsv(pc_summary, paste0(out_base, ".txt"))


# Find PC with highest number of eGenes (but allow ties)
# Then among ties, pick one randomly
if (all(pc_summary$num_eGenes == 0)) {
  message("All PCs have 0 eGenes — writing empty file for optimal_PC.")
  write("", file = paste0(out_base, "_optimal_PC.txt"))
} else {
  selected_pc <- pc_summary %>%
    filter(num_eGenes == max(num_eGenes)) %>%
    arrange(num_PCs) %>%
    slice(1) %>%
    pull(num_PCs)

  cat("Selected optimal PC:", selected_pc, "\n")
  writeLines(as.character(selected_pc), paste0(out_base, "_optimal_PC.txt"))
}

cat("Done. Plot and summary saved in", out_dir, "\n")
