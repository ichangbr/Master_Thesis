# _ ____ ___ 
# | |___ |==]

# Most affected gene from IL10 cascade

#Initialization
setwd("~/Master/TFM/code/Characterization/")
library(tidyverse)
library(pheatmap)

cats <- readRDS("results/new_categories.rds")

IL10_genes <- cats$REACTOME_INTERLEUKIN_10_SIGNALING
cts_clean <- read.csv("../analysis_counts/results/cts_clean.csv")
rownames(cts_clean) <- cts_clean$X
cts_clean <- cts_clean %>% select(-c(X, YC_1669, YC_1671, YC_1067, YC_1070))

# Visualize differeces
pheatmap(cts_clean[IL10_genes,],
         scale = "row",
         cluster_rows = F,
         cluster_cols = F)

results <- read.csv("../analysis_counts/results/edgeR_gene_table.csv")

# Extract gene with most difference

rownames(results) <- results$X
results <- results %>% select(logFC) %>%
  filter(., rownames(.) %in% IL10_genes) %>%
  filter(logFC == max(logFC)) %>% 
  rownames()

# Save gene name
print(results)
fileconn <- file("./results/IL10_most_affected")
writeLines(results, fileconn)
close(fileconn)
