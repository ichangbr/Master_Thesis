# _ ____ ___ 
# | |___ |==]

# Convert mouse genes to human

# Initialization

setwd("~/Master/TFM/code/Characterization")
set.seed(1)

library(tidyverse)
library(rjson)
library(pheatmap)

# Load Data

new_cats <- readRDS("./results/new_categories.rds")

mouse_to_human <- function(x){

  if (!exists("mouse_human_genes")) { # Loads data first time the function is used
    print("LOADING DATA")
    assign("mouse_human_genes",
           read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t"), # 2022-07-11 08:00 	13M
           envir = .GlobalEnv)
  }

  # Create provisional df for human and mouse (READABILITY)
  mouse_genes <- mouse_human_genes %>% filter(Common.Organism.Name == "mouse, laboratory")
  human_genes <- mouse_human_genes %>% filter(Common.Organism.Name == "human")

  # Get only genes in database
  in_data = x %in% mouse_genes$Symbol
  x <- x[in_data]

  # Get class keys
  class_keys <- mouse_genes %>% filter(Symbol %in% x)
  class_keys <- class_keys[["DB.Class.Key"]]

  # Get gene human gene sybols
  new_gene_names <- human_genes %>% filter(DB.Class.Key %in% class_keys)
  new_gene_names <- unique(new_gene_names[,"Symbol"])

  return(unique(new_gene_names))
}

# Convert to human genes
human_genes <- list()
for (name in names(new_cats)){
  human_genes[[name]] <- mouse_to_human(new_cats[[name]])
}

# Convert gene most differentially expressed gene from IL10 cascade to human
saveRDS(human_genes, "./results/human_categories.rds")

IL10_most <- read.table(file = "./results/IL10_most_affected", sep = ",")[["V1"]]
fileconn <- file("./results/IL10_most_affected_human")
writeLines(mouse_to_human(IL10_most), fileconn)
close(fileconn)
