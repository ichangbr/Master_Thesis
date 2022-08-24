# _ ____ ___ 
# | |___ |==]

# HUMAN GSVA

# Initialization
setwd("~/Master/TFM/code/Characterization/")
library(GSVA)
library(tidyverse)
library(pheatmap)
library(factoextra)
library(stringr)
library(rjson)

# Load and clean data
df <- read.delim("./Data/RNAseq_gene_expression_562.txt", sep = "\t")
rownames(df) <- df$Gene
df <- df %>% select(-c(Gene, Accession, Gene_ID))

patient_data <- read.delim("./Data/Schmitz.csv", sep = ";")
patient_ID <- patient_data %>% pull(dbGaP.subject.ID)

df <- df[,names(df) %in% patient_ID]
write.csv(df, file = "./Data/RNAseq_expression_clean.csv")
human_cats <- readRDS("./results/human_categories.rds")

# Generate GSVA
set.seed(1)
gsva_human <- gsva(as.matrix(df), human_cats, kcdf = "Gaussian", parallel.sz = 20, verbose=T)
row_names <- str_extract(rownames(gsva_human), "(?<=^[^_]{0,50}_)\\w+")
prov <- gsva_human
rownames(prov) <- gsub("_"," ",row_names)
order_rows <- fromJSON(file = "./results/row_order.json")
pheatmap(prov[order_rows,],
         treeheight_row = 0,
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         fontsize = 5)

write.csv(gsva_human, "./results/gsva_human.csv")
