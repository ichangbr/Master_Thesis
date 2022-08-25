# _ ____ ___ 
# | |___ |==]

# Results Analysis

setwd("~/Master/TFM/code/Characterization/")
library(biomaRt)
library(tidyverse)
library(factoextra)

# Read data
rnaseq <- read.csv("./Data/RNAseq_expression_clean.csv")
rownames(rnaseq) <- rnaseq$X
rnaseq <- rnaseq %>% select(-c(X))

cats <- readRDS("./results/new_categories.rds")
IFG_genes <- cats$HALLMARK_INTERFERON_GAMMA_RESPONSE
IL10_genes <- cats$REACTOME_INTERLEUKIN_10_SIGNALING

results <- read.csv("./results/Schmitz_w_results.csv")
pBIC <- results$dbGaP.subject.ID[results$category == 1]
pBIC10 <- results$dbGaP.subject.ID[results$category == 2]

edgeR <- read.csv("../analysis_counts/results/edgeR_gene_table.csv")
IFG_res <- edgeR %>% filter(X %in% IFG_genes)
IL10_res <- edgeR %>% filter(X %in% IL10_genes)

IFG_genes <- IFG_res$X[order(desc(abs(IFG_res$logFC)))[1:20]]
IL10_genes <- IL10_res$X[order(desc(abs(IL10_res$logFC)))[1:20]]

genes <- unique(c(IFG_genes, IL10_genes))

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

genes <- mouse_to_human(genes)
genes <- intersect(genes, rownames(rnaseq))

# Select patients from categories
rnaseq_pBIC <- t(rnaseq[genes, pBIC])
rnaseq_pBIC10 <- t(rnaseq[genes, pBIC10])

wilcox <- function(gene){
  a <- wilcox.test(as.numeric(rnaseq_pBIC[,gene]), as.numeric(rnaseq_pBIC10[,gene]))
  return(a$p.value)
}

p_vals <- c()
sig_genes <- c()
for (gene in genes) {
  p <- wilcox(gene)
  p_vals <- c(p_vals, p)
  if (p < 0.05){
    sig_genes <- c(sig_genes, gene)
  }
}

sum(p_vals < 0.05)

# Get gene description

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genedesc <- getBM(attributes = c("hgnc_symbol", "external_gene_name", "description"), filters = "external_gene_name", values = sig_genes ,mart = ensembl)

genedesc <- genedesc %>% select(-c(hgnc_symbol))

saveRDS(genedesc, "./results/gene_description.rds")

# Save genes
library(jsonlite)

p_vals <- p_vals[p_vals < 0.05]
names(p_vals) <- sig_genes

write_json(sig_genes, "./results/sig_genes.json")
saveRDS(p_vals, "./results/pvalues_genes.rds")


# Visualize results
library(pheatmap)

categorization = results$category[results$dbGaP.subject.ID %in% results$dbGaP.subject.ID[results$category %in% c(1,2)]]
names = results$dbGaP.subject.ID[results$dbGaP.subject.ID %in% results$dbGaP.subject.ID[results$category %in% c(1,2)]]

rnaseq_filtered <- t(rbind(rnaseq_pBIC, rnaseq_pBIC10))

cat_col <- data.frame(ID = names,
                      category = categorization)
cat_col <- cat_col[match(colnames(rnaseq_filtered[sig_genes,]), cat_col$ID), ]
rownames(cat_col) <- cat_col$ID
cat_col <- cat_col %>% select(-c(ID))
cat_col$category <- with(cat_col,
                         ifelse(category == 1, "pBIC",
                                ifelse(category == 2, "pBIC10",
                                       "Control")))

a <- pheatmap(rnaseq_filtered[sig_genes, ],
         cluster_rows = F,
         cluster_cols = T,
         fontsize = 5,
         show_colnames = F,
         cutree_cols = 2,
         annotation_col = cat_col
         )

rnaseq_filtered <- t(rnaseq_filtered[sig_genes,])
set.seed(1)
km <- kmeans(rnaseq_filtered, centers = 2)

km$cluster[km$cluster == 1] <-  "pBIC10"
km$cluster[km$cluster == 2] <-  "pBIC"

b <- fviz_cluster(km, rnaseq_filtered, labelsize = 0, shape = 16,
             show.clust.cent = F) +
  labs(title = "")
b
mean(km$cluster==cat_col$category)
library(cowplot)

plot_grid(as.grob(a),b, labels = c("A","B"), label_size = 10)
