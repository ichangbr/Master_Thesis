# _ ____ ___ 
# | |___ |==]

# pBIC10 & pBIC characterization

setwd("~/Master/TFM/code/Characterization/")
set.seed(1)

library(GSVA)
library(tidyverse)
library(edgeR)
library(reshape2)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(factoextra)
library(jsonlite)
library(cowplot)
# Read and clean data
cts_clean <- read.csv("~/Master/TFM/code/analysis_counts/results/cts_clean.csv")
rownames(cts_clean) <- cts_clean$X
#cts_clean <- cts_clean %>% select(-c(X,YC_1669,YC_1067,YC_1671,YC_1070))
cts_clean <- cts_clean %>% select(-c(X))

# Check if data has no immunoglobulins
sum(str_detect(rownames(cts_clean), "Ig[a h k l]{1}[d e g m v]{0,1}\\w+")) 

# Pre-processing
# Remove rows with less than 2 columns with less than 10 counts
keep <- rowSums(cts_clean >= 10) >= 2
cts_clean <- cts_clean[keep,]

# Generate logCPMs
cpms <- cpm(cts_clean, log = T) 

# violin plot of data
cts_melt <- as.data.frame(cpms)
cts_melt["X"] <- rownames(cts_melt)
cts_melt <- melt(cts_melt, id.vars = "X")
cts_melt$variable <- as.factor(cts_melt$variable)

a <- cts_melt %>% ggplot(aes(x = variable, y = value)) +
  geom_violin() +
  scale_y_log10("Counts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        legend.position = "none") +
  labs(x = "Samples")

b <- cts_melt %>% ggplot(aes(x = value, color = variable)) +
  geom_density() +
  theme(legend.position = "none")

plot_grid(a,b, labels = c("A", "B"), label_size = 10)

# Generate lists of genes for the selected categories
selected_cats <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING",
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                   "HALLMARK_GLYCOLYSIS",
                   "HALLMARK_MTORC1_SIGNALING",
                   "GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE",
                   "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
                   "GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE",
                   "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
                   "GOCC_PHAGOCYTIC_VESICLE",
                   "GOCC_T_CELL_RECEPTOR_COMPLEX",
                   "GOCC_RECEPTOR_COMPLEX",
                   "GOMF_ANTIGEN_BINDING",
                   "GOMF_PEPTIDE_ANTIGEN_BINDING",
                   "GOMF_IMMUNE_RECEPTOR_ACTIVITY",
                   "GOMF_SIGNALING_RECEPTOR_REGULATOR_ACTIVITY",
                   "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
                   "KEGG_SELENOAMINO_ACID_METABOLISM",
                   "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
                   "KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES",
                   "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
                   "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
                   "REACTOME_INTERLEUKIN_10_SIGNALING",
                   "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
                   "REACTOME_PD_1_SIGNALING"
)

# selected_cats <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING",
#                    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#                    "GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE",
#                    "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
#                    "GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE",
#                    "GOCC_PHAGOCYTIC_VESICLE",
#                    "GOCC_T_CELL_RECEPTOR_COMPLEX",
#                    "GOMF_ANTIGEN_BINDING",
#                    "GOMF_PEPTIDE_ANTIGEN_BINDING",
#                    "GOMF_IMMUNE_RECEPTOR_ACTIVITY",
#                    "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
#                    "REACTOME_INTERLEUKIN_10_SIGNALING",
#                    "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
#                    "REACTOME_TRANSLATION"
# )

# Generating gene lists for the categories

gsea.results <- readRDS("../ClusterProfiler/results/gsea_results.rds")

gene_lists_all <- list()

for (name in names(gsea.results)){
  gene_lists_all[[name]] <- geneInCategory(gsea.results[[name]])
}

new_cats <- list()
for (cat in selected_cats){
  gene_set <- str_extract(cat,"^([^_])+") # Capturing the first portion of the string
  if (gene_set == "HALLMARK"){
    gene_set <- "H"
  } else if (gene_set == "GOBP") {
    gene_set <- "GO.BP"
  } else if (gene_set == "GOCC") {
    gene_set <- "GO.CC"
  } else if (gene_set == "GOMF") {
    gene_set <- "GO.MF"
  } else if (gene_set == "KEGG") {
    gene_set <- "C2.KEGG"
  } else if (gene_set == "REACTOME") {
    gene_set <- "C2.react"
  }
  genes_prov <- gene_lists_all[[gene_set]]
  new_cats[[cat]] <- genes_prov[[cat]]
}


set.seed(1)
gsva_res <- gsva(cpms, new_cats, kcdf = "Gaussian", parallel.sz = 20, verbose=T)
# Heatmap results
prov <- gsva_res
row_names <- str_extract(rownames(gsva_res), "(?<=^[^_]{0,50}_)\\w+")
row_names <- gsub("_", " ", row_names)
rownames(prov) <- row_names
heatmap <- pheatmap(prov,
         treeheight_row = 0,
         show_rownames = T,
         show_colnames = T,
         scale = "row",
         fontsize = 5,
         cutree_cols = 3
         )
write_json(heatmap$tree_row$order, "./results/row_order.json")
# k-means to confirm the clustering works to determine the group
k <- 3
km <- kmeans(t(gsva_res), k)
fviz_cluster(km, data = t(gsva_res)) +
  labs(title = "") +
  theme(legend.position = "none")

# Save important files
# GSVA results
write.csv(gsva_res, "./results/gsva_results.csv")
#categories vector
write_json(selected_cats, "./results/selected_categories.json")
saveRDS(new_cats, "./results/new_categories.rds")
