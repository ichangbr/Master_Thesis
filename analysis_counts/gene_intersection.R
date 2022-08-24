# _ ____ ___ 
# | |___ |==]

# Gene Intersection

library(tidyverse)

setwd("~/Master/TFM/code/analysis_counts/")

p_val <- 0.01
lfc <- 1.5

# Load Data
library(stringr)
file_names <- list.files("./results/", pattern = "*.csv")
names_dfs <- str_match(file_names, ".*(?=_gene_table.csv)")
data <- list()

# Modify data
count <- 1
for (file in file_names){
  data[[names_dfs[count]]] <- read.csv(paste0("./results/", file))
  rownames(data[[names_dfs[count]]]) <- data[[names_dfs[count]]][, 1]
  data[[names_dfs[[count]]]] <- data[[names_dfs[count]]][, -1]
  count <- count + 1
}

# Check if clean

for (df in names(data)){
  if (
    sum(str_detect(rownames(data[[df]])
                   ,"Ig[a h k l]{1}[d e g m v]{0,1}\\w+")) == 0
  ) {
    print(TRUE)
  }
}

# Get DEGs names
DEGs <- list()
DEGs[["DESeq2"]] <- data$DESeq2 %>% filter(padj < p_val) %>% filter(abs(log2FoldChange) >= lfc) %>% rownames()
DEGs[["edgeR"]] <- data$edgeR %>%  filter(PValue < p_val) %>% filter(abs(logFC) >= lfc) %>% rownames()
DEGs[["limma"]] <- data$limma %>% filter(adj.P.Val < p_val) %>% filter(abs(logFC) >= lfc) %>% rownames()

# Get intersection
i <- 1
for (l in DEGs) {
  if (i == 1){
    intersection <- l
  } else {
    intersection <- intersect(intersection, l)
  }
  print(length(intersection))
  i <- i + 1
}

# # Save JSON
library(jsonlite)
write_json(intersection, "./results/intersection.json")

library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 = 2212,                          # Disable transparency of venn diagram
                 area2 = 2404,
                 area3 = 1837,
                 n12 = 2161,
                 n23 = 1789,
                 n13 = 1749,
                 n123 = 1744,
                 fill = c("blue", "green", "orange"),
                 alpha = 0.5,
                 lty = "blank",
                 category = c("DESeq2",
                              "edgeR",
                              "limma"))
