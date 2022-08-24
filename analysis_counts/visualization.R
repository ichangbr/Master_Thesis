# _ ____ ___ 
# | |___ |==]

# DEG visualization

setwd("~/Master/TFM/code/analysis_counts/")

library(DESeq2)
library(edgeR)
library(pheatmap)
library(rjson)
library(gplots)
library(RColorBrewer)
library(stringr)
load("DEG_workplace.RData")

counts <- read.csv("~/Master/TFM/code/Data/ncbi_m39/raw_counts.csv")
row.names(counts) <-  counts[,1]
counts <- counts[,-1]

intersection <- fromJSON(file = "./results/intersection.json")

filtered <- counts[row.names(counts) %in% intersection, ]
filtered <- filtered[, str_detect(names(counts), "pBIC")]
per_mil <- cpm(filtered, log = T)
boxplot(per_mil, notch =  T, las = 2)

dev.off()
mypalette <- brewer.pal(9, "YlOrRd")
pheatmap(per_mil, treeheight_row = 0,treeheight_col = 0, show_rownames = F,
         color = mypalette, legend = T, cellwidth = 25, main = "Heatmap of Differentially\nExpressed Genes")


graph_res <- as.data.frame(results_pBIC_pBIC10)

graph_res$significant <- ifelse(graph_res$padj < .05, "Significant", NA)
graph_res %>% ggplot(aes(baseMean, log2FoldChange, colour = significant))+
  geom_point(size = .1) +
  scale_y_continuous(limits = c(-4,4), oob = scales::squish)+
  scale_x_log10()+
  geom_hline(yintercept = 0, colour = "darkorchid4", size = .4, linetype = "longdash")+
  geom_density2d(colour = "black")+
  theme_bw() +
  theme(legend.position = 'none')

