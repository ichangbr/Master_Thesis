# _ ____ ___ 
# | |___ |==]

# Plotting Data

setwd("/home/ignacio/Master/TFM/code/ClusterProfiler/")

deseq_res <- read.csv("~/Master/TFM/code/analysis_counts/results/DESeq2_gene_table.csv")
edger_res <- read.csv("~/Master/TFM/code/analysis_counts/results/edgeR_gene_table.csv")
limma_res <- read.csv("~/Master/TFM/code/analysis_counts/results/limma_gene_table.csv")
cts <- read.csv("~/Master/TFM/code/analysis_counts/results/cts_clean.csv")

library(ggplot2)
library(reshape2)

# BOXPLOT OF SAMPLES
cts_melt <- melt(cts, id.vars = "X")
cts_melt$variable <- as.factor(cts_melt$variable)

cts_melt %>% ggplot(aes(x = variable, y = value)) +
  geom_boxplot() +
  scale_y_log10("Counts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        legend.position = "none") +
  labs(x = "Samples")

deseq_res %>% ggplot(aes(baseMean, log2FoldChange, colour = ifelse(padj < .05, "Significant", NA)))+
  geom_point(size = .7) +
  scale_y_continuous(limits = c(-4,4), oob = scales::squish)+
  scale_x_log10()+
  geom_hline(yintercept = 1.5, colour = "darkorchid4", size = .6, linetype = "longdash") +
  geom_hline(yintercept = -1.5, colour = "darkorchid4", size = .6, linetype = "longdash") +
  geom_density2d(colour = "black") +
  labs(col = "Significance") +
  scale_color_discrete(labels = c("Significant", "Not significant")) +
  theme(legend.position = "bottom") +
  labs(x = "Mean of normalized counts", y = "Logarithmic Fold Change")





