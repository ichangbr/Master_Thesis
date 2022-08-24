# _ ____ ___ 
# | |___ |==]

#Quality Control plot

setwd("./Master/TFM/code/alignment/")
library(tidyverse)
library(reshape2)
library(stringr)
library(cowplot)

phred <- read.csv("./Data/fastqc_per_base_sequence_quality_plot.csv")
names(phred)[1] <- "bp"
names(phred)[-c(1)] <- str_extract(names(phred)[-c(1)], "(?<=_).+")

phred_score <- melt(phred, id.vars = "bp") %>% ggplot(aes(bp, value, color = variable)) +
  geom_line() + theme(legend.position = "bottom") + labs(color = "Samples") +
  xlab("Position (bp)")+
  ylab("Phred score")

legend <- get_legend(phred_score)

phred_score <- phred_score + theme(legend.position = "none")

GC <- read.csv("./Data/fastqc_per_sequence_gc_content_plot.csv")
names(GC)[1] <- "GC"

GC_content <- melt(GC, id.vars = "GC") %>% ggplot(aes(GC, value, color = variable)) +
  geom_line() +
  theme(legend.position = "none") +
  xlab("% GC") +
  ylab("Percentage")

dupli <- read.csv("./Data/fastqc_sequence_duplication_levels_plot.csv")
names(dupli)[1] <- "dup"
ticks <- dupli$dup
dupli$dup <- 1:length(dupli$dup)

dupli_levs <- melt(dupli, id.vars = "dup") %>% ggplot(aes(dup, value, color = variable)) +
  geom_line() + 
  xlab("Sequence Duplication Level") +
  ylab("% of Library") +
  scale_x_continuous(breaks = 1:length(ticks), labels = ticks, guide = guide_axis(angle = 90)) +
  theme(legend.position = "none")

low_row<- plot_grid(GC_content, dupli_levs, labels = c("B", "C"), label_size = 10)

all <- plot_grid(phred_score, low_row, labels = c("A"), label_size = 10, nrow = 2)

plot_grid(all, legend, nrow = 2, rel_heights = c(1,.2))
