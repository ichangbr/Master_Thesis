# _ ____ ___ 
# | |___ |==]

# Differential Expression

library(tidyverse)
library(DESeq2)

p_val <- 0.01

setwd("/home/ignacio/Master/TFM/code/analysis_counts")
cts <- read_csv("~/Master/TFM/code/Data/ncbi_m39/raw_counts.csv")

# CHANGING ROW NAMES
cts <- as.data.frame(cts)
colnames(cts)
colnames(cts)[1] <- "names"

# names are the gene names, we don't want any duplicates so we check to
# be able to set them as row names
names <- cts$names
table(duplicated(names))

# We set the gene names as row names
rownames(cts) <- names
cts <- cts[-1]
head(cts)

# CONDITION TABLE
#Column name structure: NAME_CONDITION_NUMBER
#The following regex extracts the CONDITION from the column name
condition <- as.factor(str_extract(colnames(cts),
                                   pattern = "(?<=_)\\w+(?=_)"))
#The following regex extracts the portion CONDITION_NUMBER
names <- str_extract(colnames(cts), pattern = "(?<=_)\\w+")
colnames(cts) <- names #changing colum names to the condition and number
cond_data <- data.frame(row.names = names,
                        condition = as.factor(condition))

# Remove immunoglobulins
cts_clean <- cts[!str_detect(rownames(cts), "Ig[a h k l]{1}[d e g m v]{0,1}\\w+"),]

#Count how many NAs in each column
colSums(is.na(cts_clean)) # There are no NAs so there is no need to remove them

# DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(cts_clean),
                              colData = cond_data,
                              design = ~ condition)
# PRE-FILTERING
keep <- rowSums(counts(dds)) >= 10 #Keeps rows with more than 10 counts
dds <- dds[keep,]

# DEG analysis
dds <- DESeq(dds)
resultsNames(dds)

# Quality check
rld <- rlog(dds,blind = F) # Regularized log transformation
plotPCA(rld, intgroup = "condition")

# (WE DO SOME EXPLORATORY GRAPHS (HEATMAPS), IT IS IN THE MARKDOWN)

# Results pBIC v. pBIC10
results_pBIC_pBIC10 <- results(dds,
                               contrast = c("condition", "pBIC", "pBIC10"),
                               pAdjustMethod = "BH") #Benjamini-Hochberg

summary(results_pBIC_pBIC10)
# Get complete cases from results_pBIC_pBIC10
results_pBIC_pBIC10 <- results_pBIC_pBIC10[complete.cases(results_pBIC_pBIC10),] # THIS
# HAS THE COMPLETE TABLE OF ALL THE GENES WITH logFC and p-value

# Filtering by p-value
list_padj <- results_pBIC_pBIC10$padj
# list_padj <- list_padj[!is.na(list_padj)] # Remove NAs
table(list_padj < p_val)
names_deg <- rownames(results_pBIC_pBIC10[list_padj < p_val,]) #list names DEGs

res <- as.data.frame(results_pBIC_pBIC10[list_padj < p_val,]) #Results as df (THIS HAS
# THE COMPLETE TABLE FOR THE DEGs)

# write.csv(results_pBIC_pBIC10, "./results/DESeq2_gene_table.csv") #Save complete results table

# edgeR
# Remove Immunoglobulins

library(edgeR)

dge <- DGEList(cts_clean) #edgeR object

# Introduce conditions
dge$samples$group <- condition

#counts per million
myCPM <- cpm(cts_clean)

#find equivalent to 10 counts
cpm_breaks <- function(cpm_array, column, threshold){
  y_prov <- cpm_array[,column][cpm_array[,column] <= 3]
  x_prov <- cts_clean[,column][cpm_array[,column] <= 3]
  prov_df <- data.frame(CPM = y_prov, counts = x_prov)
  reg <-  lm(CPM~counts, data = prov_df)
  equivalent <- reg$coefficients[1] + reg$coefficients[2]*threshold
  return(as.numeric(equivalent))
}

num_columns <- length(colnames(myCPM))
mean_thresh <- 0
for (col in 1:num_columns) {
  mean_thresh <- mean_thresh + cpm_breaks(myCPM,col,10)
}
mean_thresh <- mean_thresh/num_columns

keep <- rowSums(myCPM > mean_thresh) >= 2
table(keep)

dge <- dge[keep, keep.lib.sizes = F] #filtered by equivalent to 10
dge <- calcNormFactors(dge) # calculate Normalization factors

# DESIGN MATRIX
design <- model.matrix(~0 + condition)
colnames(design) <- levels(condition)

# Quasi-likelihood F-test (WHAT DOES THIS DO)
dge_disp <- estimateDisp(dge,design = design, robust = T)

# Object for pBIC v. pBIC10
dge_fit <- glmQLFit(dge_disp,design,robust = T)
cont.mat <- makeContrasts(pBIC - pBIC10, levels = design)
pBIC_pBIC10 <- glmQLFTest(dge_fit,contrast = cont.mat)

# Save genes table
# write.csv(pBIC_pBIC10$table, "./results/edgeR_gene_table.csv")

# Top tags
summary(decideTests(pBIC_pBIC10,p.value = p_val, lfc = 1.5)) #filtered by p-value

pBIC_pBIC10.DE <- topTags(pBIC_pBIC10,
                          n = 6000, p.value = p_val)
dim(pBIC_pBIC10.DE)

# Summary of DEGs
summary(decideTests(pBIC_pBIC10, lfc = 1.5)) #filtered by log fold change

# Up regulated
gt2 <- pBIC_pBIC10.DE$table$logFC >= 1.5
pBIC_pBIC10.UPgt2 <- pBIC_pBIC10.DE[gt2,] #SAVE NAMES?
dim(pBIC_pBIC10.UPgt2)

# Down regulated
gt2 <- pBIC_pBIC10.DE$table$logFC <= -1.5
pBIC_pBIC10.DOWNgt2 <- pBIC_pBIC10.DE[gt2,] # SAVE NAMES?
dim(pBIC_pBIC10.DOWNgt2)

# limma
library(limma)
voom_obj <- voom(counts = dge, design = design, plot = F)

pBIC_pBIC10.lm <- lmFit(object = voom_obj, design = design)
pBIC_pBIC10.contrast <- contrasts.fit(fit = pBIC_pBIC10.lm, contrast = cont.mat)
pBIC_pBIC10.bayes <- eBayes(pBIC_pBIC10.contrast)

# Significance test
summary(decideTests(object = pBIC_pBIC10.bayes, lfc = 1.5, p.value = p_val))
pBIC_pBIC10.dt <- topTable(fit = pBIC_pBIC10.bayes, sort.by ="p",
                           n = 2000, lfc = 1.5, p.value = p_val)
summary(pBIC_pBIC10.dt$logFC >= 1.5)

# Up regulated
pBIC_pBIC10.UPreg <- rownames(pBIC_pBIC10.dt[pBIC_pBIC10.dt$logFC > 0,]) #sorted by p-value

# Down regulated
pBIC_pBIC10.DOWNreg <- rownames(pBIC_pBIC10.dt[pBIC_pBIC10.dt$logFC < 0,])

pBIC_pBIC10.DGE <- c(pBIC_pBIC10.UPreg, pBIC_pBIC10.DOWNreg)
pBIC_pBIC10.results <- topTable(fit = pBIC_pBIC10.bayes, sort.by = "p", # All genes sorted by p-value
                                number = 20000)
dim(pBIC_pBIC10.results)

# write.csv(pBIC_pBIC10.results, "./results/limma_gene_table.csv") # Save all genes and statistics obtained with limma

#write.csv(cts_clean, "./results/cts_clean.csv")
