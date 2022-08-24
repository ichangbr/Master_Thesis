# _ ____ ___ 
# | |___ |==]

# Profiles

set.seed(1)
library(tidyverse)
setwd("~/Master/TFM/code/ClusterProfiler/")

p_val <- 0.05
lfc <- 1.5

# Load intersection from all methods
library(rjson)
intersection <- fromJSON(file = "~/Master/TFM/code/analysis_counts/results/intersection.json")

# Download msig needed collections
library(clusterProfiler)
library(enrichplot)
#library(upsetplot)
library(msigdbr)
library(stringr)

species <- "Mus musculus"
msig <- list()
msig[["H"]] <- msigdbr(species = species, category = "H")

msig[["C2.KEGG"]] <- msigdbr(
  species = species,
  category = "C2",
  subcategory = "CP:KEGG"
)

msig[["C2.react"]] <- msigdbr(
  species = species,
  category = "C2",
  subcategory = "CP:REACTOME"
)

msig[["GO.BP"]] <- msigdbr(
  species = species,
  category = "C5",
  subcategory = "GO:BP"
)

msig[["GO.CC"]] <- msigdbr(
  species = species,
  category = "C5",
  subcategory = "GO:CC"
)

msig[["GO.MF"]] <- msigdbr(
  species = species,
  category = "C5",
  subcategory = "GO:MF"
)


msig.dfs <- list()
for (name in names(msig)){
  msig.dfs[[name]] <- msig[[name]] %>% select(gs_name, gene_symbol) %>% 
    as.data.frame()
}

# Saving the dataframes for later use
for (name in names(msig.dfs)){
  file_name = paste0("./msig_dfs/",str_replace_all(name,"\\.","_"),".csv")
  write.csv(msig.dfs[[name]], file_name, row.names = F)
}

# Over representation analysis
msig.results <- list()

for (name in names(msig.dfs)) {
  msig.results[[name]] <- enricher(
    gene = intersection,
    TERM2GENE = msig.dfs[[name]],
    pvalueCutoff = p_val
  )
}

saveRDS(msig.results, "./results/ORA_results.rds")
# Preparation for GSEA

edgeR <- read.csv("~/Master/TFM/code/analysis_counts/results/edgeR_gene_table.csv")
logfc <- edgeR$logFC
names(logfc) <- edgeR$X

logfc <- na.omit(logfc)
logfc <- sort(logfc, decreasing = T)

gsea.results <- list()
for (name in names(msig.dfs)) {
  gsea.results[[name]] <- GSEA(
    logfc,
    TERM2GENE = msig.dfs[[name]]
  )
}

saveRDS(gsea.results, "./results/gsea_results.rds")
# Plots for gsea and ORA
gsea.plots_a <- list()
for (name in names(gsea.results)){
  gsea.plots_a[["gseaplot2"]][[name]] <- gseaplot2(gsea.results[[name]], 1:7, base_size = 5)
  gsea.plots_a[["dotplot"]][[name]] <- dotplot(gsea.results[[name]], label_format = 60, font.size = 5)
  # ggsave(paste0('gsea_dotplot',name,'.png'),gsea.plots$dotplot[[name]], width = 10, height = 10)
  # ggsave(paste0('gseaplot2_',name,".png"),gsea.plots$gseaplot2[[name]], width = 10, height = 7.5)
}

library(cowplot)

plotting <- function(data){
    plt <- ggplot(data[1:length(data$Description)], showCategory=10,
         aes(NES, fct_reorder(Description, NES), fill=p.adjust)) +
    geom_col() +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),
                         guide=guide_colorbar(reverse=TRUE)) +
    xlab("Normalized Enrichment Score") +
    ylab(NULL) +
    theme(text = element_text(size = 5))
    return(plt)
}

titles <- c("MSIG Hallmarks", "KEGG Pathways", "Reactome Pathways",
            "Biological Processes Ontology", "Cellular component Ontology",
            "Molecular Function Ontology")

gsea.plots <- list()

for (n in 1:length(gsea.results)) {
  gsea.plots[[n]] <- plotting(gsea.results[[n]])
}

gsea.plots[[1]]
gsea.plots[[2]]
gsea.plots[[3]]
gsea.plots[[4]] <- ggplot(gsea.results$GO.BP[1:20], showCategory=10,
                          aes(NES, fct_reorder(Description, NES), fill=p.adjust)) +
  geom_col() +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),
                       guide=guide_colorbar(reverse=TRUE)) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  theme(text = element_text(size = 5))
gsea.plots[[4]]
gsea.plots[[5]]
gsea.plots[[6]]

plot_grid(gsea.plots_a$gseaplot2$H,gsea.plots[[1]], labels = c("A","B"), label_size = 7)

ora.plots <- list()
for (n in  1:length(msig.results))  {
  ora.plots[[n]] <- dotplot(msig.results[[n]], showCategory = length(msig.results[[n]]$Description),
                            title = paste("Overrepresented",titles[n]))
}

ora.plots[[1]]
ora.plots[[2]]
ora.plots[[3]]
ora.plots[[4]] <- dotplot(msig.results[[4]], showCategory = 10,
                          title = paste("Overrepresented",titles[4]))
ora.plots[[4]]
ora.plots[[5]] <- dotplot(msig.results[[5]], showCategory = 10,
                          title = paste("Overrepresented",titles[5]))
ora.plots[[5]]
ora.plots[[6]] <- dotplot(msig.results[[6]], showCategory = 10,
                          title = paste("Overrepresented",titles[6]))
ora.plots[[6]]

