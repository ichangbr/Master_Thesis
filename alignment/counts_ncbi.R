# _ ____ ___ 
# | |___ |==]

#ncbi_39 feature counts

library(Rsubread)
library(tidyverse)
library(stringr)

setwd("/home/ignacio/Master/TFM/code/Data/ncbi_m39/")

bam_files <- list.files("./bam/", pattern = ".bam$")

gff_file_name <- list.files(pattern = ".gff.*$") #get GFF or GFF3 file name

fc <- featureCounts(files=paste0("./bam/",bam_files), annot.ext= gff_file_name,
                    isGTFAnnotation = TRUE, GTF.featureType = "exon", GTF.attrType = "gene",
                    #isPairedEnd = FALSE, # single-end
                    nthreads = 14, strandSpecific = 0, countMultiMappingReads = TRUE,
                    useMetaFeatures = TRUE)

colnames(fc[[1]]) <- str_extract(colnames(fc[[1]]),".+(?=\\.)") #remove ".bam" from name

counts <- fc$counts

write.csv(counts, file = "raw_counts.csv") #write csv
