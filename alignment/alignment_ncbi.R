library(Rsubread)
library(stringr)

#NCBI_m39

setwd("/home/ignacio/Master/TFM/code/Data/ncbi_m39/")

reads_folder <- "/home/ignacio/Master/TFM/code/Data/reads/"

file_name <- list.files(reads_folder)
# Extracting sample names from file
reads_name <- str_extract(file_name, "\\w+(?=\\.)")

# Loop for creating bam files
n <- 1
for (mouse in reads_name){
  align(index = "ncbi_m39",
        readfile1 = paste0(reads_folder,file_name[[n]]),
        output_file = paste0("./bam/",reads_name[[n]],".bam"),
        nthreads = 14
  )
  n <- n + 1
}
