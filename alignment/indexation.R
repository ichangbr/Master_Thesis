library(Rsubread)

#Create indexes for ncbi and ensembl genome
setwd("/home/ignacio/Master/TFM/code/Data/ncbi_m39/")
buildindex(basename = "ncbi_m39",
           reference="GCF_000001635.27_GRCm39_genomic.fna")

setwd("/home/ignacio/Master/TFM/code/Data/ensembl_m39/")
buildindex(basename = "ensembl_m39",
           reference = "Mus_musculus.GRCm39.dna.primary_assembly.fa")
