# this script is to have the giant file with all replicates, and split it in files per chromosome
#load required libraries

library(data.table)
library(tidyverse)
library(foreach)
library(doMC)
registerDoMC(cores = 5)

# define input and output directories

file.input.name <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/sure.indel.6rep.RDS"
dir.output <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE_combined/"

# load the giant datafile

sure.indel.combined <- readRDS(file = file.input.name)

print("large dataframe has been read")
# Retrieve the chromosome numbers as one long vector

chrom.idx <- (str_remove_all(string = sure.indel.combined$chrom, pattern = "[_paternalmaternal]"))

# Count rows per chromosome as a test
x <- table(chrom.idx)
saveRDS(object = x, file = paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/tmp/nrowperchrom.RDS"))
 
print(paste(Sys.time(), "starting loop"))

foreach(i = c(1:22,"X")) %dopar%{
  
  print(paste("starting analysis of chromosome", i, "at", Sys.time()))
 
  # 1. Generate a new dataframe for chromosome i only
  
  sure.indel.chrom <- sure.indel.combined[chrom.idx == i,]
  print(paste(Sys.time(),"finished dataframe generation", i))
  
  # 2. Save the new dataframe to a new r object
  
  file.output.name <- paste0(dir.output, "sure.indel.combined.chrom.", i, ".RDS")
  saveRDS(object = sure.indel.chrom, file = file.output.name)
}


print(paste(Sys.time(), "end loop"))
  

