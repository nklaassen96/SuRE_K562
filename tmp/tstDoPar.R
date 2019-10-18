# Test the foreach loop

#install.packages("foreach")
library(foreach)

library(tidyverse)

#install.packages("doMC")
library(doMC)
registerDoMC(cores = 20)

dir.input <- "/DATA/usr/n.klaassen/projects/SuRE_K562/tmp/sure.test."

foreach(k=1:2) %dopar% {
  fname <- paste0(dir.input, k, ".RDS")
  counts <- readRDS(fname)
  
  counts.df <- as.data.frame(counts)
  #OR
  #fname <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE42-1/combined/sure.42-1.counts.indel.all.RDS"
  #counts <- load(file = fname)
  
  
  # Remove Rows with SNP_ID == "."
  counts.df <- counts.df[counts.df$SNP_ID != ".",]
  
  # Get the indel.ids of the indels that have reads for reference and alternative alleles and
  # construct a new dataframe with only these indels. 
  
  indel.ids.vector <- names(which(tapply(counts.df$SNP_VAR, counts.df$SNP_ID, function(x) {all(c(0,1) %in% x)} ) == TRUE))
  indel.ids.vector <- indel.ids.vector #Generate subset
  
  indel.df <- counts.df[counts.df$SNP_ID %in% indel.ids.vector,]
  
  
  
  # Generate an empty results dataframe for the indels 
  
  results.indel <- data.frame(matrix(nrow = length(indel.ids.vector), ncol = 17))
  colnames(results.indel) <- c("SNP_ID", "chrom", "pos", "strand", "ref.element.count", "alt.element.count", "K562.cDNA.ref.mean","K562.cDNA.alt.mean", "K562.cDNA.ref.median", "K562.cDNA.alt.median", "HepG2.cDNA.ref.mean","HepG2.cDNA.alt.mean","HepG2.cDNA.ref.median", "HepG2.cDNA.alt.median", "K562.wilcoxon.pvalue", "HepG2.wilcoxon.pvalue", "location.annotation")
  
  for (i in 1:length(indel.ids.vector)){
    if (i %% 1 == 0){print(i)}
    snp.idx <- which(indel.df$SNP_ID == indel.ids.vector[i])
    
    ref <- which(indel.df[snp.idx, "SNP_VAR"] == 0)
    alt <- which(indel.df[snp.idx, "SNP_VAR"] == 1)
    
    # Construct results dataframe
    results.indel[i, "SNP_ID"] <- indel.df[snp.idx[1], "SNP_ID"]
    results.indel[i, "chrom"]  <- str_remove_all(string = indel.df[snp.idx[1], "chrom"], pattern = "[_maternalpaternal]")
    results.indel[i, "pos"]    <- indel.df[snp.idx[1], "SNP_ABS_POS"]
    results.indel[i, "strand"] <- indel.df[snp.idx[1], "strand"]
    
    results.indel[i,"ref.element.count"] <- length(ref)
    results.indel[i,"alt.element.count"] <- length(alt)
    
    results.indel[i,"K562.cDNA.ref.mean"] <- mean(indel.df[snp.idx[ref], "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.alt.mean"] <- mean(indel.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.ref.median"] <- median(indel.df[snp.idx[ref], "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.alt.median"] <- median(indel.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])
    
    
    results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(indel.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(indel.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.ref.median"] <- median(indel.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.alt.median"] <- median(indel.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])
    
    results.indel[i,"K562.wilcoxon.pvalue"] <- wilcox.test(indel.df[snp.idx[ref], "cDNA.K562.norm.ipcr"], indel.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])$p.value
    results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(indel.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"], indel.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])$p.value
    
    
  }
  
  foutput <- paste0(dir.input, k, ".output.RDS")
  saveRDS(results.indel, file = foutput)
  
} #End forEach









#####  STOP ABSOLUTELY HERE ######
# 
# 
# 
# 
# system.time(
# foreach(i = 1:length(snp.ids.vector)) %dopar% {
# 
#     snp.idx <- which(counts$SNP_ID == snp.ids.vector[i])
#     
#     ref <- which(counts[snp.idx, SNP_VAR] == 0)
#     alt <- which(counts[snp.idx, SNP_VAR] == 0)
# 
#     # Construct results datatable
#     results.indel[i, "SNP_ID"] <- counts[snp.idx[1],]$SNP_ID
#     results.indel[i, "chrom"]  <- str_remove_all(string = counts[snp.idx[1],]$chrom, pattern = "[_maternalpaternal]")
#     results.indel[i, "pos"]    <- counts[snp.idx[1],]$SNP_ABS_POS
#     results.indel[i, "strand"] <- counts[snp.idx[1],]$strand
# 
#     results.indel[i,"ref.element.count"] <- length(ref)
#     results.indel[i,"alt.element.count"] <- length(alt)
#   system.time({
#     ref.cDNA.k <- counts.df[snp.idx[ref],]$cDNA.K562.norm.ipcr
#     ref.cDNA.h <- counts[snp.idx[ref],]$cDNA.HepG2.norm.ipcr
#     alt.cDNA.k <- counts[snp.idx[alt],]$cDNA.K562.norm.ipcr
#     alt.cDNA.h <- counts[snp.idx[alt],]$cDNA.HepG2.norm.ipcr
#     
#     results.indel[i,"K562.cDNA.ref.mean"] <- mean(ref.cDNA.k)
#     results.indel[i,"K562.cDNA.alt.mean"] <- mean(alt.cDNA.k)
#     results.indel[i,"K562.cDNA.ref.median"] <- median(ref.cDNA.k)
#     results.indel[i,"K562.cDNA.alt.median"] <- median(alt.cDNA.k)
#     
#     results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(ref.cDNA.h)
#     results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(alt.cDNA.h)
#     results.indel[i,"HepG2.cDNA.ref.median"] <- median(ref.cDNA.h)
#     results.indel[i,"HepG2.cDNA.alt.median"] <- median(alt.cDNA.h)
#  })
#     results.indel[i,"K562.wilcoxon.pvalue"] <- wilcox.test(counts[snp.idx[ref], cDNA.K562.norm.ipcr], counts[snp.idx[alt], cDNA.K562.norm.ipcr])$p.value
#     results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(counts[snp.idx[ref], cDNA.HepG2.norm.ipcr], counts[snp.idx[alt], cDNA.HepG2.norm.ipcr])$p.value
# 
# }  )
# 
# 
# 
# 
# #### test some stuff
# 
# 
# foreach(i = snp.ids.vector, .combine = 'rbind') %dopar%{
#   
#   i
#   which(counts$SNP_ID ==  i)
#   
# }
# 
# 
# ## test2
# 
# indel.list <- foreach(i = 1:length(snp.ids.vector)) %dopar% {
#   
#   snp.idx <- which(counts$SNP_ID == snp.ids.vector[i])
# 
#   snp.df <- counts[snp.idx]
# }
# indel.list <- indel.list[2:3]
# 
# lapply(indel.list, function(x){
#  
#   snp.idx <- which(counts$SNP_ID == snp.ids.vector[i])
#   
#   x <- counts[snp.idx]
#   
#   ref <- which(x[, SNP_VAR] == 0)
#   alt <- which(x[, SNP_VAR] == 0)
#   
#   
#   
#   # Construct results datatable
#   results.indel[i, "SNP_ID"] <- x[1]$SNP_ID
#   results.indel[i, "chrom"]  <- str_remove_all(string = x[1]$chrom, pattern = "[_maternalpaternal]")
#   results.indel[i, "pos"]    <- x[1]$SNP_ABS_POS
#   results.indel[i, "strand"] <- x[1]$strand
#   
#   results.indel[i,"ref.element.count"] <- length(ref)
#   results.indel[i,"alt.element.count"] <- length(alt)
#   
#   results.indel[i,"K562.cDNA.ref.mean"] <- mean(x[ref,]$cDNA.K562.norm.ipcr)
#   results.indel[i,"K562.cDNA.alt.mean"] <- mean(x[alt,]$cDNA.K562.norm.ipcr)
#   results.indel[i,"K562.cDNA.ref.median"] <- median(x[ref,]$cDNA.K562.norm.ipcr)
#   results.indel[i,"K562.cDNA.alt.median"] <- median(x[alt,]$cDNA.K562.norm.ipcr)
#   
#   results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(x[ref,]$cDNA.HepG2.norm.ipcr)
#   results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(x[alt,]$cDNA.HepG2.norm.ipcr)
#   results.indel[i,"HepG2.cDNA.ref.median"] <- median(x[ref,]$cDNA.HepG2.norm.ipcr)
#   results.indel[i,"HepG2.cDNA.alt.median"] <- median(x[alt,]$cDNA.HepG2.norm.ipcr)
#   
#   
#   results.indel[i,"K562.wilcoxon.pvalue"] <-  wilcox.test(x[ref, cDNA.K562.norm.ipcr],  x[alt, cDNA.K562.norm.ipcr])$p.value
#   results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(x[ref, cDNA.HepG2.norm.ipcr], x[alt, cDNA.HepG2.norm.ipcr])$p.value
# })
