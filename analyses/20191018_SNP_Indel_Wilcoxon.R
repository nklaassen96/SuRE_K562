### TO DO:
### variant annotation is not yet correct as the multiple alleles lead to differen granges. 
### plus and minus strand is NOT GOOD ANNOTATED

# Script to calculate p-values between alternative and reference sequence

# Load the required libraries
library(data.table)
library(tidyverse)
library(foreach)
library(doMC)
registerDoMC(cores = 6)

# Loop through all files (1 file per chromosome)

dir.input  = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE_combined/sure.snp_indel.combined.chrom."
dir.output = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_pvalue/"

print(paste("start function at", Sys.time()))



foreach(k=c(1:22,"X")) %dopar% {
  
  # open the file and remove unnecesary columns
  
  file.name <- paste0(dir.input, k, ".RDS")  
  counts.df <- readRDS(file = file.name)
  counts.df[,c("count", "SNP_PARENT", "cDNA.K562.B1", "cDNA.K562.B2", "cDNA.K562.B3", "cDNA.HepG2.B1", "cDNA.HepG2.B2", "library", "cDNA.K562.sum.norm", "cDNA.HepG2.sum.norm")] <- NULL
  gc()
  
  

  # 1. order the elements based on snp.id and remake the datatable to a dataframe because the function is much quicker that way
  counts.df <- as.data.frame(counts.df[order(counts.df$SNP_ID),])
  
  
  gc()
  
  # 2. locate the position where a new snp starts by finding duplicate snp ids. Locate the 
  # first and last snp position to generate blocks
  
  new.snp <- which(!duplicated(counts.df$SNP_ID))
  new.snp.position <- new.snp[seq(from=1, to = length(new.snp), by = 1000)]
  last.snp.position <- c(new.snp.position-1, nrow(counts.df))

  
  
  # loop through the blocks of snps that are generated in this loop
  
  all.results.indel <- NULL
  
  for (b in c(1:length(new.snp.position))){
    
  print(paste0(Sys.time(), " start block ",b,"/",length(new.snp.position), " chromosome ", k))
  
  reduced.counts.df <- counts.df[new.snp.position[b]:last.snp.position[b+1],]
    
  
  
   
  
  
  # Get the indel.ids of the indels that have reads for reference and alternative alleles and
  # construct a new dataframe with only these indels. 
  
  snp_indel.ids.vector <- names(which(tapply(reduced.counts.df$SNP_VAR, reduced.counts.df$SNP_ID, function(x) {all(c(0,1) %in% x)} ) == TRUE))
  print(length(snp_indel.ids.vector))
  
  snp_indel.df <- reduced.counts.df[reduced.counts.df$SNP_ID %in% snp_indel.ids.vector,]
  
  
  
  # Generate an empty results dataframe for the indels 
  
  results.indel <- data.frame(matrix(nrow = length(snp_indel.ids.vector), ncol = 20))
  colnames(results.indel) <- c("SNP_ID","snp.type", "ref.seq", "alt.seq", "chrom", "pos.hg19", "ref.element.count", "alt.element.count", "K562.cDNA.ref.mean","K562.cDNA.alt.mean", "K562.cDNA.ref.median", "K562.cDNA.alt.median", "HepG2.cDNA.ref.mean","HepG2.cDNA.alt.mean","HepG2.cDNA.ref.median", "HepG2.cDNA.alt.median", "K562.wilcoxon.pvalue", "HepG2.wilcoxon.pvalue", "K562.wilcoxon.pvalue.random", "HepG2.wilcoxon.pvalue.random")
  

  
  ############################################
  # Start actual loop through all the Indels and SNPS #
  ############################################
  
  
 
  for (i in 1:length(snp_indel.ids.vector)){
    #if (i %% 250 == 0){print(i)}
    
    snp.idx <- which(snp_indel.df$SNP_ID == snp_indel.ids.vector[i])
    
    ref <- which(snp_indel.df[snp.idx, "SNP_VAR"] == 0)
    alt <- which(snp_indel.df[snp.idx, "SNP_VAR"] == 1)
    
    ref.random <- sample(c(ref,alt), size = length(ref))
    alt.random <- c(ref,alt)[!c(ref,alt)%in%ref.random]
    
    # Construct results dataframe
    results.indel[i,"SNP_ID"]  <- snp_indel.df[snp.idx[1], "SNP_ID"]
    results.indel[i,"chrom"]   <- str_remove_all(string = snp_indel.df[snp.idx[1], "chrom"], pattern = "[_maternalpaternal]")
    results.indel[i,"pos.hg19"]     <- snp_indel.df[snp.idx[1], "SNP_ABS_POS_hg19"]
    results.indel[i,"snp.type"]<- snp_indel.df[snp.idx[1], "SNP_TYPE"]
    
    results.indel[i,"ref.seq"] <- snp_indel.df[snp.idx[ref[1]], "SNP_SEQ"]
    results.indel[i,"alt.seq"] <- snp_indel.df[snp.idx[alt[1]], "SNP_SEQ"]
    
    results.indel[i,"ref.element.count"] <- length(ref)
    results.indel[i,"alt.element.count"] <- length(alt)
    
    results.indel[i,"K562.cDNA.ref.mean"] <- mean(snp_indel.df[snp.idx[ref], "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.alt.mean"] <- mean(snp_indel.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.ref.median"] <- median(snp_indel.df[snp.idx[ref], "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.alt.median"] <- median(snp_indel.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])
    
    
    results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(snp_indel.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(snp_indel.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.ref.median"] <- median(snp_indel.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.alt.median"] <- median(snp_indel.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])
    
    results.indel[i,"K562.wilcoxon.pvalue"] <-  wilcox.test(snp_indel.df[snp.idx[ref], "cDNA.K562.norm.ipcr"], snp_indel.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])$p.value
    results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(snp_indel.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"], snp_indel.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])$p.value
    results.indel[i,"K562.wilcoxon.pvalue.random"] <- wilcox.test(snp_indel.df[snp.idx[ref.random], "cDNA.K562.norm.ipcr"], snp_indel.df[snp.idx[alt.random], "cDNA.K562.norm.ipcr"])$p.value
    results.indel[i,"HepG2.wilcoxon.pvalue.random"] <- wilcox.test(snp_indel.df[snp.idx[ref.random], "cDNA.HepG2.norm.ipcr"], snp_indel.df[snp.idx[alt.random], "cDNA.HepG2.norm.ipcr"])$p.value
    
  
  } # end for-loop through 1000 snps
  
  all.results.indel <- rbind(all.results.indel, results.indel)
  
  } # end for-loop through snp-blocks
  
  file.name <- paste0(dir.output, "sure.snp_indel.dataframe.pvalue.chrom.", k,".RDS")
  saveRDS(object = all.results.indel, file = file.name)
  print(paste("end function chrom",k, Sys.time()))
  
} # end foreach-loop through files


#After generation of all the individual files, load them all and concatonate them
print(paste( Sys.time(), "start concatonating"))
results.indel.all <- NULL

for (i in c(1:22,"X")){
  print(paste(Sys.time(), "start", i))
  # 1. Define data
  dir.input <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_pvalue/sure.snp_indel.dataframe.pvalue.chrom."
  file.name <- paste0(dir.input, i, ".RDS")
  0
  
  # 2. Load the data
  sure.df <- readRDS(file = file.name)
  
  
  # 3. Concatonate the data
  results.indel.all <- rbind(results.indel.all, sure.df)
  
}

# Add some extra columns

results.indel.all$max.elements <- apply(results.indel.all[,c("ref.element.count","alt.element.count")], 1, max)
results.indel.all$min.elements <- apply(results.indel.all[,c("ref.element.count","alt.element.count")], 1, min)
results.indel.all$max.k562.expression <- apply(results.indel.all[,c("K562.cDNA.alt.mean", "K562.cDNA.ref.mean")], 1, max)
results.indel.all$max.hepg2.expression <- apply(results.indel.all[,c("HepG2.cDNA.alt.mean", "HepG2.cDNA.ref.mean")], 1, max)

results.indel.all$k562.max <- c("ref", "alt")[apply(results.indel.all[,c("K562.cDNA.ref.mean","K562.cDNA.alt.mean")], 1, which.max)]
results.indel.all$hepg2.max <-c("ref", "alt")[apply(results.indel.all[,c("HepG2.cDNA.ref.mean","HepG2.cDNA.alt.mean")], 1, which.max)]

# Save the file

file.name <- paste0(dir.output, "sure.snp_indel.dataframe.pvalue.all.RDS")
saveRDS(object = results.indel.all, file = file.name)

print(paste(Sys.time(), "finish concatonating"))











