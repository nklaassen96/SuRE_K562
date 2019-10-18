### TO DO:
### variant annotation is not yet correct as the multiple alleles lead to differen granges. 
### plus and minus strand is NOT GOOD ANNOTATED

# Script to calculate p-values between alternative and reference sequence

# Load the required libraries
library(data.table)
library(tidyverse)
library(foreach)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantTools)
library(dplyr)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(doMC)
registerDoMC(cores = 23)

# Loop through all files (1 file per chromosome)

dir.input  = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE_combined/sure.indel.combined.chrom."
dir.output = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_pvalue/"

print(paste("start function at", Sys.time()))



foreach(k=c(1:22,"X")) %dopar% {
#for (k in c(1:22,"X")) {
    
file.name <- paste0(dir.input, k, ".RDS")  
counts.df <- as.data.frame(readRDS(file = file.name))
  



# Get the indel.ids of the indels that have reads for reference and alternative alleles and
# construct a new dataframe with only these indels. 

indel.ids.vector <- names(which(tapply(counts.df$SNP_VAR, counts.df$SNP_ID, function(x) {all(c(0,1) %in% x)} ) == TRUE))
indel.df <- counts.df[counts.df$SNP_ID %in% indel.ids.vector,]



# Generate an empty results dataframe for the indels 

results.indel <- data.frame(matrix(nrow = length(indel.ids.vector), ncol = 18))
colnames(results.indel) <- c("SNP_ID", "chrom", "pos", "ref.element.count", "alt.element.count", "K562.cDNA.ref.mean","K562.cDNA.alt.mean", "K562.cDNA.ref.median", "K562.cDNA.alt.median", "HepG2.cDNA.ref.mean","HepG2.cDNA.alt.mean","HepG2.cDNA.ref.median", "HepG2.cDNA.alt.median", "K562.wilcoxon.pvalue", "HepG2.wilcoxon.pvalue", "K562.wilcoxon.pvalue.random", "HepG2.wilcoxon.pvalue.random", "location.annotation")

##############################################
# Prepare variant (e.g. promoter) annotation #
##############################################  

# # Generate a granges file from the results.sure dataframe
# # The reduce function is used to remove similar regions
# snp.granges <- reduce(GRanges(seqnames = paste0("chr", str_remove_all(string = indel.df[, "chrom"], pattern = "[_maternalpaternal]")),
#                               ranges = IRanges(start = as.numeric(indel.df[,"SNP_ABS_POS_hg19"]),
#                                                end =   as.numeric(indel.df[,"SNP_ABS_POS_hg19"]))
#                               ))
# 
# # Locate variants
# variants.granges <- locateVariants(snp.granges, txdb, AllVariants())
# 
# # Save variants in a dataframe
# df.chr <- seqnames(variants.granges)
# df.snp.abs.pos <- start(variants.granges)
# df.snp.pos <- paste0(df.chr,":", df.snp.abs.pos)
# df.variant <- variants.granges$LOCATION
# variants.df <- data.frame(df.chr, df.snp.abs.pos, df.snp.pos, df.variant, stringsAsFactors = FALSE)
# 
# # Use the distinct function from dplyr to remove identical rows
# # These could for example contain multiple rows that contain an
# # intron for the same location due to multiple transcripts
# 
# variants.df.nonredundant <- distinct(variants.df)

############################################
# Start actual loop through all the Indels #
############################################



for (i in 1:length(indel.ids.vector)){
  if (i %% 50 == 0){print(i)}
  snp.idx <- which(indel.df$SNP_ID == indel.ids.vector[i])
  
  ref <- which(indel.df[snp.idx, "SNP_VAR"] == 0)
  alt <- which(indel.df[snp.idx, "SNP_VAR"] == 1)
  
  ref.random <- sample(c(ref,alt), size = length(ref))
  alt.random <- c(ref,alt)[!c(ref,alt)%in%ref.random]
  
  # Construct results dataframe
  results.indel[i, "SNP_ID"] <- indel.df[snp.idx[1], "SNP_ID"]
  results.indel[i, "chrom"]  <- str_remove_all(string = indel.df[snp.idx[1], "chrom"], pattern = "[_maternalpaternal]")
  results.indel[i, "pos"]    <- indel.df[snp.idx[1], "SNP_ABS_POS_hg19"]
  
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
  
  results.indel[i,"K562.wilcoxon.pvalue"] <-  wilcox.test(indel.df[snp.idx[ref], "cDNA.K562.norm.ipcr"], indel.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])$p.value
  results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(indel.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"], indel.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])$p.value
  results.indel[i,"K562.wilcoxon.pvalue.random"] <- wilcox.test(indel.df[snp.idx[ref.random], "cDNA.K562.norm.ipcr"], indel.df[snp.idx[alt.random], "cDNA.K562.norm.ipcr"])$p.value
  results.indel[i, "HepG2.wilcoxon.pvalue.random"] <- wilcox.test(indel.df[snp.idx[ref.random], "cDNA.HepG2.norm.ipcr"], indel.df[snp.idx[alt.random], "cDNA.HepG2.norm.ipcr"])$p.value
  
  # # Generate a unique snpi.id in the format "chr1:3824989"
  # snp.id <- paste0("chr", results.indel[i,"chrom"], ":", results.indel[i, "pos"])
  # 
  # 
  # # Generate a factor stating all variants for that specific snp (could be >1)
  # # Then sort the factor based on the levels and take the first (most important)
  # # variant as which we want it to classify
  # 
  # 
  # variants <- variants.df.nonredundant[variants.df.nonredundant$df.snp.pos == snp.id, "df.variant"]
  # variants.sorted <- sort(factor(variants, levels = c("promoter","spliceSite", "coding", "fiveUTR", "threeUTR","intron", "intergenic")))
  # 
  # 
  # results.indel[i,"location.annotation"] <- as.character(variants.sorted[1])
  
} # end for-loop through snps

file.name <- paste0(dir.output, "sure.indel.dataframe.pvalue.chrom.", k,".RDS")
saveRDS(object = results.indel, file = file.name)
print(paste("end function chrom",k, Sys.time()))

} # end foreach-loop through files


#After generation of all the individual files, load them all and concatonate them
print(paste("start concatonating", Sys.time()))
results.indel.all <- NULL

for (i in c(1:22,"X")){
  print(paste(Sys.time(), "start", i))
  # 1. Define data
  dir.input <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_pvalue/sure.indel.dataframe.pvalue.chrom."
  file.name <- paste0(dir.input, i, ".RDS")
    
  
  # 2. Load the data
  sure.df <- readRDS(file = file.name)
  
  
  # 3. Concatonate the data
  results.indel.all <- rbind(results.indel.all, sure.df)
  
}
file.name <- paste0(dir.output, "sure.indel.dataframe.pvalue.all.RDS")
saveRDS(object = results.indel.all, file = file.name)

print(paste("finish concatonating", Sys.time()))


















# 
# 
# ## Option 2 ## 
# system.time({
# 
#   for (i in 1:length(indel.ids.vector)){
#     
#     
#     print(i)
#     snp.idx <- which(counts.df$SNP_ID == indel.ids.vector[i])
#     
#     snp.df <- counts.df[snp.idx, ]
#     
#     ref <- which(snp.df[, "SNP_VAR"] == 0)
#     alt <- which(snp.df[, "SNP_VAR"] == 1)
#     
#     
#     
#     # Construct results datatable
#     results.indel[i, "SNP_ID"] <- snp.df[1,"SNP_ID"]
#     results.indel[i, "chrom"]  <- str_remove_all(string = snp.df[1,"chrom"] , pattern = "[_maternalpaternal]")
#     results.indel[i, "pos"]    <- snp.df[1,"SNP_ABS_POS"]
#     results.indel[i, "strand"] <- snp.df[1,"strand"]
#     
#     results.indel[i,"ref.element.count"] <- length(ref)
#     results.indel[i,"alt.element.count"] <- length(alt)
#    
#     results.indel[i,"K562.cDNA.ref.mean"] <- mean(snp.df[ref, "cDNA.K562.norm.ipcr"])
#     results.indel[i,"K562.cDNA.alt.mean"] <- mean(snp.df[alt, "cDNA.K562.norm.ipcr"])
#     results.indel[i,"K562.cDNA.ref.median"] <- median(snp.df[ref, "cDNA.K562.norm.ipcr"])
#     results.indel[i,"K562.cDNA.alt.median"] <- median(snp.df[alt, "cDNA.K562.norm.ipcr"])
#     
#     results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(snp.df[ref, "cDNA.HepG2.norm.ipcr"])
#     results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(snp.df[alt, "cDNA.HepG2.norm.ipcr"])
#     results.indel[i,"HepG2.cDNA.ref.median"] <- median(snp.df[ref, "cDNA.HepG2.norm.ipcr"])
#     results.indel[i,"HepG2.cDNA.alt.median"] <- median(snp.df[alt, "cDNA.HepG2.norm.ipcr"])
#    
#     
#     results.indel[i,"K562.wilcoxon.pvalue"] <-  wilcox.test(snp.df[ref, "cDNA.K562.norm.ipcr"],  snp.df[alt, "cDNA.K562.norm.ipcr"])$p.value
#     results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(snp.df[ref, "cDNA.HepG2.norm.ipcr"], snp.df[alt, "cDNA.HepG2.norm.ipcr"])$p.value
# 
#   }
#   
# })
# 
# 
# 
# 
