### TO DO:
### ADD RANDOMIZED P VALUES??
### change the dataset to subset from (`counts`) to a much smaller dataset containing only the usefull variants. Use the `snps` variable
### make sure everything is dataframe
### variant annotation is not yet correct as the multiple alleles lead to differen granges. 
### plus and minus strand is NOT GOOD ANNOTATED

# Script to calculate p-values between alternative and reference sequence

# Load the required libraries
library(data.table)
library(tidyverse)


# Load the required dataset

counts <- sure.indel.allchrom
counts.df <- as.data.frame(counts)
#OR
#fname <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE42-1/combined/sure.42-1.counts.indel.all.RDS"
#counts <- load(file = fname)


# Remove Rows with SNP_ID == "."
counts.df <- counts.df[counts.df$SNP_ID != ".",]

# Get the indel.ids of the indels that have reads for reference and alternative alleles and
# construct a new dataframe with only these indels. 

indel.ids.vector <- names(which(tapply(counts.df$SNP_VAR, counts.df$SNP_ID, function(x) {all(c(0,1) %in% x)} ) == TRUE))
indel.df <- counts.df[counts.df$SNP_ID %in% indel.ids.vector,]



# Generate an empty results dataframe for the indels 

results.indel <- data.frame(matrix(nrow = length(indel.ids.vector), ncol = 17))
colnames(results.indel) <- c("SNP_ID", "chrom", "pos", "strand", "ref.element.count", "alt.element.count", "K562.cDNA.ref.mean","K562.cDNA.alt.mean", "K562.cDNA.ref.median", "K562.cDNA.alt.median", "HepG2.cDNA.ref.mean","HepG2.cDNA.alt.mean","HepG2.cDNA.ref.median", "HepG2.cDNA.alt.median", "K562.wilcoxon.pvalue", "HepG2.wilcoxon.pvalue", "location.annotation")

##############################################
# Prepare variant (e.g. promoter) annotation #
##############################################  

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantTools)
library(dplyr)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Generate a granges file from the results.sure dataframe
snp.granges <- reduce(GRanges(seqnames = str_remove_all(string = indel.df[, "chrom"], pattern = "[_maternalpaternal]"),
                              ranges = IRanges(start = as.numeric(indel.df[,"SNP_ABS_POS_hg19"]),
                                               end =   as.numeric(indel.df[,"SNP_ABS_POS_hg19"]))
                              ))

# Locate variants
variants.granges <- locateVariants(snp.granges, txdb, AllVariants())

# Save variants in a dataframe
df.chr <- seqnames(variants.granges)
df.snp.abs.pos <- start(variants.granges)
df.snp.pos <- paste0(df.chr,":", df.snp.abs.pos)
df.variant <- variants.granges$LOCATION
variants.df <- data.frame(df.chr, df.snp.abs.pos, df.snp.pos, df.variant)

# Use the distinct function from dplyr to remove identical rows
# These could for example contain multiple rows that contain an 
# intron

variants.df.nonredundant <- distinct(variants.df)

############################################
# Start actual loop through all the Indels #
############################################


system.time(
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
)































## Option 2 ## 
system.time({

  for (i in 1:length(indel.ids.vector)){
    
    
    print(i)
    snp.idx <- which(counts.df$SNP_ID == indel.ids.vector[i])
    
    snp.df <- counts.df[snp.idx, ]
    
    ref <- which(snp.df[, "SNP_VAR"] == 0)
    alt <- which(snp.df[, "SNP_VAR"] == 1)
    
    
    
    # Construct results datatable
    results.indel[i, "SNP_ID"] <- snp.df[1,"SNP_ID"]
    results.indel[i, "chrom"]  <- str_remove_all(string = snp.df[1,"chrom"] , pattern = "[_maternalpaternal]")
    results.indel[i, "pos"]    <- snp.df[1,"SNP_ABS_POS"]
    results.indel[i, "strand"] <- snp.df[1,"strand"]
    
    results.indel[i,"ref.element.count"] <- length(ref)
    results.indel[i,"alt.element.count"] <- length(alt)
   
    results.indel[i,"K562.cDNA.ref.mean"] <- mean(snp.df[ref, "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.alt.mean"] <- mean(snp.df[alt, "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.ref.median"] <- median(snp.df[ref, "cDNA.K562.norm.ipcr"])
    results.indel[i,"K562.cDNA.alt.median"] <- median(snp.df[alt, "cDNA.K562.norm.ipcr"])
    
    results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(snp.df[ref, "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(snp.df[alt, "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.ref.median"] <- median(snp.df[ref, "cDNA.HepG2.norm.ipcr"])
    results.indel[i,"HepG2.cDNA.alt.median"] <- median(snp.df[alt, "cDNA.HepG2.norm.ipcr"])
   
    
    results.indel[i,"K562.wilcoxon.pvalue"] <-  wilcox.test(snp.df[ref, "cDNA.K562.norm.ipcr"],  snp.df[alt, "cDNA.K562.norm.ipcr"])$p.value
    results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(snp.df[ref, "cDNA.HepG2.norm.ipcr"], snp.df[alt, "cDNA.HepG2.norm.ipcr"])$p.value

  }
  
})




