# Test the foreach loop

install.packages("foreach")
library(foreach)

install.packages("doMC")
library(doMC)
registerDoMC(cores = 20)
system.time(
foreach(i = 1:length(snp.ids.vector)) %dopar% {

    snp.idx <- which(counts$SNP_ID == snp.ids.vector[i])
    
    ref <- which(counts[snp.idx, SNP_VAR] == 0)
    alt <- which(counts[snp.idx, SNP_VAR] == 0)

    # Construct results datatable
    results.indel[i, "SNP_ID"] <- counts[snp.idx[1],]$SNP_ID
    results.indel[i, "chrom"]  <- str_remove_all(string = counts[snp.idx[1],]$chrom, pattern = "[_maternalpaternal]")
    results.indel[i, "pos"]    <- counts[snp.idx[1],]$SNP_ABS_POS
    results.indel[i, "strand"] <- counts[snp.idx[1],]$strand

    results.indel[i,"ref.element.count"] <- length(ref)
    results.indel[i,"alt.element.count"] <- length(alt)
  system.time({
    ref.cDNA.k <- counts.df[snp.idx[ref],]$cDNA.K562.norm.ipcr
    ref.cDNA.h <- counts[snp.idx[ref],]$cDNA.HepG2.norm.ipcr
    alt.cDNA.k <- counts[snp.idx[alt],]$cDNA.K562.norm.ipcr
    alt.cDNA.h <- counts[snp.idx[alt],]$cDNA.HepG2.norm.ipcr
    
    results.indel[i,"K562.cDNA.ref.mean"] <- mean(ref.cDNA.k)
    results.indel[i,"K562.cDNA.alt.mean"] <- mean(alt.cDNA.k)
    results.indel[i,"K562.cDNA.ref.median"] <- median(ref.cDNA.k)
    results.indel[i,"K562.cDNA.alt.median"] <- median(alt.cDNA.k)
    
    results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(ref.cDNA.h)
    results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(alt.cDNA.h)
    results.indel[i,"HepG2.cDNA.ref.median"] <- median(ref.cDNA.h)
    results.indel[i,"HepG2.cDNA.alt.median"] <- median(alt.cDNA.h)
 })
    results.indel[i,"K562.wilcoxon.pvalue"] <- wilcox.test(counts[snp.idx[ref], cDNA.K562.norm.ipcr], counts[snp.idx[alt], cDNA.K562.norm.ipcr])$p.value
    results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(counts[snp.idx[ref], cDNA.HepG2.norm.ipcr], counts[snp.idx[alt], cDNA.HepG2.norm.ipcr])$p.value

}  )




#### test some stuff


foreach(i = snp.ids.vector, .combine = 'rbind') %dopar%{
  
  i
  which(counts$SNP_ID ==  i)
  
}


## test2

indel.list <- foreach(i = 1:length(snp.ids.vector)) %dopar% {
  
  snp.idx <- which(counts$SNP_ID == snp.ids.vector[i])

  snp.df <- counts[snp.idx]
}
indel.list <- indel.list[2:3]

lapply(indel.list, function(x){
 
  snp.idx <- which(counts$SNP_ID == snp.ids.vector[i])
  
  x <- counts[snp.idx]
  
  ref <- which(x[, SNP_VAR] == 0)
  alt <- which(x[, SNP_VAR] == 0)
  
  
  
  # Construct results datatable
  results.indel[i, "SNP_ID"] <- x[1]$SNP_ID
  results.indel[i, "chrom"]  <- str_remove_all(string = x[1]$chrom, pattern = "[_maternalpaternal]")
  results.indel[i, "pos"]    <- x[1]$SNP_ABS_POS
  results.indel[i, "strand"] <- x[1]$strand
  
  results.indel[i,"ref.element.count"] <- length(ref)
  results.indel[i,"alt.element.count"] <- length(alt)
  
  results.indel[i,"K562.cDNA.ref.mean"] <- mean(x[ref,]$cDNA.K562.norm.ipcr)
  results.indel[i,"K562.cDNA.alt.mean"] <- mean(x[alt,]$cDNA.K562.norm.ipcr)
  results.indel[i,"K562.cDNA.ref.median"] <- median(x[ref,]$cDNA.K562.norm.ipcr)
  results.indel[i,"K562.cDNA.alt.median"] <- median(x[alt,]$cDNA.K562.norm.ipcr)
  
  results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(x[ref,]$cDNA.HepG2.norm.ipcr)
  results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(x[alt,]$cDNA.HepG2.norm.ipcr)
  results.indel[i,"HepG2.cDNA.ref.median"] <- median(x[ref,]$cDNA.HepG2.norm.ipcr)
  results.indel[i,"HepG2.cDNA.alt.median"] <- median(x[alt,]$cDNA.HepG2.norm.ipcr)
  
  
  results.indel[i,"K562.wilcoxon.pvalue"] <-  wilcox.test(x[ref, cDNA.K562.norm.ipcr],  x[alt, cDNA.K562.norm.ipcr])$p.value
  results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(x[ref, cDNA.HepG2.norm.ipcr], x[alt, cDNA.HepG2.norm.ipcr])$p.value
})
