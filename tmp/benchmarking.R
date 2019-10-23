#tst for speeds

library(microbenchmark)

# 1 BEST ; dont use others (at leastt 800% slower)
microbenchmark(mean(counts.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"]),
               mean(counts.df[snp.idx[ref], ]$cDNA.HepG2.norm.ipcr),
               mean(counts[snp.idx[ref], `cDNA.HepG2.norm.ipcr`]),
               mean(counts[snp.idx[ref], ]$cDNA.HepG2.norm.ipcr),
               
               times = 100
               )

# 1,2,3 = GOOD ; # 4 = BAD (40% slower)
microbenchmark(which(counts.df[,"SNP_ID"] == snp.ids.vector[1]),
               which(counts.df$SNP_ID == snp.ids.vector[1]),
               which(counts$SNP_ID == snp.ids.vector[1]),
               which(counts[,`SNP_ID`] == snp.ids.vector[1]),
               
               times = 10)

microbenchmark({
  
  # 1. dataframe
  system.time({
  snp.idx <- which(counts.df$SNP_ID == snp.ids.vector[i])
  
  ref <- which(counts.df[snp.idx, "SNP_VAR"] == 0)
  alt <- which(counts.df[snp.idx, "SNP_VAR"] == 0)
  
  # Construct results datatable
  results.indel[i, "SNP_ID"] <- counts.df[snp.idx[1], "SNP_ID"]
  results.indel[i, "chrom"]  <- str_remove_all(string = counts.df[snp.idx[1], "chrom"], pattern = "[_maternalpaternal]")
  results.indel[i, "pos"]    <- counts.df[snp.idx[1], "SNP_ABS_POS"]
  results.indel[i, "strand"] <- counts.df[snp.idx[1], "strand"]
  
  results.indel[i,"ref.element.count"] <- length(ref)
  results.indel[i,"alt.element.count"] <- length(alt)
  
  results.indel[i,"K562.cDNA.ref.mean"] <- mean(counts.df[snp.idx[ref], "cDNA.K562.norm.ipcr"])
  results.indel[i,"K562.cDNA.alt.mean"] <- mean(counts.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])
  results.indel[i,"K562.cDNA.ref.median"] <- median(counts.df[snp.idx[ref], "cDNA.K562.norm.ipcr"])
  results.indel[i,"K562.cDNA.alt.median"] <- median(counts.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])
  
  
  results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(counts.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"])
  results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(counts.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])
  results.indel[i,"HepG2.cDNA.ref.median"] <- median(counts.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"])
  results.indel[i,"HepG2.cDNA.alt.median"] <- median(counts.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])
  
  results.indel[i,"K562.wilcoxon.pvalue"] <- wilcox.test(counts.df[snp.idx[ref], "cDNA.K562.norm.ipcr"], counts.df[snp.idx[alt], "cDNA.K562.norm.ipcr"])$p.value
  results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(counts.df[snp.idx[ref], "cDNA.HepG2.norm.ipcr"], counts.df[snp.idx[alt], "cDNA.HepG2.norm.ipcr"])$p.value
  })
},{
  system.time({
  # 2. datatable
  snp.idx <- which(counts$SNP_ID == snp.ids.vector[i])
  
  ref <- which(counts[snp.idx, `SNP_VAR`] == 0)
  alt <- which(counts[snp.idx, `SNP_VAR`] == 0)
  
  # Construct results datatable
  results.indel[i, "SNP_ID"] <- counts[snp.idx[1], `SNP_ID`]
  results.indel[i, "chrom"]  <- str_remove_all(string = counts[snp.idx[1], `chrom`], pattern = "[_maternalpaternal]")
  results.indel[i, "pos"]    <- counts[snp.idx[1], `SNP_ABS_POS`]
  results.indel[i, "strand"] <- counts[snp.idx[1], `strand`]
  
  results.indel[i,"ref.element.count"] <- length(ref)
  results.indel[i,"alt.element.count"] <- length(alt)
  
  results.indel[i,"K562.cDNA.ref.mean"] <- mean(counts[snp.idx[ref], `cDNA.K562.norm.ipcr`])
  results.indel[i,"K562.cDNA.alt.mean"] <- mean(counts[snp.idx[alt], `cDNA.K562.norm.ipcr`])
  results.indel[i,"K562.cDNA.ref.median"] <- median(counts[snp.idx[ref], `cDNA.K562.norm.ipcr`])
  results.indel[i,"K562.cDNA.alt.median"] <- median(counts[snp.idx[alt], `cDNA.K562.norm.ipcr`])
  
  
  results.indel[i,"HepG2.cDNA.ref.mean"] <- mean(counts[snp.idx[ref], `cDNA.HepG2.norm.ipcr`])
  results.indel[i,"HepG2.cDNA.alt.mean"] <- mean(counts[snp.idx[alt], `cDNA.HepG2.norm.ipcr`])
  results.indel[i,"HepG2.cDNA.ref.median"] <- median(counts[snp.idx[ref], `cDNA.HepG2.norm.ipcr`])
  results.indel[i,"HepG2.cDNA.alt.median"] <- median(counts[snp.idx[alt], `cDNA.HepG2.norm.ipcr`])
  
  results.indel[i,"K562.wilcoxon.pvalue"] <- wilcox.test(counts[snp.idx[ref], `cDNA.K562.norm.ipcr`], counts[snp.idx[alt], `cDNA.K562.norm.ipcr`])$p.value
  results.indel[i,"HepG2.wilcoxon.pvalue"] <- wilcox.test(counts[snp.idx[ref], `cDNA.HepG2.norm.ipcr`], counts[snp.idx[alt], `cDNA.HepG2.norm.ipcr`])$p.value
})
  
  })
