
dir.input.com <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE42-1/combined/"
fls <- list.files(path = dir.input.com)

counts.all.chrom <- NULL

for (counts.per.chrom in fls){
  
  print(paste0("start ", counts.per.chrom, " at ", Sys.time()))
  
  fdir <- paste0(dir.input.com, counts.per.chrom)
  
  counts.chrom <- readRDS(fdir)
  
  counts.all.chrom <- rbind(counts.all.chrom, counts.chrom)
  
}

fname <- paste0(dir.input.com, "counts.all.chrom.RDS")

saveRDS(counts.all.chrom, file = fname)

counts.all.chrom.backup <- counts.all.chrom
colnames(counts.all.chrom[,5])


#reads per billion = read / sumreads * 1 billion

counts.all.chrom$count <- round(counts.all.chrom$count / sum(counts.all.chrom$count) * 1e9, digits = 0)
counts.all.chrom[,14]  <- round(counts.all.chrom[,14] / sum(counts.all.chrom[,14]) * 1e9, digits = 0)
counts.all.chrom[,15]  <- round(counts.all.chrom[,15] / sum(counts.all.chrom[,15]) * 1e9, digits = 0)
counts.all.chrom[,16]  <- round(counts.all.chrom[,16] / sum(counts.all.chrom[,16]) * 1e9, digits = 0)
counts.all.chrom[,17]  <- round(counts.all.chrom[,17] / sum(counts.all.chrom[,17]) * 1e9, digits = 0)
counts.all.chrom[,18]  <- round(counts.all.chrom[,18] / sum(counts.all.chrom[,18]) * 1e9, digits = 0)

### test om te kijken of t goed gaat met die colom namen

tst32 <- tail(counts.all.chrom)
tst23 <- fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE42-2/equal/1.bedpe.gz", nrows = 100)
tst23 <- tst23[75:79,]

colnames(tst32)[14:18] <- c("cDNA.K562.B1", "cDNA.K562.B2", "cDNA.K562.B3", "cDNA.HepG2.B1", "cDNA.HepG2.B2")
colnames(tst23)[14:18] <- c("cDNA.HepG2.B1", "cDNA.HepG2.B2", "cDNA.K562.B1", "cDNA.K562.B2", "cDNA.K562.B3")



















######## TEST WHAT WENT WRONG WITH THE STR.SPLIT FUNCTION. CONCLUSION: NOTHING WENT WRONG ####

'
list.abspos <- strsplit(sure.BCcounts.tst$SNP_ABS_POS, split = ",")
list.seq <- strsplit(sure.BCcounts.tst$SNP_SEQ, split = ",")

strings.q <- sapply(list.abspos, length)
strings.s <- sapply(list.seq, length)

#strings.q <- str_count(sure.BCcounts.tst$SNP_ABS_POS, pattern = ",")
#strings.s <- str_count(sure.BCcounts.tst$SNP_SEQ, pattern = ",")
which(strings.q != strings.s)
'


