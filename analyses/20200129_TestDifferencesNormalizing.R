
# With this script I want to test the hypothesis that there is a difference 
# between normalizing with the total amount or reads or with the amount of reads
# after removing all reads not containing any variant. 

# define chr1 of SuRE 43-1 and then the equal 

file.name <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE43-1/equal/1.bedpe.gz"

# load read totals (ALL)

read.totals.file <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE43-1/read.totals.RDS"
sum.counts <-readRDS(read.totals.file)[,1]
sum.cdna1 <- readRDS(read.totals.file)[,2]
sum.cdna2 <- readRDS(read.totals.file)[,3]
sum.cdna3 <- readRDS(read.totals.file)[,4]
sum.cdna4 <- readRDS(read.totals.file)[,5]
sum.cdna5 <- readRDS(read.totals.file)[,6]

# load the data

sure.all <- fread(file.name, header = TRUE, sep = "\t", stringsAsFactors = FALSE, drop = c("BC", "start", "end", "SNP_ABS_POS", "SNP_SUBTYPE", "start_hg19", "end_hg19"))
sure.variants <- sure.all[sure.all$SNP_ID != ""]

# There is no need to split the rows, as it is about the read totals. 

# Normalize ALL

sure.all$count <- round(sure.all$count / sum(sure.all$count) * 1e9, digits = 2)
sure.all[,9]   <- round(sure.all[,9] / sum(sure.all[,9]) * 1e9, digits = 1)
sure.all[,10]  <- round(sure.all[,10] / sum(sure.all[,10]) * 1e9, digits = 1)
sure.all[,11]  <- round(sure.all[,11] / sum(sure.all[,11]) * 1e9, digits = 1)
sure.all[,12]  <- round(sure.all[,12] / sum(sure.all[,12]) * 1e9, digits = 1)
sure.all[,13]  <- round(sure.all[,13] / sum(sure.all[,13]) * 1e9, digits = 1)

# Normalize VARIANTS 

sure.variants$count <- round(sure.variants$count / sum(sure.variants$count) * 1e9, digits = 2)
sure.variants[,9]   <- round(sure.variants[,9] / sum(sure.variants[,9]) * 1e9, digits = 1)
sure.variants[,10]  <- round(sure.variants[,10] / sum(sure.variants[,10]) * 1e9, digits = 1)
sure.variants[,11]  <- round(sure.variants[,11] / sum(sure.variants[,10]) * 1e9, digits = 1)
sure.variants[,12]  <- round(sure.variants[,12] / sum(sure.variants[,10]) * 1e9, digits = 1)
sure.variants[,13]  <- round(sure.variants[,13] / sum(sure.variants[,10]) * 1e9, digits = 1)

# name columns VARIANTS

colnames(sure.variants)[9:13] <- c("cDNA.K562.B1", "cDNA.K562.B2", "cDNA.K562.B3", "cDNA.HepG2.B1", "cDNA.HepG2.B2")

# name columns ALL

colnames(sure.all)[9:13] <- c("cDNA.K562.B1", "cDNA.K562.B2", "cDNA.K562.B3", "cDNA.HepG2.B1", "cDNA.HepG2.B2")


# Generate SuRE signal VARIANTS

sure.variants$cDNA.K562.sum.norm <- (sure.variants$cDNA.K562.B1 + sure.variants$cDNA.K562.B2 + sure.variants$cDNA.K562.B3) / 3
sure.variants$cDNA.HepG2.sum.norm <- (sure.variants$cDNA.HepG2.B1 + sure.variants$cDNA.HepG2.B2) / 2

sure.variants$cDNA.K562.norm.ipcr <- sure.variants$cDNA.K562.sum.norm / sure.variants$count
sure.variants$cDNA.HepG2.norm.ipcr <- sure.variants$cDNA.HepG2.sum.norm / sure.variants$count


# Generate SuRE signal ALL

sure.all$cDNA.K562.sum.norm <- (sure.all$cDNA.K562.B1 + sure.all$cDNA.K562.B2 + sure.all$cDNA.K562.B3) / 3
sure.all$cDNA.HepG2.sum.norm <- (sure.all$cDNA.HepG2.B1 + sure.all$cDNA.HepG2.B2) / 2

sure.all$cDNA.K562.norm.ipcr <- sure.all$cDNA.K562.sum.norm / sure.all$count
sure.all$cDNA.HepG2.norm.ipcr <- sure.all$cDNA.HepG2.sum.norm / sure.all$count


# Generate plot

sure.all.filtered <- sure.all[sure.all$SNP_ID != ""]

sum(sure.all.filtered$SNP_ID == sure.variants$SNP_ID) == nrow(sure.variants)

sample.vector <- sample(seq(1,nrow(sure.variants)), size = 1000)
sample.all <- sure.all.filtered[sample.vector]
sample.variants <- sure.variants[sample.vector]

#plot(sure.all.filtered$cDNA.K562.norm.ipcr, sure.variants$cDNA.K562.norm.ipcr)
plot(sample.all$cDNA.K562.norm.ipcr, sample.variants$cDNA.K562.norm.ipcr)
abline(a = 0, b = 1)

# test to see what the differences are
test <- cbind(sure.variants[sample.vector], sure.all.filtered[sample.vector])
