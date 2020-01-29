# test to see what percentage of fragments contain an indel

library(data.table)
library(doMC)
registerDoMC(cores = 10)

# sure.22.eq <- fread("data/interim/SuRE_Indels_gDNA_Count/SuRE43-1/equal/X.bedpe.gz")
# sure.22.p <-  fread("data/interim/SuRE_Indels_gDNA_Count/SuRE43-1/paternal/X.bedpe.gz")
# sure.22.m <-  fread("data/interim/SuRE_Indels_gDNA_Count/SuRE43-1/maternal/X.bedpe.gz")
# 
# sure.22.combined <- rbind(sure.22.eq, sure.22.p, sure.22.m)
# 
# sure.22.combined.variant <- sure.22.combined[sure.22.combined$SNP_ID != ""]
# 
# 
# sum(grepl("indel", sure.22.combined.variant$SNP_TYPE))/nrow(sure.22.combined.variant)

## I would like to do this in a systematic way. 
dir.vec <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/"
rep.vec <- c(rep("SuRE42-1/",69), rep("SuRE42-2/",69),rep("SuRE43-1/",69),rep("SuRE43-2/",69),rep("SuRE44-1/",69),rep("SuRE44-2/",69),rep("SuRE45-1/",69),rep("SuRE45-2/",69))
par.vec <- c(rep("equal/",23), rep("paternal/",23), rep("maternal/", 23))
chr.vec <- c(1:22,"X")
end.vec <- ".bedpe.gz"

# combine all

file.vec <- paste0(dir.vec, rep.vec, par.vec, chr.vec, end.vec)
print("start loop")
fraction.indel.reads <- foreach(i=1:length(file.vec)) %dopar% {
  
  file.index <- file.vec[i]
  
  sure.reads <- fread(file.index)
  sure.reads.variant <- sure.reads[sure.reads$SNP_ID != ""]
  c(sum(grepl("indel", sure.reads$SNP_TYPE)), nrow(sure.reads.variant), nrow(sure.reads))
}
print("done with loop")

saveRDS(object = fraction.indel.reads, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/fraction.indel.containing.reads.20200109.RDS")
# This generates a list of length(file.vec) objects. for eacht object 2 numbers
# the first one is the amount of variant-containing reads with indels. the second
# one is the total amount of variant containing reads. 
# 
#indel.reads <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/fraction.indel.containing.reads.RDS")
indel.reads <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/fraction.indel.containing.reads.20200109.RDS")


indel.reads.unlist <- unlist(indel.reads)
indel.reads.vec <- indel.reads.unlist[seq(1, length(indel.reads.unlist), by = 3)]
variant.reads.vec <- indel.reads.unlist[seq(2, length(indel.reads.unlist), by = 3)]
all.reads.vec <- indel.reads.unlist[seq(3, length(indel.reads.unlist), by = 3)]

sum(indel.reads.vec)
sum(variant.reads.vec)
sum(all.reads.vec)

sum(indel.reads.vec)/sum(all.reads.vec)
sum(variant.reads.vec)/sum(all.reads.vec)

fraction.indel.reads <- sum(indel.reads)/sum(all.reads)


## some stuff for the report

snps <- table(`All.10-1000.elements.20191113`$snp.type)[2]
indels <- table(`All.10-1000.elements.20191113`$snp.type)[1]

k.snp.raqtl <- table(`K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113`$snp.type)[2]
k.indel.raqtl <- table(`K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113`$snp.type)[1]

h.snp.raqtl <- table(`HepG2.raqtl.10-1000.elements.4.min.SuREexpr.20191113`$snp.type)[2]
h.indel.raqtl <- table(`HepG2.raqtl.10-1000.elements.4.min.SuREexpr.20191113`$snp.type)[1]

k.snp.raqtl/snps
h.snp.raqtl/snps

k.indel.raqtl/indels
h.indel.raqtl/indels
