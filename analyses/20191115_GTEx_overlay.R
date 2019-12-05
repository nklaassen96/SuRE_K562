# this script is used to compare with the gtex database. 

library(data.table)

# data was downloaded from https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar on 15-11-2019 and the file
# GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz was extracted from the tar.file with 
# `tar -xvf GTEx_Analysis_v8_eQTL.tar GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz`



# the data is aligned to b38 (similar to hg38) but the current snps are aligned to hg19. Therefore I download the older version 

gtex <- fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/GTEx/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
gtex.lookup <- fread("/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")

head(gtex.lookup)

# I only need the first and seventh row to link SNP.ID to variant.id

gtex.lookup[,c(2,3,4,5,6,8)] <- NULL

head(gtex.lookup)
tail(gtex)


match.idx <- match(gtex$variant_id, gtex.lookup$variant_id)

gtex$SNP_ID <- gtex.lookup[match.idx,rs_id_dbSNP151_GRCh38p7]

#saveRDS(object = gtex, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/GTEx/GTEx.Wholeblood.v8.SNP_ID.annotated.RDS")

gtex.annotated <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/GTEx/GTEx.Wholeblood.v8.SNP_ID.annotated.RDS")


# Some snps have as snp.id a ".". We dont need those rows
gtex.annotated <- gtex.annotated[gtex.annotated$SNP_ID != ".",]

{
# SNPs have a strange format, lets change that

splitted.snp.id <- unlist(strsplit(gtex$variant_id, split = "_"))
gtex$chr <- splitted.snp.id[seq(1,length(splitted.snp.id), by = 5)]
gtex$pos <- splitted.snp.id[seq(2,length(splitted.snp.id), by = 5)]
gtex$ref <- splitted.snp.id[seq(3,length(splitted.snp.id), by = 5)]
gtex$alt <- splitted.snp.id[seq(4,length(splitted.snp.id), by = 5)]

# want to remove all the snps, so where ref and alt only have1 nt

var.size <-nchar(gtex$ref)+nchar(gtex$alt) 
table(var.size)

# 9% seems to be an indel

gtex.indels <- gtex[which(var.size != 2),]
}
# I figured that there is a SNP-ID lookup table. So lets do that instead. So we can do it also for the V8 database as mapping doesnt matter. This is still
# a nice way to remove 90% of the data that I will not use anyway. 




raqtl.hepg2 <- readRDS(file="/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/HepG2.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
raqtl.k562 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
all.variants <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")

raqtl.k562.indel <- raqtl.k562[raqtl.k562$snp.type == "indel",]

sum(raqtl.k562$SNP_ID %in% gtex.annotated$SNP_ID)
# 13% of raqtl variants is signifcant eQTL
sum(raqtl.k562.indel$SNP_ID %in% gtex.annotated$SNP_ID)
# 9% for the indels 

# Annotate all variants with GTEx data
all.variants$gtex.slope <- gtex.annotated[match(all.variants$SNP_ID, gtex.annotated$SNP_ID),slope]
all.variants$gtex.tss <- abs(gtex.annotated[match(all.variants$SNP_ID, gtex.annotated$SNP_ID),tss_distance])
all.variants.gtex <- all.variants[!is.na(all.variants$gtex.slope),]

# Annotate K562 raQTLs with GTEx data
match.idx <- match(raqtl.k562$SNP_ID, gtex.annotated$SNP_ID)
raqtl.k562$gtex.slope <- gtex.annotated[match.idx,slope]
raqtl.k562$gtex.tss <- abs(gtex.annotated[match.idx,tss_distance])
raqtl.k562.gtex <- raqtl.k562[!is.na(raqtl.k562$gtex.slope),]

# Annotate HepG2 raQTLs with GTEx data
raqtl.hepg2$gtex.slope <- gtex.annotated[match(raqtl.hepg2$SNP_ID, gtex.annotated$SNP_ID),slope]
raqtl.hepg2$gtex.tss <- abs(gtex.annotated[match(raqtl.hepg2$SNP_ID, gtex.annotated$SNP_ID),tss_distance])
raqtl.hepg2.gtex <- raqtl.hepg2[!is.na(raqtl.hepg2$gtex.slope),]

# Check concordance K562
table(raqtl.k562.gtex$k562.max == "ref" & sign(raqtl.k562.gtex$gtex.slope) == -1) #concordant
table(raqtl.k562.gtex$k562.max == "alt" & sign(raqtl.k562.gtex$gtex.slope) == -1) #disconcordant
table(raqtl.k562.gtex$k562.max == "alt" & sign(raqtl.k562.gtex$gtex.slope) == 1) #concordant
table(raqtl.k562.gtex$k562.max == "ref" & sign(raqtl.k562.gtex$gtex.slope) == 1) #disconcordant

# Annotate concordance K562
 
conc.1 <- which(raqtl.k562.gtex$k562.max == "ref" & sign(raqtl.k562.gtex$gtex.slope) == -1) #concordant
dis.1 <- which(raqtl.k562.gtex$k562.max == "alt" & sign(raqtl.k562.gtex$gtex.slope) == -1) #disconcordant
conc.2 <- which(raqtl.k562.gtex$k562.max == "alt" & sign(raqtl.k562.gtex$gtex.slope) == 1) #concordant
dis.2 <- which(raqtl.k562.gtex$k562.max == "ref" & sign(raqtl.k562.gtex$gtex.slope) == 1) #disconcordant

raqtl.k562.gtex[c(conc.1, conc.2),"gtex.concordance"] <- "concordant"
raqtl.k562.gtex[c(dis.1, dis.2),"gtex.concordance"] <- "disconcordant"

# annotate concordance All Variants

conc.1.all <- which(all.variants.gtex$k562.max == "ref" & sign(all.variants.gtex$gtex.slope) == -1) #concordant
dis.1.all <- which(all.variants.gtex$k562.max == "alt" & sign(all.variants.gtex$gtex.slope) == -1) #disconcordant
conc.2.all <- which(all.variants.gtex$k562.max == "alt" & sign(all.variants.gtex$gtex.slope) == 1) #concordant
dis.2.all <- which(all.variants.gtex$k562.max == "ref" & sign(all.variants.gtex$gtex.slope) == 1) #disconcordant

all.variants.gtex[c(conc.1.all, conc.2.all),"gtex.concordance"] <- "concordant"
all.variants.gtex[c(dis.1.all, dis.2.all),"gtex.concordance"] <- "disconcordant"

# Annotate concordance HepG2

conc.1.hepg2 <- which(raqtl.hepg2.gtex$hepg2.max == "ref" & sign(raqtl.hepg2.gtex$gtex.slope) == -1) #concordant
dis.1.hepg2 <- which(raqtl.hepg2.gtex$hepg2.max == "alt" & sign(raqtl.hepg2.gtex$gtex.slope) == -1) #disconcordant
conc.2.hepg2 <- which(raqtl.hepg2.gtex$hepg2.max == "alt" & sign(raqtl.hepg2.gtex$gtex.slope) == 1) #concordant
dis.2.hepg2 <- which(raqtl.hepg2.gtex$hepg2.max == "ref" & sign(raqtl.hepg2.gtex$gtex.slope) == 1) #disconcordant

raqtl.hepg2.gtex[c(conc.1.hepg2, conc.2.hepg2),"gtex.concordance"] <- "concordant"
raqtl.hepg2.gtex[c(dis.1.hepg2, dis.2.hepg2),"gtex.concordance"] <- "disconcordant"



hist(log10(raqtl.k562.gtex$gtex.tss), breaks = 100, col = 1)
hist(raqtl.k562.gtex$gtex.tss, breaks = 100, col = 1)

tst <- raqtl.k562.gtex[raqtl.k562.gtex$gtex.tss < 1000,]

table(tst$k562.max == "ref" & sign(tst$gtex.slope) == -1) #concordant
table(tst$k562.max == "alt" & sign(tst$gtex.slope) == 1) #concordant





#K562
# I want to calculate (for various distances to the tss) the concordance of these raqtls vs the concordance
# in the whole set (all.variants.gtex), However I am not sure if this is the right way to go about this. 

# I need to remove the raqtls from the all.variants.gtex dataframe
which(!all.variants.gtex$SNP_ID %in% raqtl.k562.gtex$SNP_ID)
all.variants.gtex.filtered <- all.variants.gtex[which(!all.variants.gtex$SNP_ID %in% raqtl.k562.gtex$SNP_ID),]

tss.distance.vector <- 10^seq(3,6,length.out = 30)
p.value.vector <- NULL
odds.ratio.vector <- NULL

for (i in c(1:length(tss.distance.vector))){
  
  print(i)
  
  # Create a df with only the raQTLs of K562 within the selected TSS distance
  df.raqtl <- raqtl.k562.gtex[raqtl.k562.gtex$gtex.tss < tss.distance.vector[i],]
  
  # Create a df with all variants within the selected TSS distance
  df.all <- all.variants.gtex.filtered[all.variants.gtex.filtered$gtex.tss < tss.distance.vector[i],]
  
  f <- fisher.test(rbind(table(df.raqtl$gtex.concordance),table(df.all$gtex.concordance)), alternative = "g") 
  p.value.vector[i] <- f$p.value
  odds.ratio.vector[i] <- f$estimate
}

plot(log10(tss.distance.vector), odds.ratio.vector, ylim = c(1,1.6), xlab = "log10(maximum distance to TSS [bp])", ylab = "concordance odds ratio", pch = 19)
abline(h=1, lty=2)




###  HepG2  ###



# Preparation 
'
gtex.hepg2 <- fread("/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/GTEx/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz")
gtex.lookup <- fread("/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
gtex.lookup[,c(2,3,4,5,6,8)] <- NULL
match.idx <- match(gtex.hepg2$variant_id, gtex.lookup$variant_id)
gtex.hepg2$SNP_ID <- gtex.lookup[match.idx,rs_id_dbSNP151_GRCh38p7]
saveRDS(object = gtex.hepg2, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/GTEx/GTEx.Liver.v8.SNP_ID.annotated.RDS")
'

gtex.hepg2 <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/GTEx/GTEx.Liver.v8.SNP_ID.annotated.RDS")

all.variants <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")
raqtl.hepg2 <- readRDS(file="/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/HepG2.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")

# Annotate all variants with GTEx data
all.variants$gtex.slope <- gtex.hepg2[match(all.variants$SNP_ID, gtex.hepg2$SNP_ID),slope]
all.variants$gtex.tss <- abs(gtex.hepg2[match(all.variants$SNP_ID, gtex.hepg2$SNP_ID),tss_distance])
all.variants.gtex <- all.variants[!is.na(all.variants$gtex.slope),]

# Annotate HepG2 raQTLs with GTEx data
raqtl.hepg2$gtex.slope <- gtex.hepg2[match(raqtl.hepg2$SNP_ID, gtex.hepg2$SNP_ID),slope]
raqtl.hepg2$gtex.tss <- abs(gtex.hepg2[match(raqtl.hepg2$SNP_ID, gtex.hepg2$SNP_ID),tss_distance])
raqtl.hepg2.gtex <- raqtl.hepg2[!is.na(raqtl.hepg2$gtex.slope),]


# annotate concordance All Variants

conc.1.all <- which(all.variants.gtex$hepg2.max == "ref" & sign(all.variants.gtex$gtex.slope) == -1) #concordant
dis.1.all <- which(all.variants.gtex$hepg2.max == "alt" & sign(all.variants.gtex$gtex.slope) == -1) #disconcordant
conc.2.all <- which(all.variants.gtex$hepg2.max == "alt" & sign(all.variants.gtex$gtex.slope) == 1) #concordant
dis.2.all <- which(all.variants.gtex$hepg2.max == "ref" & sign(all.variants.gtex$gtex.slope) == 1) #disconcordant

all.variants.gtex[c(conc.1.all, conc.2.all),"gtex.concordance"] <- "concordant"
all.variants.gtex[c(dis.1.all, dis.2.all),"gtex.concordance"] <- "disconcordant"

# Annotate concordance HepG2

conc.1.hepg2 <- which(raqtl.hepg2.gtex$hepg2.max == "ref" & sign(raqtl.hepg2.gtex$gtex.slope) == -1) #concordant
dis.1.hepg2 <- which(raqtl.hepg2.gtex$hepg2.max == "alt" & sign(raqtl.hepg2.gtex$gtex.slope) == -1) #disconcordant
conc.2.hepg2 <- which(raqtl.hepg2.gtex$hepg2.max == "alt" & sign(raqtl.hepg2.gtex$gtex.slope) == 1) #concordant
dis.2.hepg2 <- which(raqtl.hepg2.gtex$hepg2.max == "ref" & sign(raqtl.hepg2.gtex$gtex.slope) == 1) #disconcordant

raqtl.hepg2.gtex[c(conc.1.hepg2, conc.2.hepg2),"gtex.concordance"] <- "concordant"
raqtl.hepg2.gtex[c(dis.1.hepg2, dis.2.hepg2),"gtex.concordance"] <- "disconcordant"




all.variants.gtex.filtered.hepg2 <- all.variants.gtex[which(!all.variants.gtex$SNP_ID %in% raqtl.hepg2.gtex$SNP_ID),]

tss.distance.vector.h <- 10^seq(2,6,length.out = 30)
p.value.vector.h <- NULL
odds.ratio.vector.h <- NULL

for (i in c(1:length(tss.distance.vector.h))){
  
  print(i)
  
  # Create a df with only the raQTLs of K562 within the selected TSS distance
  df.raqtl <- raqtl.hepg2.gtex[raqtl.hepg2.gtex$gtex.tss < tss.distance.vector.h[i],]
  
  # Create a df with all variants within the selected TSS distance
  df.all <- all.variants.gtex.filtered.hepg2[all.variants.gtex.filtered.hepg2$gtex.tss < tss.distance.vector.h[i],]
  
  f <- fisher.test(rbind(table(df.raqtl$gtex.concordance),table(df.all$gtex.concordance)), alternative = "g") 
  p.value.vector.h[i] <- f$p.value
  odds.ratio.vector.h[i] <- f$estimate
}
png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/GTEx/hepg2.concordance.tss.png")
plot(log10(tss.distance.vector.h), odds.ratio.vector.h, ylim = c(1,2.4),xlab = "log10(maximum distance to TSS [bp])", ylab = "Concordance odds ratio", pch = 19, main = "HepG2")
abline(h=1, lty=2)
dev.off()

## Hepg2 makes no sense because wholeblood was used as a sample. 
