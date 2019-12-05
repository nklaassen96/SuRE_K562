# try to find important GWAS indel
library(data.table)

# data downloaded fro https://www.ebi.ac.uk/gwas/docs/file-downloads and then the link https://www.ebi.ac.uk/gwas/api/search/downloads/full
# was used with wget. 

gwas <- fread("/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/GWAS-NHGRI-EBI/full.tsv")
raqtl <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")

raqtl.in.gwas <- raqtl[which(raqtl$SNP_ID %in% gwas$SNPS),]
table(raqtl.in.gwas$snp.type)

indel.id <- raqtl.in.gwas[raqtl.in.gwas$snp.type == "indel",]$SNP_ID

gwas.indel <- gwas[which(gwas$SNPS %in% indel.id)]

raqtl.indel.in.gwas <- raqtl.in.gwas[raqtl.in.gwas$SNP_ID %in% indel.id,]
