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

saveRDS(object = gtex, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/GTEx/GTEx.Wholeblood.v8.SNP_ID.annotated.RDS")

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

# I figured that there is a SNP-ID lookup table. So lets do that instead. So we can do it also for the V8 database as mapping doesnt matter. This is still
# a nice way to remove 90% of the data that I will not use anyway. 





raqtl.k562 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
all.variants <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")

sum(raqtl.k562$SNP_ID %in% gtex$SNP_ID)
# 13% of raqtl variants is signifcant eQTL