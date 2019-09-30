# I use the counts sample, but could be the whole table containing all variants.

counts <- counts.all.chrom
head(counts)
table(counts$SNP_VAR)

snp.id <- which(tapply(counts$SNP_VAR, counts$SNP_ID, function(x) {all(c(0,1) %in% x)} ) == TRUE)

snps <- counts.all.chrom[counts.all.chrom$SNP_ID %in% names(snp.id),]
