# this script is used to calculate the percentage of common variants (MAF >5%)

# Data is a gz-zipped bcf file downloaded from ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/database/organism_data/
# on 2-12-2019. Downloaded with wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/database/organism_data/b151_SNPChrPosOnRef_105.bcp.gz

# load libraries

library(data.table)



# load dbSNP database

dbsnp <- fread(file = "data/external/dbSNP/b151_SNPChrPosOnRef_105.bcp.gz", nrows = 10000)


# try another ne ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz
# this seems to be the one I need. It is from dbSNP150, although I do not know whether this is the same version as 
# the VCF files for the 4 genomes. It might be that rsIDs have changed over time, but I guess that that is ok. 

dbsnp2 <- fread(file = "data/external/dbSNP/common_all_20170710.vcf.gz")

# G5A is defined as >5% minor allele frequency in each and all populations
# G5 is defined as >5% minor allele frequency in 1+ populations


tst <- dbsnp2[1:10000,]
head(tst)

maf.5.idx <- grep("G5A", dbsnp2$INFO)
maf.5.idx.alternative <- grep(";G5;", dbsnp2$INFO)

dbsnp.maf.5 <- dbsnp2[maf.5.idx,]
dbsnp.maf.5.alternative <- dbsnp2[maf.5.idx.alternative,]

sum(dbsnp.maf.5$ID %in% `All.10-1000.elements.20191113`$SNP_ID)
sum(dbsnp.maf.5.alternative$ID %in% `All.10-1000.elements.20191113`$SNP_ID)

# I will take the >5% MAF as a measurement. Then we have

n_distinct(dbsnp.maf.5.alternative$ID)


## CONCLUSION ##

# Based on GA (MAF >5% in 1+ populations):
# 49.6% of common variants is represented by the library
# 12.8 million common variants are present in dbSNP150
# 6,367,592 variants are common
