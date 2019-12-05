
# this script is used to analyse data from the CADD (Combined Annotation Dependent Depletion)
# and see what we can do with this dataset. 

# load required libraries

library(data.table)


# The dataset was downloaded from website of CADD (https://cadd.gs.washington.edu/download).
# File = 1000 Genome variants (SNVs and InDels) from v 1.3 which was retrieved by using
# wget https://krishna.gs.washington.edu/download/CADD/v1.3/1000G.tsv.gz on the linux machine

cadd <- fread("/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/CADD/1000G_inclAnno.tsv.gz")
raqtl <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
raqtl <- raqtl[raqtl$snp.type == "indel",]

# my raqtls are annotated with SNPIds. The other ones are not. I need to match them, but it 
# looks like my own annotation of hg19 is not really good (especially for indels)

## Test CHROM 1

cadd.1 <- cadd[cadd$`#Chrom` == 1 & cadd$Type != "SNV",]
cadd.1[,8:116] <- NULL
raqtl.1 <- raqtl[raqtl$chrom == 1,]

# i do not know wat the difference is between Anc and REF
sum(cadd.1$Ref == cadd.1$Anc, na.rm = T)
sum(cadd.1$Ref != cadd.1$Anc, na.rm = T)

# visual inspectation of data looks like all my hg19 positions should be +1 to match a variant in the CADD database
# lets see if we can annotate this. 

raqtl.1$pos.hg19 <- raqtl.1$pos.hg19 + 1
match.idx <- match(raqtl.1$pos.hg19,cadd.1$Pos)

sum(is.na(match.idx))
#only 40% could be annotated

raqtl.1$ref.seq.cadd <- cadd.1[match.idx,Ref]
raqtl.1$alt.seq.cadd <- cadd.1[match.idx,Alt]
raqtl.1$cadd.tfbs <- cadd.1[match.idx,TFBS]
raqtl.1$cadd.score <- cadd.1[match.idx,RawScore]

raqtl.1[which(raqtl.1$ref.seq != raqtl.1$ref.seq.cadd),]

raqtl.1[,]
