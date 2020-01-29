
# this script is used to analyse data from the CADD (Combined Annotation Dependent Depletion)
# and see what we can do with this dataset. 

# The hypothesis is that raQTLs have a more severe score than the other variants, if that is taken into account

# load required libraries

library(data.table)


# The dataset was downloaded from website of CADD (https://cadd.gs.washington.edu/download).
# File = 1000 Genome variants (SNVs and InDels) from v 1.3 which was retrieved by using
# wget https://krishna.gs.washington.edu/download/CADD/v1.3/1000G.tsv.gz on the linux machine

cadd <- fread("/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/CADD/1000G_inclAnno.tsv.gz")
raqtl.k562 <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
#raqtl <- raqtl[raqtl$snp.type == "indel",]
all.variants <- readRDS("data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")

# my raqtls are annotated with SNPIds. The other ones are not. I need to match them, but it 
# looks like my own annotation of hg19 is not really good (especially for indels)

## Test CHROM 1

cadd.1 <- cadd[cadd$`#Chrom` == 1,]
cadd.1[,8:114] <- NULL

all.variants.1 <- all.variants[all.variants$chr == 1,]
raqtl.1 <- raqtl.k562[raqtl.k562$chrom == 1,]

# i do not know wat the difference is between Anc and REF
sum(cadd.1$Ref == cadd.1$Anc, na.rm = T)
sum(cadd.1$Ref != cadd.1$Anc, na.rm = T)

# visual inspectation of data looks like all my hg19 positions should be +1 to match a variant in the CADD database
# lets see if we can annotate this. 

raqtl.1$pos.hg19 <- raqtl.1$pos.hg19 + 1
match.idx.raqtl <- match(raqtl.1$pos.hg19,cadd.1$Pos)

# maybe improve the matching to include the allele

sum(is.na(match.idx.raqtl))
#only 40% could be annotated

raqtl.1$ref.seq.cadd <- cadd.1[match.idx.raqtl,Ref]
raqtl.1$alt.seq.cadd <- cadd.1[match.idx.raqtl,Alt]
#raqtl.1$cadd.tfbs <- cadd.1[match.idx.raqtl,TFBS]
raqtl.1$cadd.score <- cadd.1[match.idx.raqtl,RawScore]

raqtl.1[which(raqtl.1$ref.seq != raqtl.1$ref.seq.cadd),]
raqtl.1[which(raqtl.1$alt.seq != raqtl.1$alt.seq.cadd),]



# also annotate all variants

all.variants.1$pos.hg19 <- all.variants.1$pos.hg19 + 1
match.idx.all.variants <- match(all.variants.1$pos.hg19,cadd.1$Pos)

# maybe improve the matching to include the allele

  
  sum(is.na(match.idx.all.variants))
#only 40% could be annotated

all.variants.1$ref.seq.cadd <- cadd.1[match.idx.all.variants,Ref]
all.variants.1$alt.seq.cadd <- cadd.1[match.idx.all.variants,Alt]
#all.variants.1$cadd.tfbs <- cadd.1[match.idx.all.variants,TFBS]
all.variants.1$cadd.score <- cadd.1[match.idx.all.variants,RawScore]

all.variants.1[which(all.variants.1$ref.seq != all.variants.1$ref.seq.cadd),]
all.variants.1[which(all.variants.1$alt.seq != all.variants.1$alt.seq.cadd),]

#Lets compare them!

summary(raqtl.1$cadd.score)
summary(all.variants.1$cadd.score)
boxplot(raqtl.1$cadd.score)
boxplot(all.variants.1$cadd.score, add = T)

hist(raqtl.1$cadd.score, breaks = 100, col = 3)
var.sample <- sample(all.variants.1$cadd.score, size = length(discard(raqtl.1$cadd.score, is.na)))
hist(var.sample, add = T, col = alpha(2, 0.5), breaks = 150)
plot(-log10(raqtl.1$K562.wilcoxon.pvalue), raqtl.1$cadd.score)
     