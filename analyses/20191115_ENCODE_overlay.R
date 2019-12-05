# This script will be used to overlay the ENCODE element annotations with the identified raQTLs and their position. 
# Importantly this is only for K562!

# load required libraries
library(data.table)
library(GenomicRanges)


# ENCODE have made annotations of genomic regions which I downloaded with  wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationCombinedK562.bed.gz
# on the LINUX machine server. downloaded on 15-11-2019

encode.annotation <- fread("/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/ENCODE/wgEncodeAwgSegmentationCombinedK562.bed.gz")

# lets inspect this a little bit more

table(encode.annotation$V4)

# So there are 7 different annotations. 
# Column 2 and 3 look exactly similar to 8 and 9, is this the case?

sum(encode.annotation$V2 != encode.annotation$V7)
sum(encode.annotation$V3 != encode.annotation$V8)

# Yes, so these columns can be removed. Does column #5 only contain "1000"?

table(encode.annotation$V5)

#Yes, so that one can also be removed. What is column #9?

table(encode.annotation$V9)

# this is the color in RGB terms, not interesting, throw away

encode.annotation <- encode.annotation[,1:4]

# to overlay this with the indels, I need to transform then to GRanges objects and then determine the overlap. 

annotation.granges <- GRanges(seqnames = encode.annotation$V1, ranges = IRanges(start = encode.annotation$V2, end = encode.annotation$V3))
annotation.granges$element <- encode.annotation$V4







## Load my indels and raqtls

all.variants <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")


#Load the raQTL sets
raqtl.k562 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
' # Test for a loop-variant. Later I made it working vectorized.
for (i in c(1:nrow(raqtl.hepg2.indel))){
  
  print(i)
  
  #first i make a granges object for that specific variant
  granges.idx <- GRanges(seqnames = paste0("chr",raqtl.hepg2.indel[i,"chrom"]), ranges = IRanges(start = raqtl.hepg2.indel[i,"pos.hg19"],
                        end = raqtl.hepg2.indel[i,"pos.hg19"]))
  
  #then i find the overlaps within the annotated granges
  f <- findOverlaps(query = annotation.granges, subject = granges.idx)
  
  # the queryhits gives the row from which the annotation.granges was found
  q <- queryHits(f)
  
  if(length(q) > 0){raqtl.hepg2.indel[i,"encode.annotation"] <- annotation.granges[q]$element} else {}
}
'

' ### This was a test to see if i get the same results. 
# This works, but I have to make it vectorized for speed. 
granges.all <- GRanges(seqnames = paste0("chr", raqtl.hepg2.indel$chrom), ranges = IRanges(start = raqtl.hepg2.indel$pos.hg19, end = raqtl.hepg2.indel$pos.hg19))
# first lets check with how many annotated elements the granges overlaps. should be only 1 or 0. this is the case. 
table(countOverlaps(granges.all, annotation.granges))

f <- findOverlaps(query = granges.all, subject = annotation.granges)
# contains the rows of the annotated variants
s <- subjectHits(f)
q <- queryHits(f)

raqtl.hepg2.indel[q,"annotation.test"] <- annotation.granges[s]$element

# I think this also works, lets check this.

table(raqtl.hepg2.indel$encode.annotation == raqtl.hepg2.indel$annotation.test)
'

## (1) All variants ~7,200,000

granges.all <- GRanges(seqnames = paste0("chr", all.variants$chrom), ranges = IRanges(start = all.variants$pos.hg19, end = all.variants$pos.hg19))
table(countOverlaps(granges.all,annotation.granges))

idx.2.overlaps <- which(countOverlaps(granges.all,annotation.granges) == 2)
annotation.granges[idx.2.overlaps]

#some elements overlap with 2 annotations. lets see how we solve this. There is no error. but there has to be right? I accept that there might be a very small percentage 2200 of 7 milion wrong annotations

f.all <- findOverlaps(query = granges.all, subject = annotation.granges, select = "all")
# contains the rows of the annotated variants
s.all <- subjectHits(f.all)
q.all <- queryHits(f.all)

all.variants[q.all,"encode.annotation"] <- as.factor(annotation.granges[s.all]$element)
all.variants.indel <- all.variants[all.variants$snp.type == "indel",]

## (2) K562 raQTL all ~30,862

granges.k562 <- GRanges(seqnames = paste0("chr", raqtl.k562$chrom), ranges = IRanges(start = raqtl.k562$pos.hg19, end = raqtl.k562$pos.hg19))
table(countOverlaps(granges.k562, annotation.granges))
#some elements overlap with 2 annotations. lets see how we solve this. There is no error. but there has to be right?

f.k562 <- findOverlaps(query = granges.k562, subject = annotation.granges, select = "all")
# contains the rows of the annotated variants
s.k562 <- subjectHits(f.k562)
q.k562 <- queryHits(f.k562)

raqtl.k562[q.k562,"encode.annotation"] <- as.factor(annotation.granges[s.k562]$element)
raqtl.k562.indel <- raqtl.k562[raqtl.k562$snp.type == "indel",]

## Now is is time to generate some figures

col.k562 <- "steelblue"



levels(raqtl.k562$encode.annotation) == levels(all.variants$encode.annotation)

fraction.elements.all <- table(all.variants$encode.annotation)/sum(table(all.variants$encode.annotation))
fraction.elements.k562 <- table(raqtl.k562$encode.annotation)/sum(table(raqtl.k562$encode.annotation))

## fig. 1 all variant enrichtments
png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/ENCODE.annotation/Fig.1.all.variants.encode.enrichtment.png")
fold.enrichment <- fraction.elements.k562/fraction.elements.all
par(mar=c(8,5,4,2))
b <- barplot(main = "All variants",ylim = c(0,12),xaxt= "n",fold.enrichment, ylab = "Fold enrichment", col = col.k562)
text(pos = 2, cex = 1, x=b+0.3, y = -0.4,c("CTCF", "Enhancer", "Promoter flanking", "Repressed", "Transcribed", "Transcription start site", "Weak enhancer"), xpd = T, srt = 45 )
abline(h = 1, lty = 3)
dev.off()

fraction.elements.indel.all <- table(all.variants.indel$encode.annotation)/sum(table(all.variants.indel$encode.annotation))
fraction.elements.indel.k562 <- table(raqtl.k562.indel$encode.annotation)/sum(table(raqtl.k562.indel$encode.annotation))

## fig. 2 indel variant enrichtments
png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/ENCODE.annotation/Fig.2.indel.encode.enrichtment.png")
fold.enrichment.indel <- fraction.elements.indel.k562/fraction.elements.indel.all
par(mar=c(8,5,4,2))
b <- barplot(main = "Indels", ylim = c(0,12), xaxt= "n",col = col.k562,fold.enrichment.indel, ylab = "Fold enrichment", add = F)
text(pos = 2, cex = 1, x=b+0.3, y = -0.4,c("CTCF", "Enhancer", "Promoter flanking", "Repressed", "Transcribed", "Transcription start site", "Weak enhancer"), xpd = T, srt = 45 )
abline(h = 1, lty = 3)
dev.off()

raqtl.k562.snp <- raqtl.k562[raqtl.k562$snp.type == "snp",]
all.variants.snp <- all.variants[all.variants$snp.type == "snp",]

fraction.elements.snp.all <- table(all.variants.snp$encode.annotation)/sum(table(all.variants.snp$encode.annotation))
fraction.elements.snp.k562 <- table(raqtl.k562.snp$encode.annotation)/sum(table(raqtl.k562.snp$encode.annotation))

## fig. 3 snp variant enrichment
png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/ENCODE.annotation/Fig.3.snp.encode.enrichtment.NOT.NEEDED.png")
fold.enrichment.snp <- fraction.elements.snp.k562/fraction.elements.snp.all
par(mar=c(8,5,4,2))
b <- barplot(ylim = c(0,12), xaxt= "n",col = col.k562,fold.enrichment.snp, ylab = "Fold enrichment",col = alpha(3, 0.5), add = T)
text(pos = 2, cex = 1, x=b+0.3, y = -0.4,c("CTCF", "Enhancer", "Promoter flanking", "Repressed", "Transcribed", "Transcription start site", "Weak enhancer"), xpd = T, srt = 45 )
abline(h = 1, lty = 3)
dev.off()

# It looked like there was a significant difference in the fold enrichment for promoter flanking regions.
# however, there are only 3 indels that overlap with these PF regions. 

barplot(fraction.elements.indel.all, ylim = c(0,0.01))
barplot(fraction.elements.snp.all, add = T, col = alpha(2, 0.5))

barplot(fraction.elements.indel.k562, ylim = c(0,0.1))
barplot(fraction.elements.snp.k562, col = alpha(2, 0.5), add = T)

## Fig. 4 SNP enrichtment compared to indels.

table(all.variants$encode.annotation, useNA = "always")

# first only select the variants that could be annotated (~95%)
all.variants.annotated <- all.variants[!is.na(all.variants$encode.annotation),]
ind.pf <- sum(all.variants.annotated$snp.type == "indel" & all.variants.annotated$encode.annotation == "PF")/sum(all.variants.annotated$encode.annotation =="PF") #snp fraction in PF
ind.e <- sum(all.variants.annotated$snp.type == "indel" & all.variants.annotated$encode.annotation == "E")/sum(all.variants.annotated$encode.annotation =="E")
ind.r <- sum(all.variants.annotated$snp.type == "indel" & all.variants.annotated$encode.annotation == "R")/sum(all.variants.annotated$encode.annotation =="R")
ind.ctcf <- sum(all.variants.annotated$snp.type == "indel" & all.variants.annotated$encode.annotation == "CTCF")/sum(all.variants.annotated$encode.annotation =="CTCF")
ind.we <- sum(all.variants.annotated$snp.type == "indel" & all.variants.annotated$encode.annotation == "WE")/sum(all.variants.annotated$encode.annotation =="WE")
ind.tss <- sum(all.variants.annotated$snp.type == "indel" & all.variants.annotated$encode.annotation == "TSS")/sum(all.variants.annotated$encode.annotation =="TSS")
ind.t <- sum(all.variants.annotated$snp.type == "indel" & all.variants.annotated$encode.annotation == "T")/sum(all.variants.annotated$encode.annotation =="T")
ind.all <- sum(all.variants.annotated$snp.type == "indel")/nrow(all.variants.annotated)

indel.fractions <- c(ind.ctcf, ind.e, ind.pf, ind.r, ind.t, ind.tss, ind.we)

png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/ENCODE.annotation/Fig.4.indel.fractions.per.annotation.png")
par(mar=c(8,5,4,2))
b <- barplot(indel.fractions, ylab = "Fraction of Indels", ylim = c(0,0.14), col = col.k562)
text(pos = 2, cex = 1, x=b+0.4, y = -0.004,c("CTCF", "Enhancer", "Promoter flanking", "Repressed", "Transcribed", "Transcription start site", "Weak enhancer"), xpd = T, srt = 45 )
abline(h = ind.all, lty = 3)
text(b[c(1,2,4,5,7)], y = indel.fractions[c(1,2,4,5,7)]+0.004, labels = "*",cex = 2 )
dev.off()

#looks pretty good, now calculate p values a loop

annotation.vector <- c("CTCF", "E", "PF", "R", "T", "TSS", "WE")
p.value.vector <- rep(NA,7)
names(p.value.vector) <- annotation.vector
odds.ratio.vector <- rep(NA,7)
names(odds.ratio.vector) <- annotation.vector

for (i in c(1:7)){
  
  annotation.id <- annotation.vector[i]
  print(annotation.id)
  specific.variant <- table(all.variants.annotated[all.variants.annotated$encode.annotation == annotation.id,]$snp.type) #indel/snp for annotation.id
  not.specific.variant <- table(all.variants.annotated[all.variants.annotated$encode.annotation != annotation.id,]$snp.type) #indel/snp for everything that is not annotation.id
  contingency.matrix <- rbind(specific.variant, not.specific.variant)
  p.value.vector[i]<- fisher.test(contingency.matrix)$p.value
  odds.ratio.vector[i] <- fisher.test(contingency.matrix)$estimate
}


