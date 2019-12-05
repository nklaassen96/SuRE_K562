
# load libraries

library(data.table)
library(tidyverse)


# load raQTL set

raqtl.k562 <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
raqtl.hepg2 <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/HepG2.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
all.variants <- readRDS("data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")

col.k562 <- "steelblue"
col.hepg2 <- "magenta4"





# load SNP2TFBS file

#data is a gz-zipped txt file downloaded from ftp://ccg.vital-it.ch/snp2tfbs/mapped_files/ at 28-11-2019. file  = snp2tfbs_JASPAR_CORE_2014_vert.txt.gz

snp2tfbs.df.2 <- read.delim("/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/snp2tfbs_JASPAR_CORE_2014_vert.txt.gz", header = FALSE, as.is = TRUE)

snp2tfbs.df.filtered <- snp2tfbs.df.2[grepl("rs", snp2tfbs.df.2$V1),]
# alternative rows are removed. 

df <- snp2tfbs.df.filtered
df$V10 <- as.numeric(df$V10)


df$ref.stronger.binding <- NA
# look at the scoredifferences column
df[df$V10 > 0,"ref.stronger.binding"] <- 0
df[df$V10 < 0,"ref.stronger.binding"] <- 1

# some of them have dots (no PWM)
df[df$V8 == ".","ref.stronger.binding"] <- 0
df[df$V9 == ".","ref.stronger.binding"] <- 1
sum(table(df$ref.stronger.binding))

# some have the same value, remove those
df <- df[df$V8 != df$V9,]

table(df$ref.stronger.binding)

raqtl.k562$ref.stronger.binding <- df[match(raqtl.k562$SNP_ID, df$V1),"ref.stronger.binding"]
raqtl.hepg2$ref.stronger.binding <- df[match(raqtl.hepg2$SNP_ID, df$V1),"ref.stronger.binding"]
all.variants$ref.stronger.binding <- df[match(all.variants$SNP_ID, df$V1), "ref.stronger.binding"]


#split files for indels.
all.variants.indel <- all.variants[all.variants$snp.type == "indel",]
raqtl.k562.indel <- raqtl.k562[raqtl.k562$snp.type == "indel",]
raqtl.hepg2.indel <- raqtl.hepg2[raqtl.hepg2$snp.type == "indel",]

all.variants.snp <- all.variants[all.variants$snp.type == "snp",]
raqtl.k562.snp <- raqtl.k562[raqtl.k562$snp.type == "snp",]
raqtl.hepg2.snp <- raqtl.hepg2[raqtl.hepg2$snp.type == "snp",]



table(raqtl.k562$ref.stronger.binding)
sum(
sum(raqtl.k562$k562.max == "ref" & raqtl.k562$ref.stronger.binding == 1, na.rm = TRUE), #conc
sum(raqtl.k562$k562.max == "alt" & raqtl.k562$ref.stronger.binding == 0, na.rm = TRUE), #conc
sum(raqtl.k562$k562.max == "ref" & raqtl.k562$ref.stronger.binding == 0, na.rm = TRUE), #disconc
sum(raqtl.k562$k562.max == "alt" & raqtl.k562$ref.stronger.binding == 1, na.rm = TRUE) #disconc
)

#so we covered everything already with the first two concordance







## SNPs and Indels combined

snpindel.conc.k562 <- sum(
  sum(raqtl.k562$k562.max == "ref" & raqtl.k562$ref.stronger.binding == 1, na.rm = TRUE), 
  sum(raqtl.k562$k562.max == "alt" & raqtl.k562$ref.stronger.binding == 0, na.rm = TRUE)) / sum(table(raqtl.k562$ref.stronger.binding))

snpindel.conc.hepg2 <- sum(
  sum(raqtl.hepg2$hepg2.max == "ref" & raqtl.hepg2$ref.stronger.binding == 1, na.rm = T),
  sum(raqtl.hepg2$hepg2.max == "alt" & raqtl.hepg2$ref.stronger.binding == 0, na.rm = T)) / sum(table(raqtl.hepg2$ref.stronger.binding))

# For the concordance of all I chose a different strategy than Joris (he compares with K562 only). I generate concordance with hepg2 and concordance with k562

snpindel.conc.all.k562 <- sum(
  sum(all.variants$k562.max == "ref" & all.variants$ref.stronger.binding == 1, na.rm = TRUE),
  sum(all.variants$k562.max == "alt" & all.variants$ref.stronger.binding == 0, na.rm = TRUE))
snpindel.conc.all.hepg2 <- sum(
  sum(all.variants$hepg2.max == "ref" & all.variants$ref.stronger.binding == 1, na.rm = TRUE),
  sum(all.variants$hepg2.max == "alt" & all.variants$ref.stronger.binding == 0, na.rm = TRUE))
snpindel.conc.all <- (snpindel.conc.all.hepg2+snpindel.conc.all.k562)/(2*sum(table(all.variants$ref.stronger.binding)))





## SNP only

snp.conc.k562 <- sum(
  sum(raqtl.k562.snp$k562.max == "ref" & raqtl.k562.snp$ref.stronger.binding == 1, na.rm = TRUE), 
  sum(raqtl.k562.snp$k562.max == "alt" & raqtl.k562.snp$ref.stronger.binding == 0, na.rm = TRUE)) / sum(table(raqtl.k562.snp$ref.stronger.binding))

snp.conc.hepg2 <- sum(
  sum(raqtl.hepg2.snp$hepg2.max == "ref" & raqtl.hepg2.snp$ref.stronger.binding == 1, na.rm = T),
  sum(raqtl.hepg2.snp$hepg2.max == "alt" & raqtl.hepg2.snp$ref.stronger.binding == 0, na.rm = T)) / sum(table(raqtl.hepg2.snp$ref.stronger.binding))

# For the concordance of all I chose a different strategy than Joris (he compares with K562 only). I generate concordance with hepg2 and concordance with k562

snp.conc.all.k562 <- sum(
  sum(all.variants.snp$k562.max == "ref" & all.variants.snp$ref.stronger.binding == 1, na.rm = TRUE),
  sum(all.variants.snp$k562.max == "alt" & all.variants.snp$ref.stronger.binding == 0, na.rm = TRUE))
snp.conc.all.hepg2 <- sum(
  sum(all.variants.snp$hepg2.max == "ref" & all.variants.snp$ref.stronger.binding == 1, na.rm = TRUE),
  sum(all.variants.snp$hepg2.max == "alt" & all.variants.snp$ref.stronger.binding == 0, na.rm = TRUE))
snp.conc.all <- (snp.conc.all.hepg2+snp.conc.all.k562)/(2*sum(table(all.variants.snp$ref.stronger.binding)))




## Indel only

indel.conc.k562 <- sum(
  sum(raqtl.k562.indel$k562.max == "ref" & raqtl.k562.indel$ref.stronger.binding == 1, na.rm = TRUE), 
  sum(raqtl.k562.indel$k562.max == "alt" & raqtl.k562.indel$ref.stronger.binding == 0, na.rm = TRUE)) / sum(table(raqtl.k562.indel$ref.stronger.binding))

indel.conc.hepg2 <- sum(
  sum(raqtl.hepg2.indel$hepg2.max == "ref" & raqtl.hepg2.indel$ref.stronger.binding == 1, na.rm = T),
  sum(raqtl.hepg2.indel$hepg2.max == "alt" & raqtl.hepg2.indel$ref.stronger.binding == 0, na.rm = T)) / sum(table(raqtl.hepg2.indel$ref.stronger.binding))

# For the concordance of all I chose a different strategy than Joris (he compares with K562 only). I generate concordance with hepg2 and concordance with k562

indel.conc.all.k562 <- sum(
  sum(all.variants.indel$k562.max == "ref" & all.variants.indel$ref.stronger.binding == 1, na.rm = TRUE),
  sum(all.variants.indel$k562.max == "alt" & all.variants.indel$ref.stronger.binding == 0, na.rm = TRUE))
indel.conc.all.hepg2 <- sum(
  sum(all.variants.indel$hepg2.max == "ref" & all.variants.indel$ref.stronger.binding == 1, na.rm = TRUE),
  sum(all.variants.indel$hepg2.max == "alt" & all.variants.indel$ref.stronger.binding == 0, na.rm = TRUE))
indel.conc.all <- (indel.conc.all.hepg2+indel.conc.all.k562)/(2*sum(table(all.variants.indel$ref.stronger.binding)))

## Figures

# Figure 4d concordance indels new method

png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/SNP2TFBS/Fig.4d.Motif.concordance.indels.0.not.ignored.png")
par(mar=c(4,10,20,1)+.1)
bp <- barplot(c(indel.conc.hepg2, indel.conc.k562, indel.conc.all), names.arg = c("HepG2 Indel raQTLs ", "K562 Indel raQTLs ", "All Indels tested"),
              horiz = TRUE, 
              xlab = "Fraction concordance SuRE and motif disruption", 
              col = c(col.hepg2, col.k562, "black"), 
              las = 1, 
              cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,1))
barplot(add = T,rep(0.5,3), names.arg = NULL,
        horiz = TRUE, 
        col = c("#D199D1", "#B5CDE1", "#999999"), 
        las = 1, axisnames = FALSE,axes = FALSE,
        cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,1), )
text(x = 0.25, y = bp[,1]-0.02, "Expected by chance", col = 1, cex = 1)
text(c(indel.conc.hepg2, indel.conc.k562, indel.conc.all)[1]+0.02, y = bp[1,], labels = "*",cex = 2)
text(c(indel.conc.hepg2, indel.conc.k562, indel.conc.all)[2]+0.02, y = bp[2,], labels = "*",cex = 2)
dev.off()

# Figure 4c concordance snps new method

png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/SNP2TFBS/Fig.4c.Motif.concordance.snp.0.not.ignored.png")
par(mar=c(4,10,20,1)+.1)
bp <- barplot(c(snp.conc.hepg2, snp.conc.k562, snp.conc.all), names.arg = c("HepG2 SNP raQTLs ", "K562 SNP raQTLs ", "All SNPs tested"),
              horiz = TRUE, 
              xlab = "Fraction concordance SuRE and motif disruption", 
              col = c(col.hepg2, col.k562, "black"), 
              las = 1, 
              cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,1))
barplot(add = T,rep(0.5,3), names.arg = NULL,
        horiz = TRUE, 
        col = c("#D199D1", "#B5CDE1", "#999999"), 
        las = 1, axisnames = FALSE,axes = FALSE,
        cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,1), )
text(x = 0.25, y = bp[,1]-0.02, "Expected by chance", col = 1, cex = 1)
text(c(snp.conc.hepg2, snp.conc.k562, snp.conc.all)[1]+0.02, y = bp[1,], labels = "*",cex = 2)
text(c(snp.conc.hepg2, snp.conc.k562, snp.conc.all)[2]+0.02, y = bp[2,], labels = "*",cex = 2)
dev.off()


## Generate a nice example figure

snp.id <- "rs200054879"

# Indel is located on chr15, so I should load that specific file

reads.chr.15 <- readRDS("data/interim/SuRE_Indels_gDNA_Count/SuRE_combined/sure.snp_indel.combined.chrom.15.RDS")

reads.chr.15.filtered <- reads.chr.15[reads.chr.15$SNP_ID == snp.id,]
d <- reads.chr.15.filtered



png(file = "data/processed/Figures/Indel_Example/Fig.1.sure.signal.per.allele.png")
#par(mar=c(4,6,4,0),)
par(mfcol = c(1,2), mar=c(4,4,4,0))
layout(matrix(c(2,1), 1, 2), widths = c(1.9,1))
stripchart(xlim = c(0.8,1.7),main = snp.id, at = c(1,1.5),axes = F, ylim = c(-2,12),list(log2(d[d$SNP_VAR == 0,cDNA.K562.norm.ipcr]+0.5), log2(d[d$SNP_VAR == 1,cDNA.K562.norm.ipcr]+0.5)), vertical = TRUE, method = "jitter", pch = 19, cex =1, col = c("gray", "black"))
axis(2, at=log2(c(0.5,1.5,5.5,50.5,1000.5)), labels=c('N.D.',1,5,50,1000),cex.axis=1.5,las = 2)
axis(1, at= c(1,1.5), labels = c("C", "CA"), cex.axis  =1.5)     
mean.ref <- mean(d[d$SNP_VAR == 0,cDNA.K562.norm.ipcr])
mean.alt <- mean(d[d$SNP_VAR == 1,cDNA.K562.norm.ipcr])
lines(x=c(0.85,1.15), y = c(log2(mean.ref)+0.5,log2(mean.ref)+0.5), col = "red", lwd = 5)
lines(x=c(1.35,1.65), y = c(log2(mean.alt)+0.5,log2(mean.alt)+0.5), col = "red", lwd = 5)
#title(ylab = 'Normalized expression', cex.lab = 1.5, line = 4)

dev.off()

## It would be very nice to generate this figure with the genomic sites included. 
# For each replicate: 3 files (pat, mat, eq)
# I reran the analysis (with the start end positions) for chrom 15 ONLY. Results may differ as downsampling is random
# but it is the best i can do right now

rep1 <- readRDS(file = "tmp/example.chr15/42_1combined.sure..counts.snp_indel.15.RDS")
rep2 <- readRDS(file = "tmp/example.chr15/42_2combined.sure..counts.snp_indel.15.RDS")
rep3 <- readRDS(file = "tmp/example.chr15/43_1combined.sure..counts.snp_indel.15.RDS")
rep4 <- readRDS(file = "tmp/example.chr15/43_2combined.sure..counts.snp_indel.15.RDS")
rep5 <- readRDS(file = "tmp/example.chr15/44_1combined.sure..counts.snp_indel.15.RDS")
rep6 <- readRDS(file = "tmp/example.chr15/44_2combined.sure..counts.snp_indel.15.RDS")
rep7 <- readRDS(file = "tmp/example.chr15/45_1combined.sure..counts.snp_indel.15.RDS")
rep8 <- readRDS(file = "tmp/example.chr15/45_2combined.sure..counts.snp_indel.15.RDS")

rep1 <- rep1[grepl(snp.id, rep1$SNP_ID),]
rep2 <- rep2[grepl(snp.id, rep2$SNP_ID),]
rep3 <- rep3[grepl(snp.id, rep3$SNP_ID),]
rep4 <- rep4[grepl(snp.id, rep4$SNP_ID),]
rep5 <- rep5[grepl(snp.id, rep5$SNP_ID),]
rep6 <- rep6[grepl(snp.id, rep6$SNP_ID),]
rep7 <- rep7[grepl(snp.id, rep7$SNP_ID),]
rep8 <- rep8[grepl(snp.id, rep8$SNP_ID),]

all.reps <- rbind(rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep8)

#generate a random order to prevent overplotting
all.reps.sampled <- sample_n(all.reps, size = nrow(all.reps))


#add jitter to ease plotting
jitter <- sample(seq(from=-0.50,to=0.50,by=0.001),nrow(all.reps.sampled),replace=TRUE)


# the SNP_ABS_POS_hg19 is not the same as the SNP_POS_ABS, therefore I have to 
# calculate the distance from this point


png("data/processed/Figures/Indel_Example/Fig.1b.sure.signal.per.allele.coverage.png")
par(mar=c(4,6,4,0))
plot(ann = F, cex.lab = 1.5, axes = F, x=1, y=2, ylim = c(-2,12),  main = snp.id, xlim = c(97034658-400,97034658+400), cex.axis = 1.5 )
arrows(x0 = all.reps.sampled[all.reps.sampled$SNP_VAR == 0,]$start_hg19, x1 = all.reps.sampled[all.reps.sampled$SNP_VAR == 0,]$end_hg19, 
       y0 = log2(all.reps.sampled[all.reps.sampled$SNP_VAR == 0,]$cDNA.K562.norm.ipcr+0.5)+jitter[1:249], lwd = 2, col = "gray", length = 0)
arrows(x0 = all.reps.sampled[all.reps.sampled$SNP_VAR == 1,]$start_hg19, x1 = all.reps.sampled[all.reps.sampled$SNP_VAR == 1,]$end_hg19, 
       y0 = log2(all.reps.sampled[all.reps.sampled$SNP_VAR == 1,]$cDNA.K562.norm.ipcr+0.5)+jitter[250:289], lwd = 2, col = 1, length = 0)
axis(2, at=log2(c(0.5,1.5,5.5,50.5,1000.5)), labels=c('N.D.',1,5,50,1000),cex.axis=1.5,las = 2)
axis(1, at = 97034658+c(-400,-200,0,200,400), labels = c("-400", "-200", "0", "200", "400"), cex.axis  =1.5)
legend("topleft", legend = c(paste0("C (mean = ", round(mean(all.reps.sampled[all.reps.sampled$SNP_VAR == 0, "cDNA.K562.norm.ipcr"]), digits = 1), " ; N = ", nrow(all.reps[all.reps$SNP_VAR == 0,]),")"),
                         paste0("CA (mean = 9.0 ; N = ", nrow(all.reps[all.reps$SNP_VAR == 1,]),")")), fill = c("gray", 1), bty = "n")
title(ylab = 'Normalized expression', cex.lab = 1.5, line = 4)
title(xlab = "Chr15:97,034,658 +/- 400 bp", cex.lab = 1.5)
title(main = snp.id)
text(x = 97034658+170, y = 11.6, labels = expression("P < 1.5 x 10"^-8), col = 1)
text(x = 97034658+50, y = 11.5, labels = "]", cex = 1.3)
dev.off()



# Next thing I want to combine the two plots into one First do this, then go through the plots

png("data/processed/Figures/Indel_Example/Fig.1c.sure.signal.per.allele.combined.png")
dev.off()

wilcox.test(d[d$SNP_VAR == 0]$cDNA.K562.norm.ipcr, d[d$SNP_VAR==1]$cDNA.K562.norm.ipcr)
