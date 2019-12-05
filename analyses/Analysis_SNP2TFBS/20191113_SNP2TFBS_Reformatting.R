# this script is used to compare all variants and the raqtl specific ones 
# in the SNP2TFBS library

library(data.table)
library(tidyverse)


## Data preparation

# First load the datasets of (1) all variants, (2) K562 raQTLs and (3) HepG2 raQTLs
all.variants <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")
raqtl.k562 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
raqtl.hepg2 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/HepG2.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")

# These datasets still contain snps and indels, but I am interested in the indels only
all.variants <- all.variants[all.variants$snp.type == "indel",]
raqtl.k562 <- raqtl.k562[raqtl.k562$snp.type == "indel",]
raqtl.hepg2 <- raqtl.hepg2[raqtl.hepg2$snp.type == "indel",]

# I Also want to do the analysis for the snps only
all.variants <- all.variants[all.variants$snp.type == "snp",]
raqtl.k562 <- raqtl.k562[raqtl.k562$snp.type == "snp",]
raqtl.hepg2 <- raqtl.hepg2[raqtl.hepg2$snp.type == "snp",]


## Use SNP2TFBS to discover motif disruption


#  

# 1. SNP2TFBS requires a .txt file with the SNP_IDs. Therefore I Write table to txt
write.table(x = all.variants$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SNP2TFBS/snp.new.id.all.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = raqtl.k562$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SNP2TFBS/snp.new.id.raqtl.562.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = raqtl.hepg2$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SNP2TFBS/snp.new.id.raqtl.hepg2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# These written tables can be transfered to the Windows PC with WinSCP

# 2. do this as input on https://ccg.epfl.ch/snp2tfbs/snpselect.php

# 3. download the snp2tfbs files

download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_43056.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.raqtl.k562.20191114.txt")
download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_44625.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.raqtl.hepg2.20191114.txt")
download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_44897.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.all.20191114.txt")

 #for SNPS:

download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_41221.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.snp.raqtl.k562.20191118.txt")
download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_42153.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.snp.raqtl.hepg2.20191118.txt")
download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_40204.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.snp.all.20191118.txt")

# 4. import the downloaded files as Robjects
tfbs.raqtl.k562 <-  fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.raqtl.k562.20191114.txt", select = c(7,1,2,4,5,6))
tfbs.raqtl.hepg2 <- fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.raqtl.hepg2.20191114.txt", select = c(7,1,2,4,5,6))
tfbs.all <-         fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.all.20191114.txt", select = c(7,1,2,4,5,6))

 #for SNPS:

tfbs.raqtl.k562 <-  fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.snp.raqtl.k562.20191118.txt", select = c(7,1,2,4,5,6))
tfbs.raqtl.hepg2 <- fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.snp.raqtl.hepg2.20191118.txt", select = c(7,1,2,4,5,6))
tfbs.all <-         fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.snp.all.20191118.txt", select = c(7,1,2,4,5,6))




# 5. Reformat the dataframe for a row per transcription factor for the above 3 files

for (tfbs.file in c("tfbs.raqtl.k562", "tfbs.raqtl.hepg2", "tfbs.all")){
  
  df <- get(tfbs.file)
  
  x <- unlist(strsplit(df$V6, split = ";"))
  match.idx <- seq(1,length(x), by= 3)
  tf.idx <- seq(2,length(x), by=3)
  diff.score.idx <- seq(3,length(x), by = 3)
  
  col.match <- as.numeric(str_remove_all(string = x[match.idx], pattern = "MATCH="))
  col.tf <- str_remove_all(x[tf.idx], pattern = "TF=")
  col.diff.score <- str_remove_all(x[diff.score.idx], pattern = "ScoreDiff=")
  
  df$V6 <- NULL
  df$m <- col.match
  df$t <- col.tf
  df$d <- col.diff.score
  colnames(df) <- c("snp.id", "chrom","snp2tfbs.pos", "ref", "alt", "tfbs.match", "transcription.factor", "score.differences")
  df <- separate_rows(df, transcription.factor, score.differences, sep = ",")
  class(df$score.differences) <- "numeric"
  
  assign(x = tfbs.file, value = df)  
}
# 6. perform fisher exact test for enrichtment of motif alteration in raQTLs

for (tfbs.raqtl in c("tfbs.raqtl.k562","tfbs.raqtl.hepg2")){
  
  raqtl.motifaltering <- length(unique(get(tfbs.raqtl)$snp.id))
  nonraqtl.motifaltering <- length(unique(tfbs.all$snp.id)) - length(unique(get(tfbs.raqtl)$snp.id))
  
  if (tfbs.raqtl == "tfbs.raqtl.k562"){
    raqtl.nonmotifaltering <- length(unique(raqtl.k562$SNP_ID)) - length(unique(get(tfbs.raqtl)$snp.id))} else {
      raqtl.nonmotifaltering <- length(unique(raqtl.hepg2$SNP_ID)) - length(unique(get(tfbs.raqtl)$snp.id))}
  
  nonmotifaltering <- length(unique(all.variants$SNP_ID)) - length(unique(tfbs.all$snp.id))
  
  nonraqtl.nonmotifaltering <- nonmotifaltering - raqtl.nonmotifaltering
  
  mat <- matrix(data = c(raqtl.motifaltering,
                         nonraqtl.motifaltering,
                         raqtl.nonmotifaltering,
                         nonraqtl.nonmotifaltering), ncol = 2, nrow = 2)
  colnames(mat) <- c("motif.altering", "non-motif.altering")
  rownames(mat) <- c("raqtl", "non-raqtl")
  
  if (tfbs.raqtl == "tfbs.raqtl.k562"){enrichment.k562.pvalue <- fisher.test(mat)$p.value
  } else {enrichment.hepg2.pvalue <- fisher.test(mat)$p.value}
}

# 7. CONCORDANCE: Find the maximum scoredifferences and plug them into the raQTL / total INDEL dataframes

### For all variants

# For all variants it is a little bit more difficult. Because the concordance is measured to the SuRE signal
# which is determined for two difference cell lines. I will determine both concordances and add them together
# to determine the total concordance.
{
all.max.scorediff <- tapply(tfbs.all$score.differences, tfbs.all$snp.id, function(x){x[which.max(abs(x))]})

match.idx <- match(names(all.max.scorediff), all.variants$SNP_ID)
all.variants[match.idx,"all.max.scorediff"] <- all.max.scorediff

# First the K562 concordance

concordance.idx.1 <- which(all.variants$k562.max == "ref" & sign(all.variants$all.max.scorediff)==-1)
concordance.idx.2 <- which(all.variants$k562.max == "alt" & sign(all.variants$all.max.scorediff)==1)
non.concordance.idx.1 <- which(all.variants$k562.max == "ref" & sign(all.variants$all.max.scorediff)==1)
non.concordance.idx.2 <- which(all.variants$k562.max == "alt" & sign(all.variants$all.max.scorediff)==-1)
undetermined.idx <- which(sign(all.variants$all.max.scorediff) == 0) # maybe I should, alternatively classify this as concordant or non.concordant

all.variants[c(concordance.idx.1, concordance.idx.2),"motif.concordance.k562"] <- "concordance"
all.variants[c(non.concordance.idx.1, non.concordance.idx.2),"motif.concordance.k562"] <- "non.concordance"
all.variants[c(undetermined.idx), "motif.concordance.k562"] <- "undetermined"
all.variants[is.na(all.variants$motif.concordance.k562), "motif.concordance.k562"] <- "no motifs"

# Then for HepG2 concordance

concordance.idx.1 <- which(all.variants$hepg2.max == "ref" & sign(all.variants$all.max.scorediff)==-1)
concordance.idx.2 <- which(all.variants$hepg2.max == "alt" & sign(all.variants$all.max.scorediff)==1)
non.concordance.idx.1 <- which(all.variants$hepg2.max == "ref" & sign(all.variants$all.max.scorediff)==1)
non.concordance.idx.2 <- which(all.variants$hepg2.max == "alt" & sign(all.variants$all.max.scorediff)==-1)
undetermined.idx <- which(sign(all.variants$all.max.scorediff) == 0) # maybe I should, alternatively classify this as concordant or non.concordant

all.variants[c(concordance.idx.1, concordance.idx.2),"motif.concordance.hepg2"] <- "concordance"
all.variants[c(non.concordance.idx.1, non.concordance.idx.2),"motif.concordance.hepg2"] <- "non.concordance"
all.variants[c(undetermined.idx), "motif.concordance.hepg2"] <- "undetermined"
all.variants[is.na(all.variants$motif.concordance.hepg2), "motif.concordance.hepg2"] <- "no motifs"

}

### For K562 raQTL (1st is max value, 2nd is mean value)
{
# For every variant there might be multiple
k562.max.scorediff <- tapply(tfbs.raqtl.k562$score.differences, tfbs.raqtl.k562$snp.id, function(x){x[which.max(abs(x))]})

# this vector contains row# for which the scoredifference is.
match.idx <- match(names(k562.max.scorediff), raqtl.k562$SNP_ID)
raqtl.k562[match.idx,"k562.max.scorediff"] <- k562.max.scorediff

concordance.idx.1 <- which(raqtl.k562$k562.max == "ref" & sign(raqtl.k562$k562.max.scorediff)==-1)
concordance.idx.2 <- which(raqtl.k562$k562.max == "alt" & sign(raqtl.k562$k562.max.scorediff)==1)
non.concordance.idx.1 <- which(raqtl.k562$k562.max == "ref" & sign(raqtl.k562$k562.max.scorediff)==1)
non.concordance.idx.2 <- which(raqtl.k562$k562.max == "alt" & sign(raqtl.k562$k562.max.scorediff)==-1)
undetermined.idx <- which(sign(raqtl.k562$k562.max.scorediff) == 0) # maybe I should, alternatively classify this as concordant or non.concordant

raqtl.k562[c(concordance.idx.1, concordance.idx.2),"motif.concordance"] <- "concordance"
raqtl.k562[c(non.concordance.idx.1, non.concordance.idx.2),"motif.concordance"] <- "non.concordance"
raqtl.k562[c(undetermined.idx), "motif.concordance"] <- "undetermined"
raqtl.k562[is.na(raqtl.k562$motif.concordance), "motif.concordance"] <- "no motifs"
}

### For HepG2 raQTL
{
hepg2.max.scorediff <- tapply(tfbs.raqtl.hepg2$score.differences, tfbs.raqtl.hepg2$snp.id, function(x){x[which.max(abs(x))]})

match.idx <- match(names(hepg2.max.scorediff), raqtl.hepg2$SNP_ID)
raqtl.hepg2[match.idx,"hepg2.max.scorediff"] <- hepg2.max.scorediff

concordance.idx.1 <- which(raqtl.hepg2$hepg2.max == "ref" & sign(raqtl.hepg2$hepg2.max.scorediff)==-1)
concordance.idx.2 <- which(raqtl.hepg2$hepg2.max == "alt" & sign(raqtl.hepg2$hepg2.max.scorediff)==1)
non.concordance.idx.1 <- which(raqtl.hepg2$hepg2.max == "ref" & sign(raqtl.hepg2$hepg2.max.scorediff)==1)
non.concordance.idx.2 <- which(raqtl.hepg2$hepg2.max == "alt" & sign(raqtl.hepg2$hepg2.max.scorediff)==-1)
undetermined.idx <- which(sign(raqtl.hepg2$hepg2.max.scorediff) == 0) # maybe I should, alternatively classify this as concordant or non.concordant

raqtl.hepg2[c(concordance.idx.1, concordance.idx.2),"motif.concordance"] <- "concordance"
raqtl.hepg2[c(non.concordance.idx.1, non.concordance.idx.2),"motif.concordance"] <- "non.concordance"
raqtl.hepg2[c(undetermined.idx), "motif.concordance"] <- "undetermined"
raqtl.hepg2[is.na(raqtl.hepg2$motif.concordance), "motif.concordance"] <- "no motifs"
}


### Start generating figures

col.k562 <- "steelblue"
col.hepg2 <- "magenta4"


## Figure 1. Fractions of motif altering SNPs for (1) all indels, (2) K562 raQTL, (3) HepG2 raQTL

tf.altering.all <- length(unique(tfbs.all$snp.id)) / length(unique(all.variants$SNP_ID))
tf.altering.raqtl.k562 <- length(unique(tfbs.raqtl.k562$snp.id)) / length(unique(raqtl.k562$SNP_ID))
tf.altering.raqtl.hepg2 <- length(unique(tfbs.raqtl.hepg2$snp.id)) / length(unique(raqtl.hepg2$SNP_ID))
values <- c(tf.altering.raqtl.hepg2, tf.altering.raqtl.k562, tf.altering.all)

# quickly also calculate the fold changes
values[1]/values[3] #1.80
values[2]/values[3] #1.46

png(filename= "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/SNP2TFBS/Fig.1.Fraction.indels.affecting.TFBS.png")
par(mar=c(4,10,20,1)+.1)
bp <- barplot(values, names.arg = c("HepG2 indel raQTLs ", "K562 indel raQTLs ", "All indels tested"),
        horiz = TRUE, 
        xlab = "Fraction of Indels affecting TFBS", 
        col = c(col.hepg2, col.k562, "black"), 
        las = 1, 
        cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,0.5), )
text(values[2]+0.02, y = bp[2,], labels = "*",cex = 2 )
text(values[1]+0.02, y = bp[1,], labels = "*",cex = 2)
dev.off()

## Figure 2. Fractions of motif altering SNPs for (1) all snps, (2) K562 raQTL snps, (3) HepG2 raQTL snps

tf.altering.all <- length(unique(tfbs.all$snp.id)) / length(unique(all.variants$SNP_ID))
tf.altering.raqtl.k562 <- length(unique(tfbs.raqtl.k562$snp.id)) / length(unique(raqtl.k562$SNP_ID))
tf.altering.raqtl.hepg2 <- length(unique(tfbs.raqtl.hepg2$snp.id)) / length(unique(raqtl.hepg2$SNP_ID))
values <- c(tf.altering.raqtl.hepg2, tf.altering.raqtl.k562, tf.altering.all)

#quicly calculate fold changes
values[1]/values[3]
values[2]/values[3]

png(filename= "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/SNP2TFBS/Fig.2.Fraction.snp.affecting.TFBS.png")
par(mar=c(4,10,20,1)+.1)
bp <- barplot(values, names.arg = c("HepG2 SNP raQTLs ", "K562 SNP raQTLs ", "All SNPs tested"),
              horiz = TRUE, 
              xlab = "Fraction of SNPs affecting TFBS", 
              col = c(col.hepg2, col.k562, "black"), 
              las = 1, 
              cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,0.5), )
text(values[2]+0.02, y = bp[2,], labels = "*",cex = 2 )
text(values[1]+0.02, y = bp[1,], labels = "*",cex = 2)
dev.off()

# Figure 3. Fractions of matif altering SNPs and Indels combined


png(filename= "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/SNP2TFBS/Fig.3.Fraction.snp_indel.affecting.TFBS.png")
par(mar=c(4,10,20,1)+.1)
bp <- barplot(values, names.arg = c("HepG2 raQTLs ", "K562 raQTLs ", "All variants tested"),
              horiz = TRUE, 
              xlab = "Fraction of variants affecting TFBS", 
              col = c(alpha(col.hepg2, 0.4), alpha(col.k562, 0.4), alpha("black", 0.4)), 
              las = 1, 
              cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,0.5), )
bp <- barplot(values.snp, 
              horiz = TRUE, 
              col = c(col.hepg2, col.k562, "black"), 
              las = 1, 
              cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,0.5), add = T, density = 40)
legend("topright", c("SNPs", "Indels"), density = c(40,0))
dev.off()

#Transparent colors:
black.trans <- mixcolor(alpha = 0.6, color1 = RGB(t(col2rgb("black"))), color2 = RGB(t(col2rgb("white"))))
black.trans <- "#999999"
k562.trans <- mixcolor(alpha = 0.6, color1 = RGB(t(col2rgb(col.k562))), color2 = RGB(t(col2rgb("white"))))
k562.trans <- "#B5CDE1"
hepg2.trans <- mixcolor(alpha = 0.6, color1 = RGB(t(col2rgb(col.hepg2))), color2 = RGB(t(col2rgb("white"))))
hepg2.trans <- "#D199D1"

# Figure 4a. Motif concordance SNPs (0 = no change)

all.variants.snp <- all.variants[all.variants$snp.type == "snp",]
raqtl.k562.snp <- raqtl.k562[raqtl.k562$snp.type == "snp",]
raqtl.hepg2.snp <- raqtl.hepg2[raqtl.hepg2$snp.type == "snp",]

conc.all.snp <- sum(c(all.variants.snp$motif.concordance.k562, all.variants.snp$motif.concordance.hepg2) %in% "concordance") / sum(c(all.variants.snp$motif.concordance.k562, all.variants.snp$motif.concordance.hepg2) %in% c("concordance", "non.concordance"))
conc.k562.snp  <- sum(raqtl.k562.snp$motif.concordance == "concordance")  / sum(raqtl.k562.snp$motif.concordance %in% c("concordance", "non.concordance"))
conc.hepg2.snp <- sum(raqtl.hepg2.snp$motif.concordance == "concordance") / sum(raqtl.hepg2.snp$motif.concordance %in% c("concordance", "non.concordance"))
values.4a <- c(conc.hepg2.snp, conc.k562.snp, conc.all.snp)

png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/SNP2TFBS/Fig.4a.Motif.concordance.snps.0.ignored.png")
par(mar=c(4,10,20,1)+.1)
bp <- barplot(values.4a, names.arg = c("HepG2 SNP raQTLs ", "K562 SNP raQTLs ", "All SNPs tested"),
              horiz = TRUE, 
              xlab = "Fraction concordance SuRE and motif disruption", 
              col = c(col.hepg2, col.k562, "black"), 
              las = 1, 
              cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,1))
barplot(add = T,rep(0.5,3),
        horiz = TRUE, 
        col = c(hepg2.trans, k562.trans, black.trans), 
        las = 1, axes = F ,
        cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,1), )
text(x = 0.25, y = bp[,1]-0.02, "Expected by chance", col = 1, cex = 1)
text(values.4a[1]+0.02, y = bp[1,], labels = "*",cex = 2)
text(values.4a[2]+0.02, y = bp[2,], labels = "*",cex = 2)
dev.off()

# Figure 4b. Motif concordance Indels (0 = deletion of motif)

all.variants.indel <- all.variants[all.variants$snp.type == "indel",]
raqtl.k562.indel <- raqtl.k562[raqtl.k562$snp.type == "indel",]
raqtl.hepg2.indel <- raqtl.hepg2[raqtl.hepg2$snp.type == "indel",]

conc.all.indel <- sum(c(all.variants.indel$motif.concordance.k562, all.variants.indel$motif.concordance.hepg2) %in% "concordance") / sum(c(all.variants.indel$motif.concordance.k562, all.variants.indel$motif.concordance.hepg2) %in% c("concordance", "non.concordance"))
conc.k562.indel  <- sum(raqtl.k562.indel$motif.concordance == "concordance")  / sum(raqtl.k562.indel$motif.concordance %in% c("concordance", "non.concordance"))
conc.hepg2.indel <- sum(raqtl.hepg2.indel$motif.concordance == "concordance") / sum(raqtl.hepg2.indel$motif.concordance %in% c("concordance", "non.concordance"))
values.4b <- c(conc.hepg2.indel, conc.k562.indel, conc.all.indel)

png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/SNP2TFBS/Fig.4b.Motif.concordance.indels.0.ignored.png")
par(mar=c(4,10,20,1)+.1)
bp <- barplot(values.4b, names.arg = c("HepG2 Indel raQTLs ", "K562 Indel raQTLs ", "All Indels tested"),
              horiz = TRUE, 
              xlab = "Fraction concordance SuRE and motif disruption", 
              col = c(col.hepg2, col.k562, "black"), 
              las = 1, 
              cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,1))
barplot(add = T,rep(0.5,3), names.arg = NULL,
        horiz = TRUE, 
        col = c(hepg2.trans, k562.trans, black.trans), 
        las = 1, axisnames = FALSE,axes = FALSE,
        cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,1), )
text(x = 0.25, y = bp[,1]-0.02, "Expected by chance", col = 1, cex = 1)
text(values.4b[1]+0.02, y = bp[1,], labels = "*",cex = 2)
text(values.4b[2]+0.02, y = bp[2,], labels = "*",cex = 2)
dev.off()


# Figure 5 concordance per fraction of raqtls 

## try for K562 SNPs

raqtl.k562.snp <- raqtl.k562[raqtl.k562$snp.type == "snp" & raqtl.k562$K562.wilcoxon.pvalue < 0.00001,]
conc.k562.snp  <- sum(raqtl.k562.snp$motif.concordance == "concordance")  / sum(raqtl.k562.snp$motif.concordance %in% c("concordance", "non.concordance"))
# this works for 1 value, now in a loop

p.values <- 10^(-seq(2,10, length.out = 20))
concordance.snp <- NULL

for (i in 1:length(p.values)){
  raqtl.k562.snp.idx <- raqtl.k562[raqtl.k562$snp.type == "snp" & raqtl.k562$K562.wilcoxon.pvalue < p.values[i],]
  concordance.snp[i]  <- sum(raqtl.k562.snp.idx$motif.concordance == "concordance")  / sum(raqtl.k562.snp.idx$motif.concordance %in% c("concordance", "non.concordance"))
}
png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/SNP2TFBS/Fig.5b.concordance.cutoff.snps.png")
plot(ylim = c(0.5, 1),-log10(p.values), concordance.snp, pch = 19, ylab = "Fraction concordance SuRE and motif disruption", xlab = "-log10(p-value cutoff)", main = "K562 raQTLs (SNPs)")
dev.off()



### STATISTICS ##

## SNPs

# K562

no.raqtl.k562 <- all.variants.snp[which(!all.variants.snp$SNP_ID %in% raqtl.k562.snp$SNP_ID),]
snp.k562 <- rbind(table(no.raqtl.k562$motif.concordance.k562)[c("concordance", "non.concordance")],table(raqtl.k562.snp$motif.concordance)[c("concordance", "non.concordance")] )
fisher.test(snp.k562)

# HepG2

no.raqtl.hepg2 <- all.variants.snp[which(!all.variants.snp$SNP_ID %in% raqtl.hepg2.snp$SNP_ID),]
snp.hepg2 <- rbind(table(no.raqtl.hepg2$motif.concordance.hepg2)[c("concordance", "non.concordance")],table(raqtl.hepg2.snp$motif.concordance)[c("concordance", "non.concordance")] )
fisher.test(snp.hepg2)

## Indels

# K562

no.raqtl.k562 <- all.variants.indel[which(!all.variants.indel$SNP_ID %in% raqtl.k562.indel$SNP_ID),]
indel.k562 <- rbind(table(no.raqtl.k562$motif.concordance.k562)[c("concordance", "non.concordance")],table(raqtl.k562.indel$motif.concordance)[c("concordance", "non.concordance")] )
fisher.test(indel.k562)

# HepG2

no.raqtl.hepg2 <- all.variants.indel[which(!all.variants.indel$SNP_ID %in% raqtl.hepg2.indel$SNP_ID),]
indel.hepg2 <- rbind(table(no.raqtl.hepg2$motif.concordance.hepg2)[c("concordance", "non.concordance")],table(raqtl.hepg2.indel$motif.concordance)[c("concordance", "non.concordance")] )
fisher.test(indel.hepg2)






