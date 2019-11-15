
library(VennDiagram)

# First load the datasets of (1) all variants, (2) K562 raQTLs and (3) HepG2 raQTLs
all.variants <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")
raqtl.k562 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
raqtl.hepg2 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/HepG2.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")

### Start generating figures

col.k562 <- "steelblue"
col.hepg2 <- "magenta4"


## Figure 1. Fractions of indels forr (1) all variants, (2) K562 raQTL, (3) HepG2 raQTL



indel.all <- sum(all.variants$snp.type == "indel")/nrow(all.variants)
indel.raqtl.k562 <- sum(raqtl.k562$snp.type == "indel")/nrow(raqtl.k562)
indel.raqtl.hepg2 <- sum(raqtl.hepg2$snp.type == "indel")/nrow(raqtl.hepg2)
values <- c(indel.raqtl.hepg2, indel.raqtl.k562, indel.all)

pdf(file= "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/raQTLs/Fig.1.Fraction.of.indels.pdf")
par(mar=c(4,10,20,1)+.1)
bp <- barplot(values, names.arg = c("HepG2 raQTLs ", "K562 raQTLs ", "All variants tested"),
              horiz = TRUE, 
              xlab = "Fraction of Indels", 
              col = c(col.hepg2, col.k562, "black"), 
              las = 1, 
              cex.axis = 1.5, cex.names = 1, width = c(0.5,0.5,0.5), xpd = FALSE, xlim = c(0,0.14), )
text(values[2]+0.005, y = bp[2,], labels = "*",cex = 2 )
text(values[1]+0.005, y = bp[1,], labels = "*",cex = 2)
dev.off()

# looks like there is a significant difference. Lets figure that out with a fischer exact test
# i need to have a contingency table

table(all.variants$snp.type)
table(raqtl.k562$snp.type)

mat.k562 <-  matrix(nrow = 2, ncol = 2, c(table(all.variants$snp.type), table(raqtl.k562$snp.type)), dimnames = list(c("indel", "snp"), c("all", "k562 raqtl")))
mat.hepg2 <- matrix(nrow = 2, ncol = 2, c(table(all.variants$snp.type), table(raqtl.hepg2$snp.type)), dimnames = list(c("indel", "snp"), c("all", "hepg2 raqtl")))

# Do they look good?
mat.k562
mat.hepg2
sum(mat.k562[,"all"]) == nrow(all.variants)

# Yes, looks fine. Now perform the fischer test

fisher.test(mat.hepg2)
fisher.test(mat.k562)

# Nice, very significant outcome (I think it is the max (2.2E-16))



## Figure 2. Venn diagram of overlapping Indels

all.indel <- all.variants[all.variants$snp.type == "indel",]
raqtl.hepg2.indel <- raqtl.hepg2[raqtl.hepg2$snp.type == "indel",]
raqtl.k562.indel <- raqtl.k562[raqtl.k562$snp.type == "indel",]

png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/raQTLs/Fig.2.Venn.raQTL.indel.overlap.png")
draw.pairwise.venn(area1 = nrow(raqtl.k562.indel), 
                   area2 = nrow(raqtl.hepg2.indel), 
                   cross.area = sum(raqtl.k562.indel$SNP_ID %in% raqtl.hepg2.indel$SNP_ID),
                   category = c("K562", "HepG2"), 
                   cat.pos = c(0,0),cex = 2,cat.cex = 2.5, fill = c(col.k562, col.hepg2), alpha = c(0.8,0.8), cat.fontfamily = c("sans","sans"), fontfamily = c("sans", "sans", "sans"))
dev.off()

## Figure 3. SuRE activity for all fragments of one specific indel

# I choose a indel for which the wilcox.value in K562 is very high and the alternative is lower than the 
# reference alele

selected.snp.id <- "rs373322167"

# As this indel is on chromosome 22, I have to load the data for that specific chromosome
reads.chr5 <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE_combined/sure.snp_indel.combined.chrom.5.RDS")


# unfortunately I have already thrown away the start and end positions of the gDNA fragments. 

reads.chr5 <- reads.chr5[reads.chr5$SNP_ID == selected.snp.id,]

#change the sure.signal = 0 to NA for easier plotting

reads.chr5[which(reads.chr5$cDNA.K562.norm.ipcr ==0),]$cDNA.K562.norm.ipcr <- NA

#change back if needed
reads.chr5[is.na(reads.chr5$cDNA.K562.norm.ipcr)]$cDNA.K562.norm.ipcr <- 0

sure.reference <- reads.chr5$cDNA.K562.norm.ipcr[reads.chr5$SNP_VAR == 0]
sure.alternative <- reads.chr5$cDNA.K562.norm.ipcr[reads.chr5$SNP_VAR == 1]

# maybe i can change al 0 to 1, that will be 0

#I can make the log values in the dataframe
reads.chr5$cDNA.K562.norm.ipcr <- log10(reads.chr5$cDNA.K562.norm.ipcr)

par(mar=c(1,100,2,8)+.1)
ggplot(reads.chr5, aes(x=SNP_SEQ,y = cDNA.K562.norm.ipcr)) +
  geom_point(position = "jitter") +

  theme_classic() +
  xlab(NULL) +
  ylab("Normalized SuRE signal") +
scale_y_continuous(na.value = 0)
log10(0.4
      )

