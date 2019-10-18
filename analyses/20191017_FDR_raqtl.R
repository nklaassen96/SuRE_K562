library(tidyverse)

# Load the dataframe with the p-values per variant
file.name <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_pvalue/sure.indel.dataframe.pvalue.all.RDS"
var.df <- readRDS(file.name)

# add columns for min/max coverage and min/max expression
var.df$max.elements <- apply(var.df[,c("ref.element.count","alt.element.count")], 1, max)
var.df$min.elements <- apply(var.df[,c("ref.element.count","alt.element.count")], 1, min)
var.df$max.k562.expression <- apply(var.df[,c("K562.cDNA.alt.mean", "K562.cDNA.ref.mean")], 1, max)
var.df$max.hepg2.expression <- apply(var.df[,c("K562.cDNA.alt.mean", "K562.cDNA.ref.mean")], 1, max)


#select for >4 elements?
var.df <- var.df[var.df$ref.element.count > 4 & var.df$alt.element.count > 4]

fig.output.dir <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/"

# Figure 1a. qq plot K562
png(filename = paste0(fig.output.dir, "qqplot_k562.png"))
p.real <- sort(-log10(var.df$K562.wilcoxon.pvalue))
p.shuf <- sort(-log10(var.df$K562.wilcoxon.pvalue.random))
plot(p.shuf, p.real, cex = 0.3, xlab = "-log10(p-values shuffled variants)", ylab = "-log10(p-values real variants)", main = "qq plot K562")
abline(a = 0, b = 1)
dev.off()

# Figure 1b. qq plot HepG2
png(filename = paste0(fig.output.dir, "qqplot_hepg2.png"))
p.real <- sort(-log10(var.df$HepG2.wilcoxon.pvalue))
p.shuf <- sort(-log10(var.df$HepG2.wilcoxon.pvalue.random))
plot(p.shuf, p.real, cex = 0.3, xlab = "-log10(p-values shuffled variants)", ylab = "-log10(p-values real variants)", main = "qq plot HepG2")
abline(a = 0, b = 1)
dev.off()






summary(p.real)

hist(p.shuf, ylim = c(0,100), xlim=c(0,9), col = "red", breaks = 20)
hist(p.real, ylim = c(0,100), xlim=c(0,9), col = alpha("gray",0.5), breaks = 40, add = T)



# for a random set of p values (10^-2 till 10^-10): calculate the FDR

FDR.values <- NULL
p.values <- NULL
sign.real.values <- NULL
sign.shuf.values <- NULL


for (p in 10^-seq(2,10,by=0.2)) {
  
  p.threshold <- -log10(p)
  sign.real <- sum(p.real > p.threshold)
  sign.shuf <- sum(p.shuf > p.threshold)
  
  FDR <- round(sign.shuf / (sign.real + sign.shuf), digits = 2)
  
  FDR.values <- c(FDR.values, FDR)
  p.values <- c(p.values, p)
  sign.real.values <- c(sign.real.values, sign.real)
  sign.shuf.values <- c(sign.shuf.values, sign.shuf)
    
  }
df.hepg2 <- data.frame(p.values, sign.real.values, sign.shuf.values, FDR.values)
colnames(df.hepg2) <- c("p value threshold", "nr. of significant 'real' p values", "nr. of significant shuffled p values", "FDR")

plot(p.values, FDR.values)

# Figure 1x.  
as.data.frame(table(var.df))


