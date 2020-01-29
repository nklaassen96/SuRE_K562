
library(data.table)
library(dplyr)


sure.new <- readRDS("data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")
sure.old <- fread("data/external/SuRE_OSF_NatGen/download.txt.gz")

sure.new.samp <- sample_n(sure.new, 300000)
sure.old.samp <- sample_n(sure.old, 300000)

# Fig. 8a
png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/General Statistics/Fig.8A.qqplot.new.old.K562.png")
par(mar=c(5,4.5,4,2)+0.1)
p.real.new <- sort(-log10(sure.new.samp$K562.wilcoxon.pvalue))
p.shuf.new <- sort(-log10(sure.new.samp$K562.wilcoxon.pvalue.random))
p.real.old <- sort(-log10(sure.old.samp$k562.wilcox.p.value))
p.shuf.old <- sort(-log10(sure.old.samp$k562.wilcox.p.value.random))
plot(cex.lab = 1.5, cex.axis = 1.5, bty = "l",c(p.shuf.new, p.shuf.old), c(p.real.new, p.real.old), pch = 19, cex = 0.3, xlab = "-log10(p-values shuffled variants)", ylab = "-log10(p-values real variants)", main = "qq plot K562", col = c(rep(1, 100000), rep(2,100000)))
abline(a = 0, b = 1)
legend("topleft", legend = c("New pipeline", "Old pipeline"), fill = c(1,2), cex = 1.5)
dev.off()


# Fig. 8b
png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/General Statistics/Fig.8B.qqplot.new.old.HepG2.png")
par(mar=c(5,4.5,4,2)+0.1)
p.real.new.h <- sort(-log10(sure.new.samp$HepG2.wilcoxon.pvalue))
p.shuf.new.h <- sort(-log10(sure.new.samp$HepG2.wilcoxon.pvalue.random))
p.real.old.h <- sort(-log10(sure.old.samp$hepg2.wilcox.p.value))
p.shuf.old.h <- sort(-log10(sure.old.samp$hepg2.wilcox.p.value.random))
plot(cex.lab = 1.5, cex.axis = 1.5, bty = "l",c(p.shuf.new.h, p.shuf.old.h), c(p.real.new.h, p.real.old.h), pch = 19, cex = 0.3, xlab = "-log10(p-values shuffled variants)", ylab = "-log10(p-values real variants)", main = "qq plot HepG2", col = c(rep(1, 100000), rep(2,100000)))
abline(a = 0, b = 1)
legend("topleft", legend = c("New pipeline", "Old pipeline"), fill = c(1,2), cex = 1.5)
dev.off()

# Fig 8c qqplot K562 indel vs snp

sure.new.snp.samp <- sample_n(sure.new[sure.new$snp.type == "snp",], 300000)
sure.new.indel.samp <- sample_n(sure.new[sure.new$snp.type == "indel",], 300000)

png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/General Statistics/Fig.8C.qqplot.snp.indel.k562.png")
par(mar=c(5,4.5,4,2)+0.1)
p.real.snp.k <- sort(-log10(sure.new.snp.samp$K562.wilcoxon.pvalue))
p.shuf.snp.k <- sort(-log10(sure.new.snp.samp$K562.wilcoxon.pvalue.random))
p.real.indel.k <- sort(-log10(sure.new.indel.samp$K562.wilcoxon.pvalue))
p.shuf.indel.k <- sort(-log10(sure.new.indel.samp$K562.wilcoxon.pvalue.random))
plot(cex.lab = 1.5, cex.axis = 1.5, bty = "l",c(p.shuf.snp.k, p.shuf.indel.k), c(p.real.snp.k, p.real.indel.k), pch = 19, cex = 0.3, xlab = "-log10(p-values shuffled variants)", ylab = "-log10(p-values real variants)", main = "qq plot K562", col = c(rep(1, length(p.real.indel.k)), rep(2,length(p.real.indel.k))))
abline(a = 0, b = 1)
legend("topleft", legend = c("SNP", "Indel"), fill = c(1,2), cex = 1.5)
dev.off()

png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/General Statistics/Fig.8D.qqplot.snp.indel.hepg2.png")
par(mar=c(5,4.5,4,2)+0.1)
p.real.snp.h <- sort(-log10(sure.new.snp.samp$HepG2.wilcoxon.pvalue))
p.shuf.snp.h <- sort(-log10(sure.new.snp.samp$HepG2.wilcoxon.pvalue.random))
p.real.indel.h <- sort(-log10(sure.new.indel.samp$HepG2.wilcoxon.pvalue))
p.shuf.indel.h <- sort(-log10(sure.new.indel.samp$HepG2.wilcoxon.pvalue.random))
plot(cex.lab = 1.5, cex.axis = 1.5, bty = "l",c(p.shuf.snp.h, p.shuf.indel.h), c(p.real.snp.h, p.real.indel.h), pch = 19, cex = 0.3, xlab = "-log10(p-values shuffled variants)", ylab = "-log10(p-values real variants)", main = "qq plot HepG2", col = c(rep(1, length(p.real.indel.k)), rep(2,length(p.real.indel.k))))
abline(a = 0, b = 1)
legend("topleft", legend = c("SNP", "Indel"), fill = c(1,2), cex = 1.5)
dev.off()

hist(sure.new.snp.samp$HepG2.wilcoxon.pvalue, col = "gray", main = "SNP")
hist(sure.new.snp.samp$HepG2.wilcoxon.pvalue.random, add = T, col = alpha(2, 0.5))

hist(sure.new.indel.samp$HepG2.wilcoxon.pvalue, col = "gray", main = "Indel")
hist(sure.new.indel.samp$HepG2.wilcoxon.pvalue.random, add = T, col = alpha(2, 0.5))

hist(sure.new.snp.samp$HepG2.wilcoxon.pvalue, col = "gray", main = "SNP vs. Indel", breaks= 100, xlab = "p-value")
hist(sure.new.indel.samp$HepG2.wilcoxon.pvalue, col = alpha(2, 0.5), main = "Indel", add = T, breaks = 100)
legend("topright", legend = c("snp", "indel"), fill = c("gray", alpha(2,0.5)), cex = 1.5)
