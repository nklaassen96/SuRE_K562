
# read the data of the "somatic" variants that have reads for REF and ALT

sure.results <- readRDS("data/interim/R_Objects/results_sure_new.rds")

table(sure.results$ref.element.count)

png("data/processed/Figures/Somatic_SNPs/Fig.1a.hist.ref.element.count.somatic.png")
h <- hist(sure.results$ref.element.count, breaks = c(0:50) , col = 1, xlim = c(0,10), xaxt = "n", xlab = "Reference element count", border = 0)
axis(side = 1, at = h$mids[1:10], labels = c(1:10))
dev.off()

png("data/processed/Figures/Somatic_SNPs/Fig.1b.hist.alt.element.count.somatic.png")
h <- hist(sure.results$alt.element.count, breaks = c(0:50) , col = 1, xlim = c(0,10), xaxt = "n", xlab = "Alternative element count", border = 0)
axis(side = 1, at = h$mids[1:10], labels = c(1:10))
dev.off()

png("data/processed/Figures/Somatic_SNPs/Fig.1c.SNV.coverage.somatic.png")
par(mar = c(5,5,4,1)+0.1)
h <- hist(sure.results$ref.element.count+sure.results$alt.element.count, breaks = c(0:50) , col = 1, xlim = c(0,10), xaxt = "n", xlab = "SNV coverage", border = 0, cex.axis = 1.5, cex.lab = 1.5)
axis(side = 1, at = h$mids[1:10], labels = c(1:10), cex.axis = 1.5, cex.lab = 1.5)
text(paste0("mean = ", round(mean(sure.results$ref.element.count+sure.results$alt.element.count), digits = 1)), x = 5, y = 3500, cex = 1.5)
text(paste0("median = ", round(median(sure.results$ref.element.count+sure.results$alt.element.count), digits = 1)), x = 5, y = 3100, cex = 1.5)
dev.off()

png("data/processed/Figures/Somatic_SNPs/Fig.1d.SNV.coverage.per.variant.somatic.png")
par(mar = c(5,5,4,1)+0.1)
h <- hist(c(sure.results$ref.element.count,sure.results$alt.element.count), breaks = c(0:50) , col = 1, xlim = c(0,10), xaxt = "n", xlab = "SNV coverage per variant", border = 0,cex.axis = 1.5, cex.lab = 1.5, )
axis(side = 1, at = h$mids[1:10], labels = c(1:10),cex.axis = 1.5, cex.lab = 1.5)
text(paste0("mean = ", round(mean(c(sure.results$ref.element.count,sure.results$alt.element.count)), digits = 1)), x = 4, y = 11000, cex = 1.5)
text(paste0("median = ", round(median(c(sure.results$ref.element.count,sure.results$alt.element.count)), digits = 1)), x = 4, y = 10000, cex = 1.5)
dev.off()


barplot(table(sure.results$chr), las = 2)

#doenst work, order of chromosomes is not good, do it again but as factor

chr.vector <- factor(sure.results$chr, levels = paste0("chr",c(1:22,"X")))
sure.results$chr.factor <- chr.vector

png("data/processed/Figures/Somatic_SNPs/Fig.2.mutations.per.chrom.png", width = 700)
par(mar = c(5,5,4,1)+0.1)
b <- barplot(table(sure.results$chr.factor), ylim = c(0,1000),las = 2, ylab = NULL, col = c("gray","black")[c(1,2,2,1,1,1,1,1,2,2,1,2,2,2,1,1,2,1,1,2,1,2,2)],  names.arg =  c(1:22,"X"), xlab = "Chromosome", cex.names = 1.5, cex.lab = 1.5, cex.axis = 1.5, density = c(20,NULL)[c(2,1,2,2,2,2,2,2,2,1,2,1,2,2,2,2,1,2,2,1,2,1,2)])
arrows(x0 = b[c(3,9,13,14,23)], y0 = 900, y1 = 800, lwd = 3)
mtext("Somatic mutations", side = 2, line = 3.5, cex = 1.5)
dev.off()

stripchart(sure.results$alt.cDNA.mean, jitter = 0.5, method = "jitter", pch = 19, vertical = TRUE)

jitter <- sample(seq(-0.03,0.03, length.out = 400), replace = TRUE, size = 9176)
jitter_h <- sample(seq(-0.2,0.2, length.out = 400), replace = TRUE, size = 9176)

png("data/processed/Figures/Somatic_SNPs/Fig.3.effect.size.png")
par(mar = c(5,5,4,1)+0.1)
plot(log2(sure.results$ref.cDNA.mean/sure.results$alt.cDNA.mean)+jitter_h, -log10(sure.results$wilcox.pvalue)+jitter, pch = 19, col = alpha(1, 0.1), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, xlab = "log2(mean SuRE signal REF/ALT)", ylab = "-log10(p-value)")
abline(h=-log10(0.05), lty = 3, lwd = 3)
dev.off()
# I want to distinguish between p=1 and p=NaN
nan.vector <- is.nan(sure.results$wilcox.pvalue)
sure.results[nan.vector,"wilcox.pvalue"] 

p.log.vector <- -log10(sure.results$wilcox.pvalue)
p.log.vector[nan.vector] <- -0.2

fc.log.vector <- log2(sure.results$ref.cDNA.mean/sure.results$alt.cDNA.mean)
fc.nan <- is.nan(fc.log.vector) | is.infinite(fc.log.vector)
fc.log.vector[fc.nan] <- 0

plot.size <- sum(is.finite(log2(sure.results$ref.cDNA.mean/sure.results$alt.cDNA.mean)) & is.finite(-log10(sure.results$wilcox.pvalue)))

plot(fc.log.vector+jitter_h, p.log.vector+jitter, pch = 19, col = alpha(1, 0.01), las = 1)
axis(side = 2, at = -0.2, labels = "ND", las = 1)

stripchart(fc.log.vector, vertical = T, pch = 19, method = "jitter", col = alpha(1,0.01))
