###
#
# FIRST GENERATE A SMALL DATASET AFTER NORMALIZING
#
###

# create a selection from the normalized counts
n = 1000000
norm.sure <- sure.42_2.counts.indel.all[1:n,]

# assign random SNP_VAR values 
random.snpvar <- sample(rep(c(0,1,0), length.out = n))
norm.sure$SNP_VAR <- random.snpvar

counts <- norm.sure

###
#
# RUN THE PVALUE WILCOXON SCRIPT (input = "counts", output = "results.indel")
#
###

results.sure <- results.indel

p.real <- sort(-log10(results.sure$K562.wilcoxon.pvalue))
p.shuf <- sort(-log10(results.sure$K562.wilcoxon.pvalue.random))

plot(p.shuf, p.real)
abline(a = 0, b = 1)

raqtl.k562 <- na.omit(results.sure[results.sure$K562.wilcoxon.pvalue < 1.05,])
raqtl.hepg2 <- na.omit(results.sure[results.sure$HepG2.wilcoxon.pvalue < 1.05,])

mat <- matrix(nrow = 2, ncol = 2)
colnames(mat) <- c("promoter", "other")
rownames(mat) <- c("sign.", "not sign.")
mat[1,1] <- nrow(results.indel[results.indel$location.annotation == "promoter" & results.indel$K562.wilcoxon.pvalue <0.005,])
mat[2,1] <- nrow(results.indel[results.indel$location.annotation == "promoter" & results.indel$K562.wilcoxon.pvalue >= 0.005,])
mat[1,2] <- nrow(results.indel[results.indel$location.annotation != "promoter" & results.indel$K562.wilcoxon.pvalue <0.005,])
mat[2,2] <- nrow(results.indel[results.indel$location.annotation != "promoter" & results.indel$K562.wilcoxon.pvalue >=0.005,])
#volcano k562
fold.change.k562 <- log2(raqtl.k562$K562.cDNA.ref.mean/raqtl.k562$K562.cDNA.alt.mean)
select.finite <- is.finite(fold.change.k562)
fold.change.k562 <- fold.change.k562[select.finite]

log10pvalue.k562 <- -log10(raqtl.k562$K562.wilcoxon.pvalue)
log10pvalue.k562 <- log10pvalue.k562[select.finite]

dfk <- data.frame(fold.change = fold.change.k562, p.value = log10pvalue.k562)

ggplot(dfk, aes(x = fold.change, y = p.value)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") + 
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "gray") +
  ggtitle("raQTL", subtitle = "K562")



#volcano hepg2
fold.change.hepg2 <- log2(raqtl.hepg2$HepG2.cDNA.ref.mean/raqtl.hepg2$HepG2.cDNA.alt.mean)
select.finite <- is.finite(fold.change.hepg2)
fold.change.hepg2 <- fold.change.hepg2[select.finite]

log10pvalue.hepg2 <- -log10(raqtl.hepg2$HepG2.wilcoxon.pvalue)
log10pvalue.hepg2 <- log10pvalue.hepg2[select.finite]

#plot(fold.change.hepg2, log10pvalue.hepg2, xlab = "log2 fold change", ylab = "-log10 pvalue", cex = 0.2)
#plot(x = 2, y = 2, add = T)

dfh <- data.frame(fold.change = fold.change.hepg2, p.value = log10pvalue.hepg2)

ggplot(dfh, aes(x = fold.change, y = p.value)) +
  geom_point() +
  ggtitle("raQTL", subtitle = "HepG2") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "gray")

  #+#geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = -Inf, ymax = 1.5), fill = "red", alpha = 0.02) +
  #geom_area(x = c(-10,10), y = c(1.6, 1.6))
