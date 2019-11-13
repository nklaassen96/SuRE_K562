# This script is used to compare the new and old pipeline with eachother.
# The new pipeline has a few significant differences. 


library(data.table)

published <- fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SuRE_OSF_NatGen/download.txt.gz")
published <- published[order(published$SNP_ID),]

colnames(published)[6] <- "K562.cDNA.ref.mean"
colnames(published)[7]<- "K562.cDNA.alt.mean"
colnames(published)[8] <- "K562.wilcoxon.pvalue"
colnames(published)[9] <- "K562.wilcoxon.pvalue.random"
colnames(published)[12] <- "HepG2.wilcoxon.pvalue"
colnames(published)[13] <- "HepG2.wilcoxon.pvalue.random"

##### recreate the exact same raQTL as in the published paper (so not with my for-loop) ####
df <- published
df <- df[df$K562.cDNA.ref.mean >4 | df$K562.cDNA.alt.mean >4,]
df$max.elements <- apply(df[,c("ref.element.count","alt.element.count")], 1, max)
df$min.elements <- apply(df[,c("ref.element.count","alt.element.count")], 1, min)
df <- df[df$max.elements <1000 & df$min.elements >= 10,]

published.raqtl.k562 <- df[df$K562.wilcoxon.pvalue <= 0.006192715,]

# create my own raqtl for k562 with the var.df script in FDR

novel.raqtl.k562 <- raqtl.k562

published.raqtl.novel.noraqtl.idx <-which(!published.raqtl.k562$SNP_ID %in% novel.raqtl.k562$SNP_ID)
published.raqtl.novel.norqtl <- published.raqtl.k562[published.raqtl.novel.noraqtl.idx,]
novel.noraqtl.published.raqtl.idx <- which(novel$SNP_ID %in% published.raqtl.novel.norqtl$SNP_ID)
novel.noraqtl.published.raqtl <- novel[novel.noraqtl.published.raqtl.idx,]

saveRDS(object = pub.raqtl.k, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/tmp/pub.raqtl.k.RDS")


# same thresholds as published (Sure >4) # use the other FDR script for this.
novel <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_pvalue/sure.snp_indel.dataframe.pvalue.all.20191105.RDS")
novel <- novel[novel$snp.type == "snp",]
novel <- novel[order(novel$SNP_ID),]



novel.raqtl.k <- raqtl.k562[order(raqtl.k562$SNP_ID),]
novel.raqtl.h <- raqtl.hepg2[order(raqtl.hepg2$SNP_ID),]
#pub.raqtl.k <- raqtl.k562[order(raqtl.k562$SNP_ID),]

# determine the overlap and calculate this overlap
# first order on snp.id


overlap.k_pub <- which(published$SNP_ID %in% novel.raqtl.k$SNP_ID)
overlap.k_nov <- which(novel.raqtl.k$SNP_ID %in% published$SNP_ID)
length(overlap.k_pub)/nrow(novel.raqtl.k)*100 #percentage of novel SNP_ID in raQTLs that are present in the published DB


pub <- published[overlap.k_pub,]
nov <- novel.raqtl.k[overlap.k_nov,]

    # generate all.comb.dt df
    published <- published[!duplicated(published$SNP_ID),]
    
    
    overlap.pub.idx <- which(published$SNP_ID %in% novel$SNP_ID)
    overlap.nov.idx <- which(novel$SNP_ID %in% published$SNP_ID)
    
    overlap.pub <- published[overlap.pub.idx,]
    overlap.nov <- novel[overlap.nov.idx,]
    
    all.comb.dt <- cbind(overlap.pub, overlap.nov)
    colnames(all.comb.dt)[16:41] <- paste0("new.", colnames(all.comb.dt)[16:41])
    
    sum(all.comb.dt[[2]] != all.comb.dt[[16]]) == 0 #check, should be true




# generate a figure directory

fig.dir <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Comparison_analysis/"

sample.vector <- sample(seq(1:nrow(all.comb.dt)), size = 200000)

# -log10(p)
png(filename = paste0(fig.dir, "All.pvalue.png"))
plot(x = -log10(all.comb.dt[[8]][sample.vector]), y = -log10(all.comb.dt[[32]][sample.vector]), col = alpha(1, 0.01), cex = 0.1, xlab = "Published SNPs", ylab = "Novel SNPs", 
     main = paste0("All SNPs (n=", nrow(all.comb.dt), "): -log10(p.value)"),
     ylim = c(0,10), xlim = c(0,10)) 
text(x = 5, y = 10 , labels = paste0("sampled to n=",length(sample.vector)))
abline(a = 0, b = 1)
dev.off()

# ref.element.count
png(filename = paste0(fig.dir, "All.ref_element_count.png"))
plot(x = all.comb.dt[[4]][sample.vector], y = all.comb.dt[[22]][sample.vector], cex = 0.1, xlab = "Published SNPs", ylab = "Novel SNPs", 
     main = paste0("All SNPs (n=", nrow(all.comb.dt), "): ref.element.count"),
     ylim = c(0,1000), col = alpha(1,0.01))
text(x = 500, y = 1000 , labels = paste0("sampled to n=",length(sample.vector)))
abline(a = 0, b = 1)
dev.off()

# alt.element.count
png(filename = paste0(fig.dir, "All.alt_element_count.png"))
plot(x = all.comb.dt[[5]][sample.vector], y = all.comb.dt[[23]][sample.vector], cex = 0.1, xlab = "Published SNPs", ylab = "Novel SNPs", 
     main = paste0("All SNPs (n=", nrow(all.comb.dt),")", ": alt.element.count"),
     ylim = c(0,1000),
     xlim = c(0,1000),
     col = alpha(1,0.01))
abline(a = 0, b = 1)
text(x = 500, y = 1000 , labels = paste0("sampled to n=",length(sample.vector)))
dev.off()


# K562.ref.mean
png(filename = paste0(fig.dir, "All.K562_ref_cDNA_mean.png"))
plot(x = all.comb.dt[[6]][sample.vector], y = all.comb.dt[[24]][sample.vector], cex = 0.1, xlab = "Published SNPs", ylab = "Novel SNPs", 
     main = paste0("All SNPs (n=", nrow(all.comb.dt),")", ": K562 ref.cDNA.mean"), col = alpha(1,0.01),
     xlim = c(0,10),
     ylim = c(0,10))
abline(a = 0, b = 1)
text(x = 5, y = 10 , labels = paste0("sampled to n=",length(sample.vector)))
dev.off()

# K562.alt.mean
png(filename = paste0(fig.dir, "All.K562_alt_cDNA_mean.png"))
plot(x = all.comb.dt[[7]][sample.vector], y = all.comb.dt[[25]][sample.vector], cex = 0.1, xlab = "Published SNPs", ylab = "Novel SNPs", 
     main = paste0("All SNPs (n=", nrow(all.comb.dt),")", ": K562 alt.cDNA.mean"), col = alpha(1,0.01),
     xlim = c(0,10),
     ylim = c(0,10))
text(x = 5, y = 10 , labels = paste0("sampled to n=",length(sample.vector)))
abline(a = 0, b = 1)
dev.off()

# I notice that there seems to be a "second slope" in the ref.cDNA.mean and alt.cDNA.mean
# I want to light out one specific example and see why this is the case

# First I select all snps that have a difference bigger than 5 in alt.cdna.mean
tst.idx <- which((all.comb.dt[[25]]-all.comb.dt[[7]]) > 5)
tst.df <- all.comb.dt[tst.idx,]
tst.df

sure.reads.allrep.chr.17 <- readRDS("~/projects/SuRE_K562/data/interim/SuRE_Test_Old_Data/sure.reads.allrep.chr.17.RDS")

#histogram novel p values

nov.raqtl_pub.NA.idx <- which(!novel.raqtl.k$SNP_ID %in% published$SNP_ID)
nov.raqtl_pub.NA <- novel.raqtl.k[nov.raqtl_pub.NA.idx,]

nov.raqtl_pub.not.raqtl.idx <- which(!novel.raqtl.k$SNP_ID %in% pub.raqtl.k$SNP_ID)
nov.raqtl_pub.not.raqtl <- novel.raqtl.k[nov.raqtl_pub.not.raqtl.idx,]

nov.raqtl_pub.raqtl.idx <- which(novel.raqtl.k$SNP_ID %in% pub.raqtl.k$SNP_ID)
nov.raqtl_pub.raqtl <- novel.raqtl.k[nov.raqtl_pub.raqtl.idx,]

# plot all nov.raqtl
#hist(-log10(novel.raqtl.k$K562.wilcoxon.pvalue), breaks= 3000, xlim = c(2,10), col = 2, main = "novel raQTL K562 (SNP only)", xlab = "-log10(p.value)")

# plot "novel raqtl but not a raqtl in published" only
hist(-log10(nov.raqtl_pub.not.raqtl$K562.wilcoxon.pvalue),xlim = c(2,10), breaks = 1000, add = F, col = 3,main = paste0("novel K562 raQTL  (SNP only) (n=", nrow(novel.raqtl.k),")"), xlab = "-log10(p.value)")

# plot "novel raqtl and also raqtl in published" only
hist(-log10(nov.raqtl_pub.raqtl$K562.wilcoxon.pvalue), breaks = 600, col = alpha(2,0.5), add = T)

# plot "novel raqtl but not present in published" only
hist(-log10(nov.raqtl_pub.NA$K562.wilcoxon.pvalue), breaks = 600, add = T, col = 6)

legend("topright", c(paste0("Novel raQTL (no raQTL in published) (n=", length(nov.raqtl_pub.not.raqtl.idx), ")"),
                     paste0("Novel raQTL (raQTL in published) (n=", length(nov.raqtl_pub.raqtl$K562.wilcoxon.pvalue), ")"),
                     paste0("Novel raQTL (SNP_ID not in published) (n=", length(nov.raqtl_pub.NA$K562.wilcoxon.pvalue), ")")
                     ),
       fill = c(3,alpha(2,0.5), 6))


#### TST #### check differences pub/nov

pub.raqtl_nov.not.raqtl.idx <- which(!pub.raqtl.k$SNP_ID %in% novel.raqtl.k$SNP_ID)
pub.raqtl_nov.raqtl.idx <-     which(pub.raqtl.k$SNP_ID %in% novel.raqtl.k$SNP_ID)

# check
length(pub.raqtl_nov.raqtl.idx)+length(pub.raqtl_nov.not.raqtl.idx) == nrow(pub.raqtl.k)

pub.raqtl_nov.not.raqtl <- pub.raqtl.k[pub.raqtl_nov.not.raqtl.idx,]
pub.raqtl_nov.raqtl <-     pub.raqtl.k[pub.raqtl_nov.raqtl.idx]

summary(-log10(pub.raqtl_nov.not.raqtl$K562.wilcoxon.pvalue))
summary(-log10(pub.raqtl_nov.raqtl$K562.wilcoxon.pvalue))

hist(-log10(pub.raqtl_nov.not.raqtl$K562.wilcoxon.pvalue), col = 3, breaks = 200, xlim = c(1.8,15), xlab = "-log10(p-value)", main = paste0("Published K562 raQTLs (n =", nrow(pub.raqtl.k), ")"))
hist(-log10(pub.raqtl_nov.raqtl$K562.wilcoxon.pvalue), add = T, col = alpha(2, 0.5), breaks = 400)
legend("topright", c(paste0("published raQTL (no novel raQTL) (n=", length(pub.raqtl_nov.not.raqtl$K562.wilcoxon.pvalue), ")"),
                     paste0("published raQTL (also novel raQTL) (n=", length(pub.raqtl_nov.raqtl$K562.wilcoxon.pvalue), ")")), fill = c(3,alpha(2,0.5)))

summary(pub.raqtl_nov.not.raqtl$ref.element.count)
summary(pub.raqtl_nov.raqtl$ref.element.count)

summary(pub.raqtl_nov.not.raqtl$alt.element.count)
summary(pub.raqtl_nov.raqtl$alt.element.count)



