# This is a short script to generate overlapping VennDiagrams for the overlapping
# raqtls for K562 in the old and new pipeline by comparing the downsampling

# Load required libraries
library(VennDiagram)
library(data.table)
library(tidyverse)

# load all datasets

combined.new.old <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/All.SNP.new.old.combined.datatable.RDS")
all.new <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")



# First load the specific raQTL sets with the old downsampling vectors

new.raqtl.k562 <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
new.raqtl.k562 <- new.raqtl.k562[new.raqtl.k562$snp.type == "snp",]

old.raqtl.k562 <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/published.raqtl.k562.RDS")




# Determine the overlapping dataframe of indels

overlap.raqtl <- new.raqtl.k562[which(new.raqtl.k562$SNP_ID %in% old.raqtl.k562$SNP_ID),]

png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Comparison_analysis/Venn.K562.raQTL.overlap.png")
draw.pairwise.venn(area1 = nrow(new.raqtl.k562),
                   area2 = nrow(old.raqtl.k562),
                   cross.area = nrow(overlap.raqtl),
                   category = c("New K562 raQTLs (SNPs only)", "Published K562 raQTLs"), 
                   fill = c(2,3), 
                   cat.pos = c(-20,20),cex = rep(2,3),cat.cex = rep(1.3,2))
dev.off()

new.raqtl.old.noraqtl.idx <- which(!new.raqtl.k562$SNP_ID %in% old.raqtl.k562$SNP_ID)
old.raqtl.new.noraqtl.idx <- which(!old.raqtl.k562$SNP_ID %in% new.raqtl.k562$SNP_ID)


new.raqtl.old.noraqtl <- new.raqtl.k562[new.raqtl.old.noraqtl.idx,]
old.raqtl.new.noraqtl <- old.raqtl.k562[old.raqtl.new.noraqtl.idx,]

new.noraqtl.old.raqtl.idx <- which(all.new$SNP_ID %in% old.raqtl.new.noraqtl$SNP_ID)
new.noraqtl.old.raqtl <- all.new[new.noraqtl.old.raqtl.idx,]

# write snp.ids to .txt files for ludo for further analysis

write(new.raqtl.old.noraqtl$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/tmp/new.raqtl.k562.snp.ids.txt")
write(old.raqtl.new.noraqtl$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/tmp/old.raqtl.k562.snp.ids.txt")
write(overlap.raqtl$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/tmp/overlap.raqtl.k562.snp.ids.txt")

# I want to plot the reference element counts for the new vs old pipeline in a scatterplot for 3 distinct groups:
# (1) only new raqtl, (2) only old raqtl and (3) overlapping raqtl. 
# Create these 3 distinct dataframes with data from both pipelines

combined.new <- combined.new.old[which(combined.new.old$SNP_ID %in% new.raqtl.old.noraqtl$SNP_ID),]
combined.old <- combined.new.old[which(combined.new.old$SNP_ID %in% old.raqtl.new.noraqtl$SNP_ID),]
combined.overlap <- combined.new.old[which(combined.new.old$SNP_ID %in% overlap.raqtl$SNP_ID),]

# plot new
plot(x = combined.new[[4]], y = combined.new[[22]], cex = 0.1, xlab = "Published SNPs", ylab = "Novel SNPs", 
     main = paste0("All SNPs (n=", nrow(combined.new), "): ref.element.count"),
     ylim = c(0,1000), col = alpha(1,0.05))
text(x = 500, y = 1000 , labels = paste0("sampled to n=",length(sample.vector)))
abline(a = 0, b = 1)

# plot new
plot(x = combined.old[[4]], y = combined.old[[22]], cex = 0.1, xlab = "Published SNPs", ylab = "Novel SNPs", 
     main = paste0("All SNPs (n=", nrow(combined.new), "): ref.element.count"),
     ylim = c(0,1000), col = alpha(1,0.05))
text(x = 500, y = 1000 , labels = paste0("sampled to n=",length(sample.vector)))
abline(a = 0, b = 1)

# plot overlap
plot(x = combined.overlap[[4]], y = combined.overlap[[22]], cex = 0.1, xlab = "Published SNPs", ylab = "Novel SNPs", 
     main = paste0("All SNPs (n=", nrow(combined.new), "): ref.element.count"),
     ylim = c(0,1000), col = alpha(1,0.05))
text(x = 500, y = 1000 , labels = paste0("sampled to n=",length(sample.vector)))
abline(a = 0, b = 1)

# This type of plotting does notreally work, I want different colors for each plot. Therefore I have to combine
# the dataframe and add an additional column specifying the subset it came from

combined.new$subset <- factor("new")
combined.old$subset <- factor("old")
combined.overlap$subset <- factor("overlap")

#Make a combined datatable, but sample to the lowest number to make differences better visualizable

combined.all <- rbind(combined.old, 
                      sample_n(combined.new, size = nrow(combined.old)),
                      sample_n(combined.overlap, size = nrow(combined.old)))

combined.all <- rbind(combined.old, 
                      combined.new,
                      combined.overlap)

#plot ref.element.count per subset
{

png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Comparison_analysis/raQTL.K562.persubset.ref.element.count.png")
plot(x = combined.all$ref.element.count,
     y = combined.all$new.ref.element.count,
     col = c(alpha(2, 0.3),
             alpha(3, 0.3),
             alpha(4, 0.3))[combined.all$subset], 
     ylim = c(0,500),
     xlim = c(0,500),
     cex = 0.2,
     main = paste0("K562 raQTLs (sampled to n=",nrow(combined.old), ")"), 
     xlab = "Old ref.element.count",
     ylab = "New ref.element.count")
abline(a = 0, b = 1)
legend("topright", 
       fill = c(alpha(2, 1),
                alpha(3, 1),
                alpha(4, 1)),
       legend = c(paste0("Only old raQTL (n=",sum(combined.all$subset == "old"), ")"),
                  paste0("Only new raQTL (sampled to n=",sum(combined.all$subset == "new"), ")"),
                  paste0("Both old and new raQTL (sampled to n=",sum(combined.all$subset == "overlap"), ")")))
dev.off()

# I notice that increased ref. reads are enriched for new snps
}

#plot alt.element.count per subset
{
png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Comparison_analysis/raQTL.K562.persubset.alt.element.count.png")
plot(x = combined.all$alt.element.count,
     y = combined.all$new.alt.element.count,
     col = c(alpha(2, 0.3),
             alpha(3, 0.3),
             alpha(4, 0.3))[combined.all$subset], 
     ylim = c(0,500),
     xlim = c(0,500),
     cex = 0.2,
     main = paste0("K562 raQTLs (sampled to n=",nrow(combined.old), ")"), 
     xlab = "Old alt.element.count",
     ylab = "New alt.element.count")
abline(a = 0, b = 1)
legend("topright", 
       fill = c(alpha(2, 1),
                alpha(3, 1),
                alpha(4, 1)),
       legend = c(paste0("Only old raQTL (n=",sum(combined.all$subset == "old"), ")"),
                  paste0("Only new raQTL (sampled to n=",sum(combined.all$subset == "new"), ")"),
                  paste0("Both old and new raQTL (sampled to n=",sum(combined.all$subset == "overlap"), ")")))
dev.off()

# I notice that increased alt. reads are enrichted for new snps
}

#plot alt + ref.element.count per subset
{
png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Comparison_analysis/raQTL.K562.persubset.alt_and_ref.element.count.png")
plot(x = combined.all$alt.element.count+combined.all$ref.element.count,
     y = combined.all$new.alt.element.count+combined.all$new.ref.element.count,
     col = c(alpha(2, 0.3),
             alpha(3, 0.3),
             alpha(4, 0.3))[combined.all$subset], 
     ylim = c(0,600),
     xlim = c(0,600),
     cex = 0.2,
     main = paste0("K562 raQTLs (sampled to n=",nrow(combined.old), ")"), 
     xlab = "Old alt+ref.element.count",
     ylab = "New alt+ref.element.count")
abline(a = 0, b = 1)
legend("topright", 
       fill = c(alpha(2, 1),
                alpha(3, 1),
                alpha(4, 1)),
       legend = c(paste0("Only old raQTL (n=",sum(combined.all$subset == "old"), ")"),
                  paste0("Only new raQTL (sampled to n=",sum(combined.all$subset == "new"), ")"),
                  paste0("Both old and new raQTL (sampled to n=",sum(combined.all$subset == "overlap"), ")")))
dev.off()
}


# It seems that the new raqtls have gained the most reads. 
# It might be insightfull to generate histograms of the ref.element.difference

summary(combined.overlap$new.ref.element.count-combined.overlap$ref.element.count)
hist(combined.overlap$new.ref.element.count-combined.overlap$ref.element.count, breaks = 1000)

summary(combined.new$new.ref.element.count-combined.new$ref.element.count)
png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Comparison_analysis/raQTL.K562.boxplot.ref.element.count.difference.png")
boxplot(combined.old$new.ref.element.count-combined.old$ref.element.count,
        combined.new$new.ref.element.count-combined.new$ref.element.count, 
        combined.overlap$new.ref.element.count-combined.overlap$ref.element.count,names = c("Old raQTL", "New raQTL", "Overlaping raQTL"),col = c(2,3,4),outline = FALSE, ylab = "new.ref.element.count - old.ref.element.count")       
abline(h = 0, cex = 2, lty = 3)
dev.off()

## For the SNPs that I do not find (old.raqtl.new.noraqtl) signicant. I want to know where their p-values exactly are
# therefore i want to plot all new p-values (downsampled to the amount of unfindable snps) and only the "unfindable" p-values


all.sampled <- sample_n(all.new, size = nrow(new.noraqtl.old.raqtl))

png(filename = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Comparison_analysis/All.Hist.p.value.distribution.png")
hist(-log10(all.sampled$K562.wilcoxon.pvalue), breaks = 100, col = "gray", xlim = c(0,6), main = "P-value distribution new pipeline", xlab = "-log10(p-value)")
hist(-log10(new.noraqtl.old.raqtl$K562.wilcoxon.pvalue), breaks = 200, add = T, col = alpha(2,0.5))

# I want to add a vertical line that indicates the p-value threshold that was used in the NEW pipeline to define
# the raqtls. This p-value is the highest p-value within the completet K562 raQTL new set.

abline(v = -log10(max(new.raqtl.k562$K562.wilcoxon.pvalue)), lty = 2, cex = 3, lwd = 2)
legend("topright", c(paste0("All tested SNPs (sampled to n=", nrow(all.sampled), ")"), 
                     paste0("Old raQTL, new no raQTL (n=", nrow(new.noraqtl.old.raqtl), ")")), fill = c("gray", alpha(2,0.5)))
dev.off()
