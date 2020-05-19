library(karyoploteR)

png("data/processed/Figures/Indel_Example/Fig.2.mapping.sinks.karyoplot.png")
par(mar=c(5,4,4,2)+0.1)
kp <- plotKaryotype(genome = "hg19", chromosomes = c("chr1", "chr10"), cytobands = c(), plot.type =1)
kpPlotMarkers(kp, chr = c("chr1", "chr10"), x = c(569882, 42399363), labels = c("rs527599589", "rs12761468"), text.orientation = "horizontal", label.margin = 5, cex = 1.4, ignore.chromosome.ends = TRUE)
dev.off()


#load the combined dataset
#to lazy to actually script

All.SNP.new.old.combined.datatable <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/All.SNP.new.old.combined.datatable.RDS")

All.SNP.new.old.combined.datatable$total.elements <- All.SNP.new.old.combined.datatable$ref.element.count + All.SNP.new.old.combined.datatable$alt.element.count
All.SNP.new.old.combined.datatable$new.total.elements <- All.SNP.new.old.combined.datatable$new.ref.element.count + All.SNP.new.old.combined.datatable$new.alt.element.count
All.SNP.new.old.combined.datatable$dif.elements <- All.SNP.new.old.combined.datatable$total.elements- All.SNP.new.old.combined.datatable$new.total.elements

# Find all variants with >10.000 elements mapped to them in the NEW pipeline (that are the first 178 if you order them)
tst <- All.SNP.new.old.combined.datatable[order(All.SNP.new.old.combined.datatable$dif.elements),]
tst2 <- tst[1:178,]


dif <- sort(All.SNP.new.old.combined.datatable$total.elements- All.SNP.new.old.combined.datatable$new.total.elements)

tail(dif, 1000)
summary(dif)
sum(dif)
karyoploteR::plotKaryotype()


# Now I want to find all variants with >10.000 elements mapped to them in the old pipeline.

All.SNP.old.Joris <- readRDS("/DATA/usr/joris/git/SuRE_shared/Joris/SNP_project/all_snp_analysis_wilcox/sure.snp.df.strand.nonspecific.exact.wilcox.downsampled.JvA.180806.rds")

# check if all variants are indeed in there

summary(All.SNP.old.Joris$ref.element.count)
# Yes, they are

# I want to generate a combined dataframe using this dataset

All.Var.new <- readRDS("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.20191113.RDS")

# There are a few duplicates in All.SNP.old 

sum(duplicated(All.SNP.old.Joris$SNP_ID))

# lets look at these duplicates

duplicates <- All.SNP.old.Joris[duplicated(All.SNP.old.Joris$SNP_ID),]
table(duplicates$SNP_ID)

# most of these are actually SNP_ID == "." So there is no real worry. 

# Lets remove them 
All.SNP.old.Joris <- All.SNP.old.Joris[!duplicated(All.SNP.old.Joris$SNP_ID),]

# Are there also duplicates in the new analysis?

sum(duplicated(All.Var.new$SNP_ID))
duplicates.new <- All.Var.new[duplicated(All.Var.new$SNP_ID),]

# Yeah, there is one, with extremely high coverage 3.9 million. SNP_ID == ".". I dont know why....
# Lets remove the duplicate

All.Var.new <- All.Var.new[!duplicated(All.Var.new$SNP_ID),]

# Order the dataframes on SNP_IDs

All.Var.new <- All.Var.new[order(All.Var.new$SNP_ID),]
All.SNP.old.Joris <- All.SNP.old.Joris[order(All.SNP.old.Joris$SNP_ID),]


overlap.pub.idx <- which(All.SNP.old.Joris$SNP_ID %in% All.Var.new$SNP_ID)
overlap.nov.idx <- which(All.Var.new$SNP_ID %in% All.SNP.old.Joris$SNP_ID)

length(overlap.pub.idx) == length(overlap.nov.idx)

overlap.pub <- All.SNP.old.Joris[overlap.pub.idx,]
overlap.nov <- All.Var.new[overlap.nov.idx,]

all.comb.dt <- cbind(overlap.pub, overlap.nov)
colnames(all.comb.dt)[26:52] <- paste0("new.", colnames(all.comb.dt)[26:52])

sum(all.comb.dt[[2]] != all.comb.dt[[26]]) == 0 #check, should be true

# Add the total amount of elements

all.comb.dt$total.elements <-     all.comb.dt$ref.element.count     + all.comb.dt$alt.element.count
all.comb.dt$new.total.elements <- all.comb.dt$new.ref.element.count + all.comb.dt$new.alt.element.count
all.comb.dt$dif.elements <-       all.comb.dt$new.total.elements    - all.comb.dt$total.elements         

summary(all.comb.dt$total.elements)
summary(all.comb.dt$new.total.elements)
summary(all.comb.dt$dif.elements) # positive == more elements in the NEW pipeline. Negative == less elements in the NEW pipeline

all.comb.dt <- all.comb.dt[order(all.comb.dt$dif.elements),]

saveRDS(object = all.comb.dt, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/All.SNP.new.old.combined.full.coverage.dataframe.20200206.RDS")

all.comb.dt <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/All.SNP.new.old.combined.full.coverage.dataframe.20200206.RDS")


# I want to make a much smaller dataframe that contains>10.000 element gain or loss
all.comb.dt[10:18] <- NULL

sure.10000.diff <- all.comb.dt[abs(all.comb.dt$dif.elements) >= 10000,]
sure.10000.gain <- all.comb.dt[all.comb.dt$dif.elements >= 10000,]
sure.10000.loss <- all.comb.dt[all.comb.dt$dif.elements <= -10000,]

# Gaining regions
png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Indel_Example/Fig.3a.Hist.10000.gain.new.vs.old.elements.png")
hist(log10(sure.10000.gain$total.elements), breaks = 20, col = 2, xlim = c(0,7), ylim = c(0,160), main = paste0("Variants that gain >10.000 mapped elements (n=", nrow(sure.10000.gain), ")"), xlab = "log10(element.coverage)")
abline(v=3, lwd = 3, lty = 2)
hist(log10(sure.10000.gain$new.total.elements), col = alpha(3,0.5), add = T)
legend("topleft", legend = c("mapped elements (old)","mapped elements (new)"), fill = c(2, alpha(3,0.5)), bty = "n")
legend(x=-0.25, y=150, legend = "1000 elements", lty = 2, lwd = 3, bty= "n")
dev.off()

# Losing regions
png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Indel_Example/Fig.3b.Hist.10000.loss.new.vs.old.elements.png")
hist(log10(sure.10000.loss$total.elements),xlim = c(0,7), col  =2 , main = paste0("Variants that lose >10.000 mapped elements (n=", nrow(sure.10000.loss), ")"))
abline(v=3, lwd = 3, lty = 2)
hist(log10(sure.10000.loss$new.total.elements), add = T, col = alpha(3,0.5), breaks = 20)
legend("topleft", legend = c("nr. mapped elements (old)","nr. mapped elements (new)"), fill = c(2, alpha(3,0.5)))
dev.off()

#Karyoplot

# lets see if we can plot all the gaining SNPs on the karyoplot
pdf("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Indel_Example/Fig.3c.Karyoplot.10000.gain.locations.pdf")
png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Indel_Example/Fig.3c.Karyoplot.10000.gain.locations.png")
kp <- plotKaryotype(genome = "hg19", plot.type =1)
kpPlotMarkers(kp, chr = sure.10000.gain$chr, x = sure.10000.gain$SNPabspos, labels = "", text.orientation = "horizontal", label.margin = 0, cex = 1, ignore.chromosome.ends = TRUE)
dev.off()
# same for losing snps


png("/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Indel_Example/Fig.3d.Karyoplot.10000.loss.locations.png")
kp <- plotKaryotype(genome = "hg19", plot.type =1)
kpPlotMarkers(kp, chr = sure.10000.loss$chr, x = sure.10000.loss$SNPabspos, labels = "", text.orientation = "horizontal", label.margin = 0, cex = 1, ignore.chromosome.ends = TRUE)
dev.off()

