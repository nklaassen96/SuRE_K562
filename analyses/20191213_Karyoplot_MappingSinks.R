library(karyoploteR)

png("data/processed/Figures/Indel_Example/Fig.2.mapping.sinks.karyoplot.png")
par(mar=c(5,4,4,2)+0.1)
kp <- plotKaryotype(genome = "hg19", chromosomes = c("chr1", "chr10"), cytobands = c(), plot.type =1)
kpPlotMarkers(kp, chr = c("chr1", "chr10"), x = c(569882, 42399363), labels = c("rs527599589", "rs12761468"), text.orientation = "horizontal", label.margin = 5, cex = 1.4, ignore.chromosome.ends = TRUE)
dev.off()


#load the combined dataset
#to lazy to actually script

All.SNP.new.old.combined.datatable$total.elements <- All.SNP.new.old.combined.datatable$ref.element.count + All.SNP.new.old.combined.datatable$alt.element.count
All.SNP.new.old.combined.datatable$new.total.elements <- All.SNP.new.old.combined.datatable$new.ref.element.count + All.SNP.new.old.combined.datatable$new.alt.element.count
All.SNP.new.old.combined.datatable$dif.elements <- All.SNP.new.old.combined.datatable$total.elements- All.SNP.new.old.combined.datatable$new.total.elements

tst <- All.SNP.new.old.combined.datatable[order(All.SNP.new.old.combined.datatable$dif.elements),]
tst2 <- tst[1:178,]


dif <- sort(All.SNP.new.old.combined.datatable$total.elements- All.SNP.new.old.combined.datatable$new.total.elements)

tail(dif, 1000)
summary(dif)
sum(dif)
karyoploteR::plotKaryotype()
