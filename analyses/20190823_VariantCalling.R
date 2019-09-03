## STARTUP
{
## Library installation and data reading
{
##Install Bioconductor
    #if (!requireNamespace("BiocManager", quietly = FALSE))
    #  install.packages("BiocManager")
    #BiocManager::install()

##install TxDB database
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene",INSTALL_opts = c('--no-lock'))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantTools)
library(karyoploteR)
library(ggplot2)
library(vcfR)
#library(VariantAnnotation) #is within VariantTools?
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#read the vcf file that was downloaded
vcfK <- readVcf("~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF606RIC.vcf.gz", 
               genome = "hg19")
}

## Testing Set START ##
{
## According to https://bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf ##

  fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")
  header(vcf) #shows the header
  samples(header(vcf)) #shows subset 'samples' of the header
  geno(header(vcf)) #shows subset 'geno' of the header
  head(rowRanges(vcf),5) #retrieve GRanges object of 3 rows
  ref(vcf)[1:5] #retrieve REF sequence for the first 5 rows in a DNAStringSet object
  qual(vcf)[1:5] #retrieve QUAL scores of the first 5 rows
  alt(vcf)[1:5] #retrieve ALT sequences (could be more than 1)
  geno(vcfK) #shows al genotypes as described in FORMAT fields
  geno(vcfK)$GT #retrieves genotypList
  
  sapply(geno(vcf), class)
  info(vcf) #shows some sort of summary of all the data. Rows are indicated by rs-numbers or chromosome positions 
  rowRanges(vcf) #Gives the Granges of all rows (what is the difference with head(Rowranges)???)
  as.data.frame(rowRanges(vcf)) #generates dataframe of these rows
  vcf[1:100,] #subset few rows within a VCF
  seqlevels(vcf) <- "chr22" #set the seqlevels all to "chr22" instead of 22 (this is the required format)
  rd <- rowRanges(vcf)
  #lengths(gr$ALT) #determine the amount of ALT sequences
  #table(width(gr)) #contingency table for the length of the DNA fragments (REF or ALT??), ik vermoed REF
  expand(vcf) #generates an "expanded-VCF" object that shows a row for every alternative sequence
  #gr_q <- gr[gr$QUAL>1500]
  #tapply(gr$SNV, seqnames(gr), sum) #retrieve nr. of SNVs for all chromosomes
}
}

## Generating metadatacolumns in GRanges file
{
gr <- rowRanges(vcfK) #generetes a GRanges object from the vcf-file
gr$GT <- geno(vcfK)$GT #Add a metadatacolumn for genotype (i.e. 0/0/1/1)
gr$HET <- is_het(gr$GT, na_is_false=TRUE) #Add a metadatacolumn for heterogeneity (TRUE = heterozygote, FALSE = homozygote)
gr$SNV <- isSNV(vcfK, singleAltOnly = FALSE) #Add a metadatacolumn for SNV (TRUE = SNV, FALSE = else)
gr$NEW <- grepl("chr", names(gr)) #Add a metadatacolumn for whether the mutation is new (TRUE) or already has an rs-number (FALSE)
gr$OVERLAP <- countOverlaps(gr) > 1
}

## Annotate the variants
{
gr_all <- locateVariants(gr, txdb, AllVariants()) #locate the variants
gr_ig <- locateVariants(gr, txdb,IntergenicVariants())
gr_c <- locateVariants(gr, txdb, CodingVariants())
table(gr_all$LOCATION) #To check where the variants are ##doesnt work; way to many variants -> 3.9Mc
barplot(table(gr_all$LOCATION)/1000000,ylab="Frequency (*10^6)", las = 2)

seqlevels(txdb) <- paste0("chr", c(1:22, "X"))
select(txdb, keys = keys, columns = columns(txdb), keytype = "GENEID")

}

#find variants 
length(gr[gr$HET == TRUE & gr$SNV == TRUE & gr$NEW == TRUE & gr$OVERLAP == FALSE])

## Plot various relationships within the variations
{
#Plot SNV vs other
b <- barplot(table(gr$SNV), 
             main = "SNVs",
             ylab = "Frequency",
             ylim = c(0,3500000),
             names.arg = c("Other (e.g. Indel)", "SNV"))
text(x=b, y=table(gr$SNV)+200000, labels = as.character(table(gr$SNV)))

#Plot nr. of Overlaps
b <- barplot(table(countOverlaps(gr)), 
             main = "Overlapping variants",
             xlab = "Number of overlaps +1",
             ylab = "Frequency",
             ylim = c(0,4200000),
                          )
text(x=b, y=table(countOverlaps(gr))+200000, labels = as.character(table(countOverlaps(gr))))

#Plot NEW vs rs
b <- barplot(table(gr$NEW), 
             main = "Novel variants (no rs-number)",
             ylab = "Frequency",
             ylim = c(0,4200000),
             names.arg = c("dbSNP138", "Novel"))
text(x=b, y=table(gr$NEW)+200000, labels = as.character(table(gr$NEW)))

#Plot Heterozygote vs Homozygote
b <- barplot(table(gr$HET), 
             main = "Variant heterozygosity",
             ylab = "Frequency",
             ylim = c(0,2500000),
             names.arg = c("Homozygote", "Heterozygote")
)
text(x=b, y=table(gr$HET)+200000, labels = as.character(table(gr$HET)))
}

## Plot Karyotype
{
pp <- getDefaultPlotParams(plot.type=1)
pp$topmargin <- 15
kp <- plotKaryotype(plot.type = 1, 
                    chromosomes = c("chr1","chr2","chr3"), 
                    main = "Mutation density plot (including SNVs and other variants)",
                    plot.params = pp
)
kp <- kpDataBackground(kp)
kp <- kpAxis(kp, ymin = 0, ymax = 1, )
kp <- kpAddLabels(kp, labels = "All", col = "gray", r0 = 0.75, label.margin = 0.01, srt = 45, cex = 2)
kp <- kpAddLabels(kp, labels = "Novel", col = "red", r0 = 0.35, label.margin = 0.01, srt = 45, cex = 2)
kpPlotDensity(kp, data=gr, col = alpha("gray", 1), window.size = 1000000, border = 1)
kpPlotDensity(kp, data=gr[gr$NEW == TRUE], col = alpha("red", 0.6), window.size = 1000000, border = 1)
# ---- # Karyoplot Intergenic vs Genes
pp <- getDefaultPlotParams(plot.type=1)
pp$topmargin <- 15
kp <- plotKaryotype(plot.type = 1, 
                    chromosomes = c("chr1","chr2","chr3"), 
                    main = "Mutation density plot (including SNVs and other variants)",
                    plot.params = pp
)
kp <- kpDataBackground(kp)
kp <- kpAxis(kp, ymin = 0, ymax = 1, )
kp <- kpAddLabels(kp, labels = "Intergenic", col = "gray", r0 = 0.75, label.margin = 0.01, srt = 45, cex = 2)
kp <- kpAddLabels(kp, labels = "Intragenic", col = "red", r0 = 0.35, label.margin = 0.01, srt = 45, cex = 2)
kpPlotDensity(kp, data=gr_ig, col = alpha("gray", 1), window.size = 2000000, border = 1)
kpPlotDensity(kp, data=gr_c, col = alpha("red", 0.6), window.size = 2000000, border = 1)


gr_tst <- gr[gr$HET == TRUE & gr$SNV == TRUE & gr$OVERLAP == FALSE 
             #& seqnames(gr) == "chr1"
             ]

pp <- getDefaultPlotParams(plot.type=4)
pp$topmargin <- 15
pp$ideogramheight <- 0
pp$data1inmargin <- 2
pp$data2inmargin <- 0
pp$leftmargin <- 0.2
kp <- plotKaryotype(plot.type = 4, 
                    chromosomes = paste0("chr", c(1:22, "X")), 
                    main = "Heterozygote/SNV/No_Overlap",
                    plot.params = pp,
                    labels.plotter = NULL,
                    ideogram.plotter = NULL
)
kp <- kpAddLabels(kp, labels = "QUAL Score", srt = 90)
kp <- kpDataBackground(kp)
kp <- kpAxis(kp, ymin = 0, ymax = max(gr_tst$QUAL)/10)
kp <- kpAddChromosomeNames(kp, srt = 90)
kp <- kpAddLabels(kp, labels = paste0("All (n = ", length(gr_tst), ")"), col = "1", r0 = 0.75, label.margin = 0.01, srt = 45, cex = 2)
kp <- kpAddLabels(kp, labels = paste0("New (n = ", length(gr_tst[gr_tst$NEW == T]),")"), col = "2", r0 = 0.25, label.margin = 0.01, srt = 45, cex = 2)
kp <- kpAddCytobandsAsLine(kp)
kp <- kpPoints(kp, 
               data = gr_tst,
               #chr = "chr1", 
               #x = start(gr_tst), 
               y = gr_tst$QUAL/max(gr_tst$QUAL)*10,
               col = alpha(colour = 1, alpha = 0.01)) #All
kp <- kpPoints(kp, 
               data = gr_tst[gr_tst$NEW == T],
               #chr = "chr1", 
               #x = start(gr_tst[gr_tst$NEW == T]), 
               y = gr_tst$QUAL[gr_tst$NEW == T]/max(gr_tst$QUAL)*10,
               col = alpha(colour = 2, alpha = 0.03)) #Only NEW
}

## Trimming of GRanges
ggr <- as(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene), "GRanges")
fgr <- subsetByOverlaps(gr, ggr, type = "within")
badGR <- GRanges(paste0("chr", c(1,2,3,4,5)), IRanges(c(249250600, 3,4,5,6), c(249250645, 6,7,8,9)))         
regdb <- read.table("~/projects/SuRE_K562/data/external/RegulomeDB.dbSNP132.Category3.txt")

#Summarizing Table
{
grp <- gr[gr$OVERLAP == F]
data <- matrix(c(
               length(grp[grp$NEW == F]),
               length(grp[grp$NEW == T]),
               length(grp[grp$NEW == T & grp$SNV == F]),
               length(grp[grp$NEW == T & grp$SNV == T]),
               length(grp[grp$NEW == T & grp$SNV == T & grp$HET == T]),
               length(grp[grp$NEW == T & grp$SNV == T & grp$HET == F]))
               , nrow = 2) #First fills the columns, then the rows
colnames(data) <- c("All", "Novel", "Novel SNV")
data_percentage <- apply(data, 2, function(x){x*100/sum(x)})

b <- barplot(data_percentage, 
             main = "Variants (non-overlapping)",
             xlab = NULL,
             ylab = "Percentage",
             ylim = c(0,130),
             names.arg = colnames(data),
)
text(x=b, 
     y=c(data_percentage[1,])+3, 
     labels = data[2,], srt = 0,
     #adj = c(0.5,2)
     )
text(x=b, 
     y=c(data_percentage[1,])-3, 
     labels = data[1,], srt = 0,
     #adj = c(0.5,1.5),
     col = 0
     )
text(x=b-0.35,
     y=10,
     labels = c("dbSNP138", "Other (e.g. Indel)", "Heterozygote"),
     srt = 60,
     adj = c(0,0.5),
     col = 0)
text(x=b-0.35,
     y=10,
     labels = c("dbSNP138", "Other (e.g. Indel)", "Heterozygote"),
     srt = 60,
     adj = c(0,0.5),
     col = 0)
text(x=b,
     y=105,
     labels = c("Novel","SNV", "Homozygote"),
     srt = 60,
     adj = c(0,0.5))
}
#apply(as.matrix(b), 1, function(x){lines(c(x,x),c(100,105))})






