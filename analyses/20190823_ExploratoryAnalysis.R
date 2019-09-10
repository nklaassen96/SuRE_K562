########################################################
#####   STARTUP   ######################################
########################################################

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
vcfK <- readVcf("~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF241LKI.vcf.gz", 
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
  #tapply(gr$GT, seqnames(gr), summary) # Show genotypes per chromosome
  #find variants 
  #length(gr[gr$HET == TRUE & gr$SNV == TRUE & gr$NEW == TRUE & gr$OVERLAP == FALSE])
}
}

########################################################
#####   Generating metadatacolumns in GRanges file   ###
########################################################

{
gr <- rowRanges(vcfK) #generetes a GRanges object from the vcf-file
gr$GT <- geno(vcfK)$GT #Add a metadatacolumn for genotype (i.e. 0/0/1/1)
gr$HET <- is_het(gr$GT, na_is_false=TRUE) #Add a metadatacolumn for heterogeneity (TRUE = heterozygote, FALSE = homozygote)
gr$SNV <- isSNV(vcfK, singleAltOnly = FALSE) #Add a metadatacolumn for SNV (TRUE = SNV, FALSE = else)
gr$NEW <- !grepl("rs", names(gr)) #Add a metadatacolumn for whether the mutation is new (TRUE) or already has an rs-number (FALSE)
gr$OVERLAP <- countOverlaps(gr) > 1
}

########################################################
#####   Annotate the variants ##########################
########################################################

{
for (i in 1:2) {
  
  if (i == 1) {variant = gr} else {variant = gr[gr$NEW == T]}
#variant <- gr
#gr_all <- locateVariants(variant, txdb, AllVariants()) #locate the variants

gr_in <- locateVariants(variant, txdb,IntronVariants())  
gr_inu <- gr_in[!duplicated(gr_in$QUERYID)]
  
gr_ig <- locateVariants(variant, txdb,IntergenicVariants())
gr_igu <- gr_ig[!duplicated(gr_ig$QUERYID)] ## Unique_Intergenic (no values discarded)

gr_c <- locateVariants(variant, txdb, CodingVariants()) 
gr_cu <- gr_c[!duplicated(gr_c$QUERYID)] ## Unique_Coding (28589 values discarded)

gr_p <- locateVariants(variant, txdb, PromoterVariants())
gr_pu <- gr_p[!duplicated(gr_p$QUERYID)]

gr_3 <- locateVariants(variant, txdb, ThreeUTRVariants())
gr_3u <- gr_3[!duplicated(gr_3$QUERYID)]

gr_5 <- locateVariants(variant, txdb, FiveUTRVariants())
gr_5u <- gr_5[!duplicated(gr_5$QUERYID)]

gr_sp <- locateVariants(variant, txdb, SpliceSiteVariants())
gr_spu <- gr_sp[!duplicated(gr_sp$QUERYID)]

count <- c(length(variant),
           length(gr_inu),
           length(gr_igu), 
           length(gr_pu), 
           length(gr_cu),
           length(gr_3u),
           length(gr_5u),
           length(gr_spu))
percentage <-  c(count/count[1]*100)

var <- matrix(data=c(count, percentage), 
              ncol = 2, 
              dimnames = list(c("All", "Intron", "Intergenic", "Promoter", "Coding", "3-UTR", "5-UTR", "SpliceSite"),
                              c("Count", "Percentage")))
var
varsum <- matrix(nrow = 8, ncol = 2)
varsum[,i] <- var[,2]
if (i == 1){add = F} else {add = T}
b <- barplot(var[,2], las = 2, ylab = "Percentage", ylim = c(0,130), main = "All variants", add = add , col = i+1)
text(b, y = var[,2]+3 ,labels = round(var[,2], digits = 1), adj = c(0.5,0.5))
}
keys <- c("100033416", "100033417", "100033420")
chrom <- paste0("chr", c(1:22, "X"))

total <- sum(seqlengths(gr)) #this includes still chromosome M and Y
tst <- exons(txdb, filter=list(tx_chrom = chrom))
tstr <- reduce(tst) #Merges overlapping regions
exon <- sum(width(tstr))
select(txdb, keys = keys, columns = columns(txdb), keytype = "GENEID")
#seqlevels(txdb) <- paste0("chr", c(1:22, "X"))
#select(txdb, keys = keys, columns = columns(txdb), keytype = "GENEID")

}

###############################################################
#####   Plot various relationships within the variations   ####
###############################################################

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

###############################################################
#####   Plot Karyotype   ######################################
###############################################################

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

########################################################
#####   Summarizing Table   ############################
########################################################

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
#apply(as.matrix(b), 1, function(x){lines(c(x,x),c(100,105))})
}

########################################################
#####   QUAL SCORES  ###################################
########################################################

h1 = gr$QUAL
hist(h1, xlim = c(0, 6000), col = 3, breaks = 2000, xlab = "QUAL Score", main = "Variants")
h2 = gr[gr$NEW== T]$QUAL
hist(h2, xlim = c(0, 6000), col = 2, breaks = 2000, add = T)
legend("topright", c("All", "New"), fill = c("green","red"))

########################################################
#####   HETEROZYGOTIC MUTATION SPREAD   ################
########################################################

{
length(gr[seqnames(gr)=="chr3" & gr$HET == T])
allho <- tapply(gr$SNV[gr$HET == F], seqnames(gr[gr$HET == F]), sum) #retrieve nr. of SNVs for all chromosomes
newho <- tapply(gr$SNV[gr$HET == F & gr$NEW == T], seqnames(gr[gr$HET == F & gr$NEW == T]), sum)
allhet <- tapply(gr$SNV[gr$HET == T], seqnames(gr[gr$HET == T]), sum) #retrieve nr. of SNVs for all chromosomes
newhet <- tapply(gr$SNV[gr$HET == T & gr$NEW == T], seqnames(gr[gr$HET == T & gr$NEW == T]), sum)

barplot(allhet, las = 2, ylab = "Frequency", main = "Heterozygotic mutaions")
barplot(newhet, add = T, col = 2, las = 2)
legend("topright", legend = c("All", "New"), fill = c("gray", "red"))
box()
barplot(newhet/allhet*100, las= 2, ylab = "New mutations (%)", main = "Heterozygote SNVs")
barplot(newho/allho*100, las= 2, ylab = "New mutations (%)", main = "Homozygote SNVs")
}

########################################################
#####   FILTERING VCF   ################################
########################################################



