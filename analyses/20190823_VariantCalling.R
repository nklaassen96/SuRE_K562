#I use the data.frame file containing 69,546 mutations that are not annotated in the dbSNP138(?) database.

##Install Bioconductor
#if (!requireNamespace("BiocManager", quietly = FALSE))
#  install.packages("BiocManager")
#BiocManager::install()

##install TxDB database
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene",INSTALL_opts = c('--no-lock'))

##install VariantTools
#BiocManager::install(c("VariantTools"))

#Load required libraries
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantTools)
#library(VariantAnnotation) #is within VariantTools?

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#read the vcf file that was downloaded
vcfK <- readVcf("~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF606RIC.vcf.gz", 
               genome = "hg19")
#locCV <- locateVariants(vcf, txdb, CodingVariants())
#locIG <- locateVariants(vcf, txdb, IntergenicVariants())

## Testing Set ##
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
geno(vcf) #shows al genotypes as described in FORMAT fields
sapply(geno(vcf), class)
info(vcf) #shows some sort of summary of all the data. Rows are indicated by rs-numbers or chromosome positions
rowRanges(vcf) #Gives the Granges of all rows (what is the difference with head(Rowranges)???)
as.data.frame(rowRanges(vcf)) #generates dataframe of these rows
vcf[1:100,] #subset few rows within a VCF
seqlevels(vcf) <- "chr22" #set the seqlevels all to "chr22" instead of 22 (this is the required format)
rd <- rowRanges(vcf)
tst <- locateVariants(vcfK, txdb, AllVariants()) #locate the variants
table(tst$LOCATION) #To check where the variants are ##doesnt work; way to many variants -> 3.9M
barplot(table(tst$LOCATION)/1000000,ylab="Frequency (*10^6)", las = 2)

#find #of overlapping variations
ts <- table(countOverlaps(vcfK))
barplot(ts[2:8], xlab = "# of Overlaps", ylab = "Frequency")
ts <- as.data.frame(table(countOverlaps(vcfK)))

#try to remove overlapping variations
rr <- rowRanges(vcfK)
hits <- findOverlaps(rr, drop.self = TRUE, drop.redundant = TRUE)
ovpairs <- Pairs(rr, rr, hits = hits)
