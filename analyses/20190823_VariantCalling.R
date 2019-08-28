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

## Testing Set START ##
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
  lengths(gr$ALT) #determine the amount of ALT sequences
  table(width(gr_snv)) #contingency table for the length of the DNA fragments (REF or ALT??), ik vermoed REF
  expand(vcf) #generates an "expanded-VCF" object that shows a row for every alternative sequence
  
## Testing Set END ##

  
  
gr <- rowRanges(vcfK) #generetes a GRanges object from the vcf-file

## Differentiate between SNVs and other variants
gr_snv <- gr[isSNV(vcfK, singleAltOnly = FALSE)] #isSNV(vcfK) returns TRUE or FALSE for whether it is a SNV; only works for VCF-objects, not granges, so should be performed first?
table(lengths(gr_snv$ALT)) #shows the number of alternative sequences
table(isSNV(vcfK, singleAltOnly = FALSE)) #shows amount of SNVs compared to other variations

## Remove overlapping regions
gr_no_overlap <- gr[countOverlaps(gr) == 1] #filter the GRanges object for where the number of counts equals 1 or less
table(countOverlaps(gr))
table(countOverlaps(gr_no_overlap))  
barplot(table(countOverlaps(gr))[2:8], xlab = "# of Overlaps", ylab = "Frequency")

## Find reference-SNPs numbers or non-reference "new" variations
gr_rs <- gr[grep("rs", names(gr))]
gr_new <- gr[grep("chr", names(gr))]
length(gr_rs)+length(gr_new) == length(gr) 
  
## Annotate the variants
gr_all <- locateVariants(gr, txdb, AllVariants()) #locate the variants
gr_ig <- locateVariants(gr, txdb,IntergenicVariants())
table(gr_all$LOCATION) #To check where the variants are ##doesnt work; way to many variants -> 3.9Mc
barplot(table(gr_all$LOCATION)/1000000,ylab="Frequency (*10^6)", las = 2)

#Filter for heterozygotes


### SOME EXPLORATORY ANALYSIS

matt <- matrix(c(length(gr_snv),2,3,3088312,3020306,69553), 
       dimnames = list(c("SNV:All", "SNV:dbSNP138", "SNV:Novel"),c("This Analysis", "K562-paper")),
       nrow = 3)
diff <- mat[,1]-mat[,2]
### END



#https://www.bioconductor.org/help/course-materials/2014/BioC2014/Lawrence_Tutorial.pdf
#shows what you can do with variants and vcf-files




