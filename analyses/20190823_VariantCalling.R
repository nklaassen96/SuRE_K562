#I use the data.frame file containing 69,546 mutations that are not annotated in the dbSNP138(?) database.

#Install Bioconductor
if (!requireNamespace("BiocManager", quietly = FALSE))
  install.packages("BiocManager")
BiocManager::install()

#install TxDB database
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene"))

#Load the required UCSC transcript database
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#install VariantTools
BiocManager::install(c("VariantTools"))
library(VariantTools)

#install VariantAnnotation?
library(VariantAnnotation)

#read the vcf file that was downloaded
vcf <- readVcf("~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF606RIC.vcf.gz", 
               genome = "hg19")
locCV <- locateVariants(vcf, txdb, CodingVariants())
locIG <- locateVariants(vcf, txdb, IntergenicVariants())
