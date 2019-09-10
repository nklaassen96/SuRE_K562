library(vcfR)
library(VariantTools)
vcfK <- readVcf("~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF606RIC.vcf.gz", genome = "hg19")
gr <- rowRanges(vcfK)

#### LOOP FOR GENERATION OF TFBS ####

files <- list.files("~/projects/SuRE_K562/data/external/bed_files/bed")
out <- NULL
for (file in files){
  tab <- read.table(paste0("~/projects/SuRE_K562/data/external/bed_files/bed/", file, sep=""), stringsAsFactors = FALSE)
  tab$File <- rep(file, nrow(tab))
  if (grepl("hg19", tab$V4[1]) == FALSE){
    next
  }
  out <- rbind(out,tab)
}

SEQ   <- out$V1
START <- out$V2+1 #0-based to 1-based
END   <- out$V3

gr_tfbs <- GRanges(seqnames = SEQ, 
                   ranges = IRanges(start = START, end = END), 
                   File = out$File, 
                   strand = out$V6)  
gr_tfbs <- sortSeqlevels(gr_tfbs)
gr_tfbs <- sort(gr_tfbs)
table(countOverlaps(gr, gr_tfbs))
#OVERLAP <- subsetByOverlaps(gr, gr_tfbs)  
OVERLAP <- gr[countOverlaps(gr, gr_tfbs) > 0]

#### TRY TO ANNOTATE THE OVERLAPPING SNVS WITH THE CORRESPONDING TFBS ####

gr_new <- gr[gr$NEW == T]
f <- findOverlaps(gr_new, gr_tfbs)
q <- queryHits(f)
s <- subjectHits(f)

gr_tfbs$File[s] 

OVERLAP.REDUNDANT <- gr_new[q]
OVERLAP.REDUNDANT$FILE <- gr_tfbs$File[s]




