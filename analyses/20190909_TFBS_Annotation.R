library(vcfR)
library(VariantTools)
vcfK <- readVcf("~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF606RIC.vcf.gz", genome = "hg19")
gr <- rowRanges(vcfK)

#### LOOP FOR GENERATION OF TFBS ####

files <- list.files("~/projects/SuRE_K562/data/external/TFBS_Annotation/bed_files/bed")
out <- NULL
for (file in files){
  tab <- read.table(paste0("~/projects/SuRE_K562/data/external/TFBS_Annotation/bed_files/bed/", file, sep=""), stringsAsFactors = FALSE)
  tab$File <- rep(file, nrow(tab))
  if (grepl("hg19", tab$V4[1]) == FALSE){
    next
  }
  out <- rbind(out,tab)
}

#### GENERATE OVERLAP FILE ####

SEQ   <- out$V1
START <- out$V2+1 #0-based to 1-based
END   <- out$V3

gr_tfbs <- GRanges(seqnames = SEQ, 
                   ranges = IRanges(start = START, end = END), 
                   File = out$File, 
                   strand = out$V6)  
gr_tfbs <- sortSeqlevels(gr_tfbs)
gr_tfbs <- sort(gr_tfbs)
OVERLAP <- gr[countOverlaps(gr, gr_tfbs) > 0]

#### TRY TO ANNOTATE THE OVERLAPPING SNVS WITH THE CORRESPONDING TFBS ####

gr_new <- gr[gr$NEW == F & gr$SNV == T]
f <- findOverlaps(gr_new, gr_tfbs)
q <- queryHits(f)
s <- subjectHits(f)
length(unique(q))/length(gr_new)*100

gr_tfbs$File[s] 

OVERLAP.REDUNDANT <- gr_new[q]
OVERLAP.REDUNDANT$FILE <- gr_tfbs$File[s]

#### Try to put in the sequence ####

TF <- "MA0528.1.bed"
overlaptst <- OVERLAP.REDUNDANT[OVERLAP.REDUNDANT$FILE == TF]
overlaptst <- OVERLAP.REDUNDANT
f <- findOverlaps(overlaptst, gr_tfbs)
q <- queryHits(f)
s <- subjectHits(f)
overlaptstredundant <- overlaptst[q]
overlaptstredundant$POS <- start(overlaptst[q])-start(gr_tfbs[s])+1


#### figure out whether they are from the human genome or not

files <- list.files("~/projects/SuRE_K562/data/external/TFBS_Annotation/bed_files/bed")
tflist <- NULL
for (file in files){
  tab <- read.table(paste0("~/projects/SuRE_K562/data/external/TFBS_Annotation/bed_files/bed/", file, sep=""), stringsAsFactors = FALSE, nrows = 1)
  tab$File <- rep(file, nrow(tab))
  if (grepl("hg19", tab$V4[1]) == FALSE){
    next
  }
  tflist <- rbind(tflist,tab)
}

# retrieve the total spanning of TFBS on the human genome
ch1 <- gr_tfbs[seqnames(gr_tfbs) == "chr2"]
ch1r <- reduce(ch1)
sum(end(ch1r)-start(ch1r))
gr_tfbs_r <- reduce(gr_tfbs)
reg <- tapply(gr_tfbs_r, seqnames(gr_tfbs_r), function(x){sum(end(x)-start(x))} )

#### PWMs ####
PWM(pwm)
con <- consensusMatrix(HNF4alpha)
ppm <- apply(con, 2, function(x){x/sum(x)})
pwm <- PWM(HNF4alpha)
seqLogo::seqLogo(ppm)
