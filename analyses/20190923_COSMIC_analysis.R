## script to analyse the COSMIC database
# Data was downloaded from https://cancer.sanger.ac.uk/cosmic/download  and then `Non coding variants`

library(data.table)
library(GenomicRanges)

cosmic.noncoding <- fread(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/COSMIC/CosmicNCV.tsv.gz")

# Remove columns that are unnessesary
cosmic.noncoding.filter.rows <- cosmic.noncoding
cosmic.noncoding.filter.rows[,c(1,3:11,13,21:28)] <- NULL
head(cosmic.noncoding.filter.rows)

# Remove Insertions and deletions
filter.SNV <- nchar(cosmic.noncoding.filter.rows$WT_SEQ) == 1 & nchar(cosmic.noncoding.filter.rows$MUT_SEQ) == 1
cosmic.noncoding.filter.rows.SNV <- cosmic.noncoding.filter.rows[filter.SNV]
head(cosmic.noncoding.filter.rows.SNV)

# Generate genomic positions
genomic.locations <- paste0("chr", cosmic.noncoding.filter.rows.SNV$`genome position`)
cosmic.granges <- as(genomic.locations, "GRanges")

# Add metadata to the GRanges file
cosmic.granges$


cosmic.granges.chrom <- cosmic.granges[seqnames(cosmic.granges) %in% paste0("chr", c(1:22)),]
cosmic.granges.chrom <- sort(sortSeqlevels(cosmic.granges.chrom))
saveRDS(cosmic.granges.chrom, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/R_Objects/cosmic.granges.RDS")

### import data from K562 Paper

# Generate GRanges file with known (rs-number) SNV (one base mutations)
library(vcfR)
vcfK <- readVcf("~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF606RIC.vcf.gz", 
                genome = "hg19")
gr <- rowRanges(vcfK)
gr$SNV <- isSNV(vcfK, singleAltOnly = FALSE)
gr$NEW <- !grepl("rs", names(gr))
gr.rs.SNV <- gr[gr$SNV == T & gr$NEW == F]

# Find overlaps between K562 mutations and COSMIC database
overlap <- countOverlaps(gr.rs.SNV, cosmic.granges.chrom)
tab <- table(overlap)
sum(tab[2:length(tab)])

cosmic.overlap.snp.id <- names(overlap[overlap > 0])


