## script to analyse the COSMIC database
# Data was downloaded from https://cancer.sanger.ac.uk/cosmic/download  and then `Non coding variants` or `COSMIC Mutation Data` from the COSMIC v90 database

library(data.table)
library(GenomicRanges)
dir = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/COSMIC/"


cosmic.noncoding <- fread(file = paste0(dir, "CosmicNCV.tsv.gz"))
cosmic.coding <- fread(file = paste0(dir, "CosmicMutantExport.tsv.gz"))
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


'
***************************************
  CURRENT ANNOUNCEMENTS
***************************************
  Revised: Sept 17, 2018
****************************************
  
  This document describes the FTP repository of dbSNP data in the
sections listed below:
  
  I.   FTP SITE DESCRIPTION AND FTP UPDATE FREQUENCY
II.   DIRECTORY STRUCTURE
IV.   FTP SITE REVISION HISTORY

********************************************************************************
  
  I. FTP SITE DESCRIPTION AND FTP UPDATE FREQUENCY

Access to the NCBI FTP site is available via the web or anonymous FTP. 
The URL/host addresses are:
  
  World Wide Web:     ftp://ftp.ncbi.nlm.nih.gov/snp/
  Anonymous FTP:      host ftp.ncbi.nlm.nih.gov
cd snp

Announcements of the release of new builds and notification of corrections
to existing data content will be posted on a public mail list. To receive
these notifications, you can subscribe to the dbSNP announcement list at
https://www.ncbi.nlm.nih.gov/mailman/listinfo/dbsnp-announce.

********************************************************************************
  
  II. DIRECTORY STRUCTURE
DIRECTORIES: 
  
  /bin             Contains demo software tools for using RefSNP JSON files.

/specs           Contains the RefSnp API schema used by RefSNP JSON files.
/refsnp_specification.yaml

/latest_release  Contains the most recent release of human SNP data, in VCF and
API JSON format, along with the release notes:
  /JSON
/VCF
/release_notes.txt


SUBDIRECTORIES:
  
  /JSON              RefSNP JSON files. Refer to JSON_README.txt for details.

/VCF               RefSNP VCF files for GRC (Genome Reference Consortium) human assembly
37 (GCF_000001405.25) and 38 (GCF_000001405.38). Files are compressed
by bgzip and with the tabix index.

IV. FTP SITE REVISION HISTORY:
  Rev: Sept 17, 2018
Rewrite the readme file for redesigned dbSNP.
********************************************************************************
  
  dbSNP build 153 release notes

Organism name: Homo sapiens
Taxonomy ID: 9606

1. RefSnp (RS)

Total RS: 695448618

1.1 RS counts by location:
chr1	52293824
chr2	55992290
chr3	45826475
chr4	44074378
chr5	41348564
chr6	38673067
chr7	37008263
chr8	35112021
chr9	29024921
chr10	30877811
chr11	31641810
chr12	30588046
chr13	22584505
chr14	20572956
chr15	19243104
chr16	21175632
chr17	18769467
chr18	17863560
chr19	14340102
chr20	14679288
chr21	8797523
chr22	9140034
chrX	25851875
chrY	1297482
chrM	3319
Alt Only	31583
Patch	7959
Not On	39896
PAR	813952
Unplaced	93409
Unlocalized	186701

NOTE: Assembly term (ALT, PAR, etc.) definitions (https://www.ncbi.nlm.nih.gov/grc/help/definitions/).

1.2 RS counts by type:
  Live RS	667501404
Unsupported RS	117770
Withdrawn RS	6854924
Locationless RS	1286
Merged RS	20973234

2. SubSnp (SS)

Total SS: 1997618652
Unmapped SS: 60209

RS Type Definitions:
  Live = RS has location on reference sequence 
Merged = RS merged to existing RS due to improved clustering algorithm or possibly a change to the reference sequence that would result in identical canonical alleles (e.g. updated repeat regions).
Unsupported = No Submitted SNP (SS) matched any of the RS alleles. (Same causes as merging.)
Locationless = An older RS where the location couldnt be obtained from SS and was not available for the build.
Withdrawn = All SS that belong to the RS cluster were withdrawn.

All above RS including non-Live records have history for traceability.
'
