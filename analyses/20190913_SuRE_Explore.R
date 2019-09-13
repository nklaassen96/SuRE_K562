
dir <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/SuRE_output"
files <- list.files(dir)

sure = NULL
for (chrom in paste0("chr", c(1:13,15:18, 21:22, "X"))){
  table <- read.table(paste0(dir, "/", "SuRE-counts_", chrom, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sure <- rbind(sure, table)
  print(paste0(chrom, " is finished"))
}


gr3 <- gr[gr$SNV == T & seqnames(gr) == "chr3"]
table(grepl(",", table$SNPbase)) #How many rows have multiple alternative alleles
                       # 1505 have multiple alleles
nrow(table)            # 3131457 Total amount of rows
table(table$SNPrelpos) # 3097099 entries have no SNPrelpos
table(table$SNPbase)   # 3097099 entries have no SNPbase etc...
length(unique(table$SNP_ID)) # 3935 unique SNP_IDs 
length(unique(names(gr3)))   # 171206 unique variants  --> 2% is covered


tab1a <- table[!grepl(",", table$SNPbase),]
table(tab1a$SNPbase)   # Many bases are not the normal ACTG due to low sequence quality
table(tab1a$SNPvar)

### Total (minus chr14 and chr19)
length(unique(sure$SNP_ID))
length(unique(names(gr)))
length(unique(sure$SNP_ID))/length(unique(names(gr)))*100 # % coverage of IDs (not sure if i should take indels into account)
sureuni <- sure[!grepl(",", table$SNPbase),]