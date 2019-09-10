
#### STARTUP ####

  # This is to load the required libraries, generate a 
  # tabix file that is linked with the VCF

    library(VariantTools)
    library(vcfR)
    library(R.utils)
    fl <- "~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF606RIC.vcf.gz"
    vcf <- readVcf(fl, genome = "hg19")
    compressVcf <- bgzip(fl, tempfile())
    idx <- indexTabix(compressVcf, "vcf")
    tab <- TabixFile(compressVcf, idx)

#### WRITE .VCF ALL ####
  
  # This function is used to generate vcf-files per chromosome
  # as input for the SuRE pipeline. The gzip-function can be 
  # used to generate gzip-files from the vcf-files
  
  for (chrom in paste0("chr", c(1:22, "X"))){
    ranges <- as(seqinfo(vcf)[chrom], "GRanges")
    #param <- ScanVcfParam(which = ranges)
    param <- ScanVcfParam(which = ranges, geno = "GT")
    vcf.chr <- readVcf(tab, genome = "hg19", param = param)
    dir = paste0("~/projects/SuRE_K562/data/interim/SNPs/vcf/All_K562_NonPhased_", chrom, ".vcf")
    writeVcf(vcf.chr, dir)
    print(paste0("All_K562_NonPhased_", chrom, ".vcf", " was successfully generated"))
    gzip(dir, remove = TRUE)
    print(paste0("All_K562_NonPhased_", chrom, ".vcf", " was successfully gzipped to", " All_K562_NonPhased_", chrom, ".vcf.gz"))
    
  }
  
  # Generating .vcf.gz files from the .vcf files  
      
  for (chrom in paste0("chr", c(1:22, "X"))){
    dir = paste0("~/projects/SuRE_K562/data/interim/SNPs/vcf/All_K562_NonPhased_", chrom, ".vcf")
    gzip(dir, remove = FALSE)
    print(paste0("All_K562_NonPhased_", chrom, ".vcf", " was successfully gzipped to", " All_K562_NonPhased_", chrom, ".vcf.gz"))
  }
  
#### WRITE .TXT ALL ####
  
  # This function is used to generate tab-delimited txt 
  # files per chromosome for all mutations. These should
  # be in the format start | ref | alt
  
    gr <- granges(vcf)
    gr_no.overlap <- gr[countOverlaps(gr) == 1]
    
    # start the loop here
    
    for (chrom in paste0("chr", c(1:22, "X"))){
      gr_chrom <- gr_no.overlap[seqnames(gr_no.overlap) == chrom]
      start <- start(ranges(gr_chrom))
      ref <- as.character(gr_chrom$REF)
      altlist <- CharacterList(gr_chrom$ALT)
      alt <- unstrsplit(altlist, sep=",")
      df <- data.frame(start, ref, alt, stringsAsFactors = FALSE)
      dir = paste0("~/projects/SuRE_K562/data/interim/SNPs/vcf/SNPs/All_NoOverlap_K562_NonPhased_", chrom, ".txt")
      write.table(df, file = dir, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      print(paste0("All_NoOverlap_K562_NonPhased_",chrom,".txt", " was successfully generated"))
    }  
#### WRITE .TXT NEW ####
    
  # This function is used to generate tab-delimited txt
  # files per chromosome for the novel mutations only
    
    gr_no.overlap_new <- gr_no.overlap[!grepl("rs", names(gr_no.overlap))]
    
    for (chrom in paste0("chr", c(1:22, "X"))){
      gr_chrom <- gr_no.overlap_new[seqnames(gr_no.overlap_new) == chrom]
      start <- start(ranges(gr_chrom))
      ref <- as.character(gr_chrom$REF)
      altlist <- CharacterList(gr_chrom$ALT)
      alt <- unstrsplit(altlist, sep=",")
      df <- data.frame(start, ref, alt, stringsAsFactors = FALSE)
      dir = paste0("~/projects/SuRE_K562/data/interim/SNPs/vcf/SNPs/New_NoOverlap_K562_NonPhased_", chrom, ".txt")
      write.table(df, file = dir, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      print(paste0("New_NoOverlap_K562_NonPhased_",chrom,".txt", " was successfully generated"))
    }

#### TEST VCF WITHOUT EXTRA METADATA ####
    
    param <- ScanVcfParam(geno = "GT")
    vcf.tst <- readVcf(tab, genome = "hg19", param = param)
    dir = paste0("~/projects/SuRE_K562/data/interim/SNPs/vcf/All_", "TEST", ".vcf")
    writeVcf(vcf.tst, dir)
    