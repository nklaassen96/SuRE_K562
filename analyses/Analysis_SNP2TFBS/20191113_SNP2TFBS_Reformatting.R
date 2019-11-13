# this script is used to compare all variants and the raqtl specific ones 
# in the SNP2TFBS library

## Data preparation

# First load the datasets of (1) all variants, (2) K562 raQTLs and (3) HepG2 raQTLs
all.variants <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/All.10-1000.elements.20191113.RDS")
raqtl.k562 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/K562.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")
raqtl.hepg2 <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indel_raQTLs/HepG2.raqtl.10-1000.elements.4.min.SuREexpr.20191113.RDS")

# These datasets still contain snps and indels, but I am interested in the indels only
all.variants <- all.variants[all.variants$snp.type == "indel",]
raqtl.k562 <- raqtl.k562[raqtl.k562$snp.type == "indel",]
raqtl.hepg2 <- raqtl.hepg2[raqtl.hepg2$snp.type == "indel",]





## Use SNP2TFBS to discover motif disruption


#  

# 1. SNP2TFBS requires a .txt file with the SNP_IDs. Therefore I Write table to txt
write.table(x = all.variants$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SNP2TFBS/indel.id.all.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = raqtl.k562$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SNP2TFBS/indel.id.raqtl.562.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = raqtl.hepg2$SNP_ID, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SNP2TFBS/indel.id.raqtl.hepg2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# These written tables can be transfered to the Windows PC with WinSCP

# 2. do this as input on https://ccg.epfl.ch/snp2tfbs/snpselect.php

# 3. download the snp2tfbs files

download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_14161.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/sec.motif.alteration.indel.raqtl.k562.txt")
download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_25813.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/sec.motif.alteration.indel.raqtl.hepg2.txt")
download.file(url = "https://ccg.epfl.ch/snp2tfbs/wwwtmp/match_output_27429.txt", destfile = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/sec.motif.alteration.indel.all.txt")

# 4. save the downloaded files as Robjects
tfbs.raqtl.k562 <-  fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.raqtl.k562.txt", select = c(7,1,2,4,5,6))
tfbs.raqtl.hepg2 <- fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.raqtl.hepg2.txt", select = c(7,1,2,4,5,6))
tfbs.all <-         fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SNP2TFBS/motif.alteration.indel.all.txt", select = c(7,1,2,4,5,6))

# 5. perform fisher exact test for enrichtment of motif alteration in raQTLs

for (tfbs.raqtl in c("tfbs.raqtl.k562","tfbs.raqtl.hepg2")){
  
  raqtl.motifaltering <- length(unique(get(tfbs.raqtl)$snp.id))
  nonraqtl.motifaltering <- length(unique(tfbs.all$snp.id)) - length(unique(get(tfbs.raqtl)$snp.id))
  
  if (tfbs.raqtl == "tfbs.raqtl.k562"){
    raqtl.nonmotifaltering <- length(unique(raqtl.k562$SNP_ID)) - length(unique(get(tfbs.raqtl)$snp.id))} else {
      raqtl.nonmotifaltering <- length(unique(raqtl.hepg2$SNP_ID)) - length(unique(get(tfbs.raqtl)$snp.id))}
  
  nonmotifaltering <- length(unique(all.variants$SNP_ID)) - length(unique(tfbs.all$snp.id))
  
  nonraqtl.nonmotifaltering <- nonmotifaltering - raqtl.nonmotifaltering
  
  mat <- matrix(data = c(raqtl.motifaltering,
                         nonraqtl.motifaltering,
                         raqtl.nonmotifaltering,
                         nonraqtl.nonmotifaltering), ncol = 2, nrow = 2)
  colnames(mat) <- c("motif.altering", "non-motif.altering")
  rownames(mat) <- c("raqtl", "non-raqtl")
  
  if (tfbs.raqtl == "tfbs.raqtl.k562"){enrichment.k562.pvalue <- fisher.test(mat)$p.value
  } else {enrichment.hepg2.pvalue <- fisher.test(mat)$p.value}
}



# 6. Reformat the dataframe for a row per transcription factor for the above 3 files

for (tfbs.file in c("tfbs.raqtl.k562", "tfbs.raqtl.hepg2", "tfbs.all")){
  
  df <- get(tfbs.file)
  
  x <- unlist(strsplit(df$V6, split = ";"))
  match.idx <- seq(1,length(x), by= 3)
  tf.idx <- seq(2,length(x), by=3)
  diff.score.idx <- seq(3,length(x), by = 3)
  
  col.match <- as.numeric(str_remove_all(string = x[match.idx], pattern = "MATCH="))
  col.tf <- str_remove_all(x[tf.idx], pattern = "TF=")
  col.diff.score <- str_remove_all(x[diff.score.idx], pattern = "ScoreDiff=")
  
  df$V6 <- NULL
  df$m <- col.match
  df$t <- col.tf
  df$d <- col.diff.score
  colnames(df) <- c("snp.id", "chrom","snp2tfbs.pos", "ref", "alt", "tfbs.match", "transcription.factor", "score.differences")
  df <- separate_rows(df, transcription.factor, score.differences, sep = ",")
  class(df$score.differences) <- "numeric"
  
  assign(x = tfbs.file, value = df)  
}


# 7. Find the maximum scoredifferences and plug them into the raQTL / total INDEL dataframes

### For all variants

all.max.scorediff




### For K562 raQTL (1st is max value, 2nd is mean value)
k562.max.scorediff <- tapply(tfbs.raqtl.k562$score.differences, tfbs.raqtl.k562$snp.id, function(x){x[which.max(abs(x))]})
k562.max.scorediff <- tapply(tfbs.raqtl.k562$score.differences, tfbs.raqtl.k562$snp.id, mean)

for (i in c(1:nrow(raqtl.k562))){
  
  print(i)
  # find the snp.id and corresponding scorediff from the tfbs dataframe
  snp.idx <- raqtl.k562[i,"SNP_ID"]
  scorediff <- k562.max.scorediff[snp.idx]
  
  # if the value is present in the tfbs dataframe, plug the values in
  if (!is.na(scorediff)){
    
    raqtl.k562[i,"k562.max.scorediff"] <- scorediff
    raqtl.k562[i,"tf.match"] <- tfbs.raqtl.k562[tfbs.raqtl.k562$snp.id == snp.idx, tfbs.match][1]
    
    # if the maximum value of the sure expression is in the ref allele AND the maximum
    # difference in tfbs score is negative, this means that the sure expression is in 
    # concordance with the tfbs scores
    
    if (raqtl.k562[i, "k562.max"] == "ref" & sign(raqtl.k562[i, "k562.max.scorediff"]) == -1) { raqtl.k562[i,"motif.concordance"] <- "concordance"} 
    if (raqtl.k562[i, "k562.max"] == "ref" & sign(raqtl.k562[i, "k562.max.scorediff"]) == 1)  { raqtl.k562[i,"motif.concordance"] <- "non.concordance"} 
    if (raqtl.k562[i, "k562.max"] == "alt" & sign(raqtl.k562[i, "k562.max.scorediff"]) == -1) { raqtl.k562[i,"motif.concordance"] <- "non.concordance"} 
    if (raqtl.k562[i, "k562.max"] == "alt" & sign(raqtl.k562[i, "k562.max.scorediff"]) == 1)  { raqtl.k562[i,"motif.concordance"] <- "concordance"}
    if (sign(raqtl.k562[i, "k562.max.scorediff"]) == 0 )                                      { raqtl.k562[i,"motif.concordance"] <- "undetermined"} 
    
    
    
  } else {
    
    raqtl.k562[i, "motif.concordance"] <- "no motifs"
  }
  
  
  
  
}

### For HepG2 raQTL
hepg2.max.scorediff <- tapply(tfbs.raqtl.hepg2$score.differences, tfbs.raqtl.hepg2$snp.id, function(x){x[which.max(abs(x))]})
hepg2.max.scorediff <- tapply(tfbs.raqtl.hepg2$score.differences, tfbs.raqtl.hepg2$snp.id, mean)

for (i in c(1:nrow(raqtl.hepg2))){
  
  print(i)
  # find the snp.id and corresponding scorediff from the tfbs dataframe
  snp.idx <- raqtl.hepg2[i,"SNP_ID"]
  scorediff <- hepg2.max.scorediff[snp.idx]
  
  # if the value is present in the tfbs dataframe, plug the values in
  if (!is.na(scorediff)){
    
    raqtl.hepg2[i,"hepg2.max.scorediff"] <- scorediff
    raqtl.hepg2[i,"tf.match"] <- tfbs.raqtl.hepg2[tfbs.raqtl.hepg2$snp.id == snp.idx, tfbs.match][1]
    
    # if the maximum value of the sure expression is in the ref allele AND the maximum
    # difference in tfbs score is negative, this means that the sure expression is in 
    # concordance with the tfbs scores
    
    if (raqtl.hepg2[i, "hepg2.max"] == "ref" & sign(raqtl.hepg2[i, "hepg2.max.scorediff"]) == -1) { raqtl.hepg2[i,"motif.concordance"] <- "concordance"} 
    if (raqtl.hepg2[i, "hepg2.max"] == "ref" & sign(raqtl.hepg2[i, "hepg2.max.scorediff"]) == 1)  { raqtl.hepg2[i,"motif.concordance"] <- "non.concordance"} 
    if (raqtl.hepg2[i, "hepg2.max"] == "alt" & sign(raqtl.hepg2[i, "hepg2.max.scorediff"]) == -1) { raqtl.hepg2[i,"motif.concordance"] <- "non.concordance"} 
    if (raqtl.hepg2[i, "hepg2.max"] == "alt" & sign(raqtl.hepg2[i, "hepg2.max.scorediff"]) == 1)  { raqtl.hepg2[i,"motif.concordance"] <- "concordance"}
    if (sign(raqtl.hepg2[i, "hepg2.max.scorediff"]) == 0 )                                      { raqtl.hepg2[i,"motif.concordance"] <- "undetermined"} 
    
    
    
  } else {
    
    raqtl.hepg2[i, "motif.concordance"] <- "no motifs"
  }
  
  
  
  
}
