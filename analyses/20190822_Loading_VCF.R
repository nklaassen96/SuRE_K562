#Reading the Indel and SNV data from the K562 paper
#The comment.char argument removes the lines followed by a "#" 

K562_IndelSNV_NF <- read.table("~/projects/SuRE_K562/data/external/Encode_K562_VCF/ENCFF606RIC.vcf.gz", comment.char = "#")
colnames(K562_IndelSNV_NF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", 'INFO', 'FORMAT', 'K562_PE')

#Make a short testing dataframe with only 1000 MNV
K562_IndelSNV_NF_short <- K562_IndelSNV_NF[c(1:1000),c(1:10)]
#Generate a histogram 
#hist(K562_IndelSNV_NF_short[,6], col = 2, xlab = "QUAL", breaks = 40)
hist(K562_IndelSNV_NF[,6], col = 3, xlab = "QUAL", xlim = c(0, 6000),breaks = 1000, main = "Variant Calling Quality")

Count = c()
for (i in 1:23) {
  if (i == 23) {i <- c("X")}
        Chr_i = paste("chr", i, sep="")
        Chr_i_list <- K562_IndelSNV_NF[which(K562_IndelSNV_NF$CHROM == Chr_i),]
        MNV_Count <- nrow(Chr_i_list)
  if (i == "X") {i <- 23}
        Count[i] <- MNV_Count
} #Counting MNV in chromosome 1 till 22 plus X

#plot the nmbr of MNVs per chromosome
barplot(Count,space = 0.2,names.arg = c(1:23), col = 2, xlab = "Chromosome", ylab = "# of Multi-nucleotide Variants (MNV)")

#Generate a dataframe of SNVs (removing indels) #does not work as it does not account for multiple ALT sequences
Nt_REF <- nchar(as.character(K562_IndelSNV_NF[,4])) #Count nr. of characters of REF in a vector
SNV_REF_Only <- K562_IndelSNV_NF[which(Nt_REF == 1),] #Filter in the database for REF==1
Nt_ALT <- nchar(as.character(SNV_REF_Only[,5])) #Count in the new database nr of characters in ALT
K562_SNV <- SNV_REF_Only[which(Nt_ALT == 1),]
K562_MUT_SNV <- K562_SNV[which(nchar(as.character(K562_SNV[,3])) == 1),] #filter for rs numbers in single-nucleotide variations
K562_MUT_MNV <- K562_IndelSNV_NF[which(nchar(as.character(K562_IndelSNV_NF[,3])) == 1),] #filter for rs numbers in all variations


#Get a dataframe of SNVs or indels that appear at the same location
df <- K562_IndelSNV_NF
df$CHROMPOS <- paste(df$CHROM, df$POS)
n_occur <- data.frame(table(df$CHROMPOS))
uni <- n_occur[n_occur$Freq > 1,]

#install data.table
install.packages("data.table")
library(data.table)
library(dplyr)
#Get a dataframe of SNVs or indels that overlap eachother
#in accordance with https://stackoverflow.com/questions/40129485/overlap-ranges-in-single-dataframe
df2 <- K562_IndelSNV_NF[1:1000,]
df2$MAX <- df2$POS+nchar(as.character(df2$REF))
c <- outer(df2$MAX, df2$POS, ">")
d <- outer(df2$POS, df2$MAX, "<")
df2o <- df2 %>%
  mutate(Overlap = apply(c & d, 1, sum) > 1
)
table(df2o$Overlap)
df2o <- df2o %>% filter(df2o$Overlap == TRUE) #returns overlapping mnvs
#works for 1000; works for 10,000; works for 50,000 (5 min); does not work for 3,800,000

#Check ALT sequences that contain a "," indicating multiple alternative sequences
tst <- K562_IndelSNV_NF[grep(",", as.character(K562_IndelSNV_NF$ALT)),]
tst2 <- tst[which(nchar(as.character(tst$ALT)) == 3),]
