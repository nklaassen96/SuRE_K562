
dir <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/SuRE_output"
files <- list.files(dir)

sure = NULL
for (chrom in paste0("chr", c(1:13,15:18, 21:22, "X"))){
  table <- read.table(paste0(dir, "/", "SuRE-counts_", chrom, ".txt.gz"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sure <- rbind(sure, table)
  print(paste0(chrom, " is finished"))
}

### Total (minus chr14 and chr19 and chr20)
length(unique(sure$SNP_ID))
length(unique(names(gr)))
length(unique(sure$SNP_ID))/length(unique(names(gr)))*100 # % coverage of IDs (not sure if i should take indels into account)
sureuni <- sure[!grepl(",", sure$SNPbase),]
sure
nrow(sure)-nrow(sureuni)
table(sureuni$SNPbase)
tst <- table(sureuni$SNPbase)
tot <- sum(tst[2:11])

sum(grepl("rs", sureuni$SNP_ID))
sum(grepl("chr", sureuni$SNP_ID))


df <- as.data.frame(table(sureuni$SNP_ID), stringsAsFactors = FALSE, skip = 1)
df <- df[2:nrow(df),]
summary(df$Freq)

dfrs <- df[grepl("rs", df$Var1),]
dfchr <- df[grepl("chr", df$Var1),]
summary(dfrs$Freq)
summary(dfchr$Freq)
boxplot(dfrs$Freq, ylim = c(0,25))

hist(dfrs$Freq, xlim = c(0,25), breaks = 2000, col = alpha("red", 0.5))
hist(dfchr$Freq, xlim = c(0,25), breaks = 100, col = alpha("green", 0.5), add = F)

gr_sureuni <- GRanges(seqnames = sureuni$chr,
                   ranges = IRanges(start = sureuni$start, 
                                    end = sureuni$end), 
                   strand = sureuni$strand)


df <- as.data.frame(tapply(sureuni$SuRE23_45_B1_T1, sureuni$SNP_ID, sum), skip = 1) #SuRE-counts per SNP_ID (1st replicate only)
df <- df[2:nrow(df),]
sureuni$SuRE_B1 <- sureuni$SuRE23_45_B1_T1+sureuni$SuRE23_45_B1_T2+sureuni$SuRE23_45_B1_T3+sureuni$SuRE23_45_B1_T4
sureuni$SuRE_B2 <- sureuni$SuRE23_55_B2_T1+sureuni$SuRE23_55_B2_T2
sureuni$SuRE_B <- sureuni$SuRE_B1+sureuni$SuRE_B2
sureuni_tst <- sureuni[#sureuni$SNP_ID != "" & 
                       grepl("chr", sureuni$SNP_ID) &
                       sureuni$SNPvar != 3 &
                       sureuni$SuRE_B > 0,
                         ] #Data with: SNP_ID (rs or chr) normalized counts of >1 for both Biological replicates 

#Try to see if i can plot the surecount coverage over a part of chr1
gr_sure <- GRanges(seqnames = sure$chr,
                   ranges = IRanges(start = sure$start, 
                                    end = sure$end), 
                   strand = sure$strand)
gr_sub <- gr_sure[2700000:2900000] #100,000 - 200,000
gr_subtst <- gr_sure[100000:100100]

zoom <- toGRanges(data.frame("chr1", 149030000, 149040000))
#zoom <- toGRanges("chr1:140,000,000-160,050,000")
kp <- plotKaryotype(plot.type = 2, 
                    chromosomes = c("chr1"), 
                    main = "",
                    zoom = zoom)
#kpAddBaseNumbers(kp, tick.dist = 1000000,minor.tick.dist = 100000 )
kpAddBaseNumbers(kp)
kp <- kpPlotCoverage(kp, data = gr_sub,)
kp <- kpPlotRegions(kp, data = gr_sub, data.panel = 2, layer.margin = 0.0001)

