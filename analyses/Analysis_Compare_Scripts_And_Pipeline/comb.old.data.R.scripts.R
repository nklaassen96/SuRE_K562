library(VennDiagram)
# This script is used to compare my R scripts with the R scripts from Joris

# For this we have two different datasets:
# 1. OLD: made with the old pipeline (NO equal, maternal, paternal etc.) and the R scripts from Joris. 
#         Contains the giant p-value table out of which we can subset the exact amount of raQTLs
# 2. NEW: made with the old pipeline (NO equal, maternal, paternal etc.) and the R scripts from myself (Noud)
#         Contains the giant p-value table. 




## Make the OLD input into raQTLs for K562
  old.pvalue.df.input <- fread(input = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/external/SuRE_OSF_NatGen/download.txt.gz")
  
  # Change colnames
  colnames(old.pvalue.df.input)[6] <- "K562.cDNA.ref.mean"
  colnames(old.pvalue.df.input)[7]<- "K562.cDNA.alt.mean"
  colnames(old.pvalue.df.input)[8] <- "K562.wilcoxon.pvalue"
  colnames(old.pvalue.df.input)[9] <- "K562.wilcoxon.pvalue.random"
  colnames(old.pvalue.df.input)[12] <- "HepG2.wilcoxon.pvalue"
  colnames(old.pvalue.df.input)[13] <- "HepG2.wilcoxon.pvalue.random"
  
    # min.max elements have already been selected for this file. So all we need to do is require the max sure signal to be higher than 4
  old.p.df <- old.pvalue.df.input[old.pvalue.df.input$K562.cDNA.ref.mean >4 | old.pvalue.df.input$K562.cDNA.alt.mean >4,] # 170,118 SNPs -> equal to paper
  old.raqtl.k562 <- old.p.df[old.p.df$K562.wilcoxon.pvalue <= 0.006192715,] # 19,236 K562 raQTL -> equal to paper (-1)

## Make the NEW input into raQTLs for K562
  new.pvalue.df.input <- readRDS(file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Test_Old_Data/pvalue.dataframe.all.20191106.RDS")
  
  new.p.df <- new.pvalue.df.input[new.pvalue.df.input$max.elements <1000 & new.pvalue.df.input$min.elements >= 10,] 
  
  # check what the differences are at this point
  sum(!old.pvalue.df.input$SNP_ID %in% new.p.df$SNP_ID) # 9238 SNPs (out of 5,919,239) are present in OLD but not in NEW p.value df (before calculating raQTLs)
  sum(!new.p.df$SNP_ID %in% old.pvalue.df.input$SNP_ID) # 2401 SNPs (out of 5,912,396) are present in NEW but not in OLD p.value df (before calculating raQTLs)
  
  new.p.df <- new.p.df[new.p.df$max.k562.expression >4,]  
  
  # make raQTL 
  {
  
  k562.p.real <- sort(-log10(new.p.df$K562.wilcoxon.pvalue))  
  k562.p.shuf <- sort(-log10(new.p.df$K562.wilcoxon.pvalue.random))
  k562.FDR.values <- NULL
  k562.p.values <- NULL
  k562.sign.real.values <- NULL
  k562.sign.shuf.values <- NULL
  
  # loop through a series of p value threshold to determines at what p value threshold the FDR
  # is around 0.05
  
  
  
  for (p in 10^-seq(2,3,by=0.01)) {
    
    # generate the -log10(p) value
    print(p)
    p.threshold <- -log10(p)
    
    # count the amount of p values that would fall within the p value threshold
    
    sign.real <- sum(k562.p.real > p.threshold)
    sign.shuf <- sum(k562.p.shuf > p.threshold)
    
    # calculate the FDR (fraction of shuffled p values that fall within the p value threshold)
    
    FDR <- round(sign.shuf / (sign.real + sign.shuf), digits = 3)
    
    # summarize the results in 4 vectors
    
    k562.FDR.values <- c(k562.FDR.values, FDR)
    k562.p.values <- c(k562.p.values, p)
    k562.sign.real.values <- c(k562.sign.real.values, sign.real)
    k562.sign.shuf.values <- c(k562.sign.shuf.values, sign.shuf)
    
    k562.df <- data.frame(k562.p.values, k562.sign.real.values, k562.sign.shuf.values, k562.FDR.values)
    colnames(k562.df) <- c("p value threshold", "nr. of significant 'real' p values", "nr. of significant shuffled p values", "FDR")
  }
  
  k562.fdr.pvalue <- k562.df[k562.df$FDR == 0.050,]$`p value threshold`[1] #find the first value for which FDR = 5% (it could be possible that this value is not there)
  new.raqtl.k562 <- new.p.df[new.p.df$K562.wilcoxon.pvalue < k562.fdr.pvalue,]
  new.raqtl.k562 <- new.raqtl.k562[!is.na(new.raqtl.k562$K562.wilcoxon.pvalue),]
  }
  
## Generate one big OLD + NEW dataframe
  
  # remove 7 duplicated values
  old.pvalue.df.input <- old.pvalue.df.input[!duplicated(old.pvalue.df.input$SNP_ID),]
  
  # sort both dataframes on SNP_ID
  old.pvalue.df.input <- old.pvalue.df.input[order(old.pvalue.df.input$SNP_ID),]
  new.pvalue.df.input <- new.pvalue.df.input[order(new.pvalue.df.input$SNP_ID),]
  
  # find the overlapping snp.ids
  pub.in.nov.idx <- which(old.pvalue.df.input$SNP_ID %in% new.pvalue.df.input$SNP_ID)
  nov.in.pub.idx <- which(new.pvalue.df.input$SNP_ID %in% old.pvalue.df.input$SNP_ID)

  # combine the snp.ids into one big dataframe
  all.comb.dt <- cbind(old.pvalue.df.input[pub.in.nov.idx,], new.pvalue.df.input[nov.in.pub.idx,])
  
  #check
  sum(all.comb.dt[[2]] != all.comb.dt[[16]]) == 0
  
## Generate figures for the combined dataframe
  
  # generate a figure directory
  
  fig.dir <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/processed/Figures/Compare_Rscripts/"
  
  sample.vector <- sample(seq(1:nrow(all.comb.dt)), size = 100000)
  
  # -log10(p)
  png(filename = paste0(fig.dir, "All.pvalue.png"))
  plot(x = -log10(all.comb.dt[[8]])[sample.vector], y = -log10(all.comb.dt[[31]])[sample.vector], cex = 0.1, xlab = "Published -log10(p-value)", ylab = "Novel -log10(p-value)", main = paste0("All tested SNPs (sampled n=", length(sample.vector),")"), sub = "hi") 
  abline(a = 0, b = 1)
  dev.off()
  
  # ref.element.count
  png(filename = paste0(fig.dir, "All.ref_element_count.png"))
  plot(x = all.comb.dt[[4]][sample.vector], y = all.comb.dt[[21]][sample.vector], cex = 0.1, xlab = "Published ref_element_count", ylab = "Novel ref_element_count", main = paste0("All tested SNPs (sampled n=", length(sample.vector),")"))
  abline(a = 0, b = 1)
  dev.off()
  
  # alt.element.count
  png(filename = paste0(fig.dir, "All.alt_element_count.png"))
  plot(x = all.comb.dt[[5]][sample.vector], y = all.comb.dt[[22]][sample.vector], cex = 0.1, xlab = "Published alt_element_count", ylab = "Novel alt_element_count", main = "Alt Element count")
  abline(a = 0, b = 1)
  dev.off()
  
  
  # K562.ref.mean
  png(filename = paste0(fig.dir, "All.ref_cDNA_mean.png"))
  plot(x = all.comb.dt[[6]][sample.vector], y = all.comb.dt[[23]][sample.vector], cex = 0.1, xlab = "Published ref_cDNA_mean", ylab = "Novel ref_cDNA_mean", main = "ref_cDNA_mean",
       )
  abline(a = 0, b = 1)
  
  dev.off()
  
  # K562.alt.mean
  png(filename = paste0(fig.dir, "All.alt_cDNA_mean.png"))
  plot(x = all.comb.dt[[7]][sample.vector], y = all.comb.dt[[24]][sample.vector], cex = 0.1, xlab = "Published alt_cDNA_mean", ylab = "Novel alt_cDNA_mean", main = "alt_cDNA_mean")
  abline(a = 0, b = 1)
  dev.off()
  
  ## CDNA STRANGE VALUE DISCOVERY
  
    # I notice that there seems to be a "second slope" in the ref.cDNA.mean and alt.cDNA.mean
    # I want to light out one specific example and see why this is the case
    
    sure.reads.allrep.chr.17 <- readRDS("~/projects/SuRE_K562/data/interim/SuRE_Test_Old_Data/sure.reads.allrep.chr.17.RDS")
    sure.reads.allrep.chr.7 <- readRDS("~/projects/SuRE_K562/data/interim/SuRE_Test_Old_Data/sure.reads.allrep.chr.7.RDS")
    
    # First I select all snps that have a difference bigger than 5 in ref.cdna.mean
    tst.idx <- which((all.comb.dt[[23]]-all.comb.dt[[6]]) > 5)
    tst.df <- all.comb.dt[tst.idx,]
    tst.df[[23]]-tst.df[[6]]
    
    # select a snp that you want to figure out a bit more
    tst.snp.id <- tst.df[8,]$SNP_ID 
    
    tst.1 <- sure.reads.allrep.chr.17[sure.reads.allrep.chr.17$SNP_ID == tst.snp.id,]
    tst.2 <- sure.reads.allrep.chr.7[sure.reads.allrep.chr.7$SNP_ID == tst.snp.id,]
    
    # add a column to specify the library and variant in one column
    tst.1$lib.var <- paste0(tst.1$SNPvarInf,".", tst.1$library)
    tst.2$lib.var <- paste0(tst.2$SNPvarInf,".", tst.2$library)
    
    # calculate the mean per library per variant (ref or alt)
    tapply(tst.1$cDNA.K562.norm.ipcr, tst.1$lib.var, mean)
    tapply(tst.2$cDNA.K562.norm.ipcr, tst.2$lib.var, mean)
  
  
  # draw overlap in raqtl SNP_IDs
  
  draw.pairwise.venn(area1 = n_distinct(new.raqtl.k562$SNP_ID),
                     area2 = n_distinct(old.raqtl.k562$SNP_ID), 
                     cross.area = sum(new.raqtl.k562$SNP_ID %in% old.raqtl.k562$SNP_ID),
                     category = c("NEW", "OLD"), fill = c(2,3),cat.cex = 2 , cat.pos = c(-40,40), cex = c(rep(1.4, 3)))
  
  
  # Find the 1043 SNPs that I cannot find in my analysis
  
  unfindable.idx <- which(!old.raqtl.k562$SNP_ID %in% new.raqtl.k562$SNP_ID)
  unfindable.snp.ids <- old.raqtl.k562[unfindable.idx,]$SNP_ID
  
  # change the colnames of the new parameters to new.parameter
  unfindable.comb.dt <- all.comb.dt[all.comb.dt$SNP_ID %in% unfindable.snp.ids,]
  colnames(unfindable.comb.dt)[16:41] <- paste0("new.", colnames(all.comb.dt)[16:41])
  
  # HIST ref.element.count
  foo <- hist(unfindable.comb.dt$ref.element.count-unfindable.comb.dt$new.ref.element.count, 
              breaks = 20, 
              col = 3, 
              xaxt="n", main = "'Unfindable SNPs' (n=1043)", xlab = "Difference in reference element count (old - new pipeline)")
  axis(side = 1, at = c(foo$mids[1]-1,foo$mids)  ,labels = c(foo$breaks))
  
  # HIST alt.element.count
  foo <-hist(unfindable.comb.dt$alt.element.count-unfindable.comb.dt$new.alt.element.count, 
             breaks = 10, 
             col = 3,
             xaxt="n", 
             main = "'Unfindable SNPs' (n=1043)", xlab = "Difference in alternative element count (old - new pipeline)")
  axis(side = 1, at = c(foo$mids[1]-1,foo$mids)  ,labels = c(foo$breaks))

  #hist(-log10(unfindable.comb.dt$K562.wilcoxon.pvalue)--log10(unfindable.comb.dt$new.K562.wilcoxon.pvalue), breaks = 50, col = 3)
  
  # HIST ref.mean
  foo <- hist(unfindable.comb.dt$K562.cDNA.ref.mean-unfindable.comb.dt$new.K562.cDNA.ref.mean, 
         breaks = 50, 
         col = 3, xlim = c(-3,3), xaxt = "n", main = "'Unfindable SNPs' (n=1043)", xlab = "Difference in reference mean (old - new pipeline)")
  axis(side = 1, at = c(foo$mids[1]-1,foo$mids)  ,labels = c(foo$breaks))
  
  # HIST alt.mean
  foo <- hist(unfindable.comb.dt$K562.cDNA.alt.mean-unfindable.comb.dt$new.K562.cDNA.alt.mean, 
              breaks = 100, 
              col = 3, 
              xlim = c(-5,5))
  
  
  
  samp <- sample(new.p.df$K562.wilcoxon.pvalue, size = 1043)
  #plot all pvalues
  hist(-log10(samp), breaks = 100, col = alpha("gray", 0.5), ylim = c(0,650), xlab = "-log10(p-value)", main = "Histogram pvalue 'unfindable SNPs'")
  #plot unfindable new.pvalues
  hist(-log10(unfindable.comb.dt$new.K562.wilcoxon.pvalue), add = T, breaks = 50, col = alpha("red", 0.5))
  abline(v = -log10(k562.fdr.pvalue), col = 4, lwd = 3, lty = 2)
  legend("topright", fill = c("gray", "red"), c("All tested SNPs (sampled 1043)", "published raQTL, new analysis no raQTL (n=1043)"))
  text(x = 5, y = 200, labels = "p-value threshold", col = 4, cex = 1.5)
  
  
  
  #histogram novel p values
  
  nov.raqtl_pub.NA.idx <- which(!novel.raqtl.k$SNP_ID %in% published$SNP_ID)
  nov.raqtl_pub.NA <- novel.raqtl.k[nov.raqtl_pub.NA.idx,]
  
  nov.raqtl_pub.not.raqtl.idx <- which(!novel.raqtl.k$SNP_ID %in% pub.raqtl.k$SNP_ID)
  nov.raqtl_pub.not.raqtl <- novel.raqtl.k[nov.raqtl_pub.not.raqtl.idx,]
  
  nov.raqtl_pub.raqtl.idx <- which(novel.raqtl.k$SNP_ID %in% pub.raqtl.k$SNP_ID)
  nov.raqtl_pub.raqtl <- novel.raqtl.k[nov.raqtl_pub.raqtl.idx,]
  
  # plot all nov.raqtl
  #hist(-log10(novel.raqtl.k$K562.wilcoxon.pvalue), breaks= 3000, xlim = c(2,10), col = 2, main = "novel raQTL K562 (SNP only)", xlab = "-log10(p.value)")
  
  # plot "novel raqtl but not a raqtl in published" only
  hist(-log10(nov.raqtl_pub.not.raqtl$K562.wilcoxon.pvalue),xlim = c(2,10), breaks = 1000, add = F, col = 3,main = paste0("novel K562 raQTL  (SNP only) (n=", nrow(novel.raqtl.k),")"), xlab = "-log10(p.value)")
  
  # plot "novel raqtl and also raqtl in published" only
  hist(-log10(nov.raqtl_pub.raqtl$K562.wilcoxon.pvalue), breaks = 600, col = alpha(2,0.5), add = T)
  
  # plot "novel raqtl but not present in published" only
  hist(-log10(nov.raqtl_pub.NA$K562.wilcoxon.pvalue), breaks = 600, add = T, col = 6)
  
  legend("topright", c(paste0("Novel raQTL (no raQTL in published) (n=", length(nov.raqtl_pub.not.raqtl.idx), ")"),
                       paste0("Novel raQTL (raQTL in published) (n=", length(nov.raqtl_pub.raqtl$K562.wilcoxon.pvalue), ")"),
                       paste0("Novel raQTL (SNP_ID not in published) (n=", length(nov.raqtl_pub.NA$K562.wilcoxon.pvalue), ")")
  ),
  fill = c(3,alpha(2,0.5), 6))
  
  
  
  