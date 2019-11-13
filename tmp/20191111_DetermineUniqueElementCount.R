
library(data.table)
library(tidyverse)
library(foreach)
library(doMC)
registerDoMC(cores = 40)

#### OLD unique.element.count ####




rep1.files.vector <- paste0("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE42/SuRE-pipelineOutput/", list.files("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE42/SuRE-pipelineOutput/",pattern = "inf"))[c(1:22,24)]
rep3.files.vector <- paste0("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE43/SuRE43-pipelineOutput/",list.files("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE43/SuRE43-pipelineOutput/", pattern = "inf"))[c(1:22,24)]
rep5.files.vector <- paste0("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE44/SuRE44-pipelineOutput/",list.files("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE44/SuRE44-pipelineOutput/", pattern = "inf"))[c(1:22,24)]
rep7.files.vector <- paste0("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE45_1/SuRE45_1-pipelineOutput/",list.files("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE45_1/SuRE45_1-pipelineOutput/", pattern = "inf"))[c(1:22,24)]
rep8.files.vector <- paste0("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE45_2/SuRE45_2-pipelineOutput/",list.files("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE45_2/SuRE45_2-pipelineOutput/", pattern = "inf"))[c(1:22,24)]
rep2.files.vector <- paste0("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE42_2/pipelineOutput/", list.files("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE42_2/pipelineOutput/", pattern = "txt.gz"))
rep4.files.vector <- paste0("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE43_2/pipelineOutput/", list.files("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE43_2/pipelineOutput/", pattern = "txt.gz"))
rep6.files.vector <- paste0("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE44_2/pipelineOutput/", list.files("/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE44_2/pipelineOutput/", pattern = "txt.gz"))

file.list <- list(rep1.files.vector,
                  rep2.files.vector,
                  rep3.files.vector,
                  rep4.files.vector,
                  rep5.files.vector,
                  rep6.files.vector,
                  rep7.files.vector,
                  rep8.files.vector)

rep.vector <- c("42_1","42_2","43_1","43_2","44_1","44_2","45_1","45_2")

unique.elements.old <- for(rep.idx in c(1:8)){  
  
  # retrieve all files in a vector for this replicate
  file.vector <- unlist(file.list[rep.idx])
  
  # Genarate a string that can be used to name and identify files
  rep.str <- rep.vector[rep.idx]
  
  
  
  # First step is to normalize the data per replicate, to do that,
  # we need the total amount of ipcr counts and cdna counts. The
  # foreach loops provide a simple way to count the sum of the aforementioned
  # variables
  
  #check whether the normalizing file is already there, otherwise, generate it
  
  read.totals.file <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Test_Old_Data/unique.element.count.rep.",rep.idx,".RDS")
  
    unique.element.count <- sum(unlist(
      foreach(k = c(1:length(file.vector))) %dopar% {
        
        nrow(fread(file.vector[k], 
                  header = TRUE, 
                  sep = "\t", 
                  stringsAsFactors = FALSE, 
                  select = c("iPCR")))}
      ))
    
    assign(x = paste0("unique.element.count.", rep.str), value = unique.element.count)
    #saveRDS(object = unique.element.count, file = read.totals.file)
}

  


#### NEW unique.element.count ####
library(data.table)
library(tidyverse)
library(foreach)
library(doMC)
registerDoMC(cores = 8)

files.vector <- rep(list.files("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE42-1/equal/"), each = 3)

rep1.files.vector <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE42-1/",c("equal/", "maternal/", "paternal/"),files.vector)
rep2.files.vector <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE42-2/",c("equal/", "maternal/", "paternal/"),files.vector)
rep3.files.vector <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE43-1/",c("equal/", "maternal/", "paternal/"),files.vector)
rep4.files.vector <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE43-2/",c("equal/", "maternal/", "paternal/"),files.vector)
rep5.files.vector <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE44-1/",c("equal/", "maternal/", "paternal/"),files.vector)
rep6.files.vector <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE44-2/",c("equal/", "maternal/", "paternal/"),files.vector)
rep7.files.vector <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE45-1/",c("equal/", "maternal/", "paternal/"),files.vector)
rep8.files.vector <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/SuRE45-2/",c("equal/", "maternal/", "paternal/"),files.vector)

file.list <- list(rep1.files.vector,
                  rep2.files.vector,
                  rep3.files.vector,
                  rep4.files.vector,
                  rep5.files.vector,
                  rep6.files.vector,
                  rep7.files.vector,
                  rep8.files.vector)

rep.vector <- c("42_1","42_2","43_1","43_2","44_1","44_2","45_1","45_2")

unique.elements.old <- foreach(rep.idx = c(1:8)) %dopar% {  
  
  # retrieve all files in a vector for this replicate
  file.vector <- unlist(file.list[rep.idx])
  
  # Genarate a string that can be used to name and identify files
  rep.str <- rep.vector[rep.idx]
  
  
  
  # First step is to normalize the data per replicate, to do that,
  # we need the total amount of ipcr counts and cdna counts. The
  # foreach loops provide a simple way to count the sum of the aforementioned
  # variables
  
  #check whether the normalizing file is already there, otherwise, generate it
  
  read.totals.file <- paste0("/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Test_Old_Data/unique.element.count.rep.",rep.idx,".RDS")
  
  sum(unlist(
    foreach(k = c(1:length(file.vector))) %dopar% {nrow(fread(file.vector[k],header = TRUE, sep = "\t",stringsAsFactors = FALSE,select = 1))} # 2 = length(file.vector)
  ))
  
  #assign(x = paste0("unique.element.count.", rep.str), value = unique.element.count)
  #saveRDS(object = unique.element.count, file = read.totals.file)
}

saveRDS(object = unique.elements.old, file = "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/unique.element.count.allreps.RDS")
