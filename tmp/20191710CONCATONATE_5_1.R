
input.dir <- "/DATA/usr/n.klaassen/projects/SuRE_K562/data/interim/SuRE_Indels_gDNA_Count/"

x5 <- readRDS(file = paste0(input.dir, "sure.indel.5rep.RDS"))
x1 <- readRDS(file = paste0(input.dir, "SuRE44-2/combined/sure.44_2.counts.indel.all.RDS"))

x <- rbind(x1, x5)

saveRDS(object = x, file = paste0(input.dir, "sure.indel.6rep.RDS"))