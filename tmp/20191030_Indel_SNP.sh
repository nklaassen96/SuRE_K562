/usr/bin/time -v nice -19 R --no-save -q < ~/projects/SuRE_K562/analyses/20191018_SNP_Indel_Import.R > ~/projects/SuRE_K562/analyses/Logs/20191031IndSnp_Import.Rout 2> ~/projects/SuRE_K562/analyses/Logs/20191031IndSnp_Import.Log;
/usr/bin/time -v nice -19 R --no-save -q < ~/projects/SuRE_K562/analyses/20191018_SNP_Indel_Wilcoxon.R > ~/projects/SuRE_K562/analyses/Logs/20191031IndSnp_Wilcox.Rout 2> ~/projects/SuRE_K562/analyses/Logs/20191031IndSnp_Wilcox.Log

