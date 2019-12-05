/usr/bin/time -v nice -19 R --no-save -q < ~/projects/SuRE_K562/analyses/Analysis_Compare_Scripts_And_Pipeline/20191118_SNP_Indel_Import_NEWDOWNSAMPLE.R > ~/projects/SuRE_K562/analyses/Logs/20191118IndSnp_Import.Rout 2> ~/projects/SuRE_K562/analyses/Logs/20191118IndSnp_Import_ADJUST.Log;
/usr/bin/time -v nice -19 R --no-save -q < ~/projects/SuRE_K562/analyses/Analysis_Compare_Scripts_And_Pipeline/20191118_SNP_Indel_Wilcoxon_NEWDOWNSAMPLE.R > ~/projects/SuRE_K562/analyses/Logs/20191118IndSnp_Wilcox.Rout 2> ~/projects/SuRE_K562/analyses/Logs/20191118IndSnp_Import_ADJUST.Log;

