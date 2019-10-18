/usr/bin/time -v nice -19 R --no-save -q < ~/projects/SuRE_K562/tmp/20191710CONCATONATE_5_1.R > ~/projects/SuRE_K562/analyses/Logs/20191017scripts_1.Rout;
/usr/bin/time -v nice -19 R --no-save -q < ~/projects/SuRE_K562/analyses/20191016_SplitNormalizedCounts.R > ~/projects/SuRE_K562/analyses/Logs/20191017scripts_2.Rout;
/usr/bin/time -v nice -19 R --no-save -q < ~/projects/SuRE_K562/analyses/20190930_IndelWilcoxon.R > ~/projects/SuRE_K562/analyses/Logs/20191017scripts_3.Rout

