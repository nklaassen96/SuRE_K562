DATE="date +'%Y%m%d_%H.%M.%S'"
DIR="~/projects/SuRE_K562/tmp/20190925_OutputTest.Rout"
CONFIGFILE=paste <($DATE) <($DIR)

/usr/bin/time -v nice -19 R --no-save -q < ~/projects/SuRE_K562/analyses/20190925_Indel_ImportNormalizeReformat.R > ~/projects/SuRE_K562/tmp/20190925_OutputTest.Rout