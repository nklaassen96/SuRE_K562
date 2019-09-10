DRYRUN=""
# DRYRUN="--dryrun "

DATETAG="NK$( date +"%y%m%d_%H%M" )"
SNAKEFILE=/DATA/usr/ludo/projects/LP190816_SuRE_K562_SNP/code/SuRE-WASP-pipeline/SuRE-snakemake
CONFIG=SuRE-K562.yaml
LOG="${CONFIG%.yml}_run-${DATETAG}.log"
TARGET="bedpe_merged_smpls"
TARGET="merged_ipcr_cdna"
TARGET="reversed_liftover"
TARGET="sorted_cnt_tbls"
TARGET="bed2coverage_done"
TARGET="all"

CMD="/usr/bin/time -v nice -19 snakemake "${DRYRUN}"-prs "${SNAKEFILE}" --use-conda --resources ram=150 --configfile "${CONFIG}" --cores 25 ${TARGET} &> "${LOG}""
echo "${CMD}"
eval ${CMD}

