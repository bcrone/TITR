ROOT_PATH="$1"
RESULTS_PATH="$2"
ANCESTRY="$3"
TRAIT="$4"
TISSUE="$5"
PLINK_PATH="$6"
PLINK2_PATH="$7"

PARTITIONS=( 10 50 100 200 500 )

EUR_1KG_BFILE="${ROOT_PATH}/data/1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC"
SUMSTATS="${ROOT_PATH}/data/GWAS/${TRAIT}/${TRAIT}.PLINK.TITR"
RANGE_LIST="${ROOT_PATH}/data/UKB/phenos/range_list.expanded"

BED_PATH="${ROOT_PATH}/data/UKB/${ANCESTRY}/ukb_imp_chrALL_v3.bed"
BIM_PATH="${ROOT_PATH}/data/UKB/${ANCESTRY}/ukb_imp_chrALL_v3.bim"
FAM_PATH="${ROOT_PATH}/data/UKB/phenos/${ANCESTRY}/${ANCESTRY}.${TRAIT}.fam"
SAMPLES_PATH="${ROOT_PATH}/data/UKB/phenos/${ANCESTRY}/${ANCESTRY}.${TRAIT}.sample.IDs"

#Standard
${PLINK_PATH} --bfile ${EUR_1KG_BFILE} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump ${SUMSTATS} --clump-snp-field SNP --clump-field P --out ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.standard.clump

${PLINK2_PATH} --pgen ${BED_PATH} --pvar ${BIM_PATH} --psam ${FAM_PATH} --score ${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range ${RANGE_LIST} ${SUMSTATS}.SNP.pvalue --extract ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.standard.clump.clumped --keep ${SAMPLES_PATH} --out ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.standard

for PARTITION in ${PARTITIONS[@]}
do
	#IMPACT
	${PLINK_PATH} --bfile ${EUR_1KG_BFILE} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump ${SUMSTATS} --clump-snp-field SNP --clump-field P --extract ${IMPACT_PATH}/IMPACT707_EUR_chrALL.${PARTITION}.${TRAIT}.SNPs --out ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.IMPACT.${PARTITION}.clump

	${PLINK2_PATH} --pgen ${BED_PATH} --pvar ${BIM_PATH} --psam ${FAM_PATH} --score ${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range ${RANGE_LIST} ${SUMSTATS}.SNP.pvalue --extract ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.IMPACT.${PARTITION}.clump.clumped --keep ${SAMPLES_PATH} --out ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.IMPACT.${PARTITION}

	#SURF
	${PLINK_PATH} --bfile ${EUR_1KG_BFILE} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump ${SUMSTATS} --clump-snp-field SNP --clump-field P --extract ${SURF_PATH}/1000G_phase3_master_scores_1KG-pruned_quantile_normalized_chrALL.GENERIC.${PARTITION}.SNPs --out ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.SURF.${PARTITION}.clump

	${PLINK2_PATH} --pgen ${BED_PATH} --pvar ${BIM_PATH} --psam ${FAM_PATH} --score ${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range ${RANGE_LIST} ${SUMSTATS}.SNP.pvalue --extract ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.SURF.${PARTITION}.clump.clumped --keep ${SAMPLES_PATH} --out ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.SURF.${PARTITION}

	#TURF
	${PLINK_PATH} --bfile ${EUR_1KG_BFILE} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump ${SUMSTATS} --clump-snp-field SNP --clump-field P --extract ${TURF_PATH}/1000G_phase3_master_scores_1KG-pruned_quantile_normalized_chrALL.${TISSUE}.${PARTITION}.SNPs --out ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.TURF.${PARTITION}.clump

	${PLINK2_PATH} --pgen ${BED_PATH} --pvar ${BIM_PATH} --psam ${FAM_PATH} --score ${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range ${RANGE_LIST} ${SUMSTATS}.SNP.pvalue --extract ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.TURF.${PARTITION}.clump.clumped --keep ${SAMPLES_PATH} --out ${RESULTS_PATH}/${TRAIT}.${ANCESTRY}.TURF.${PARTITION}
done

touch ${RESULTS_PATH}/${TRAIT}.OK