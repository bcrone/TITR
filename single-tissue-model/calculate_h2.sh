#!/bin/bash

conda activate ldsc

ROOT_PATH="$1"
TRAIT="$2"
TISSUE="$3"

LDSC_PATH="${ROOT_PATH}/tools/ldsc"
WEIGHTS_PREFIX="${ROOT_PATH}/data/1KG/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
FREQ_PREFIX="${ROOT_PATH}/data/1KG/1000G_Phase3_frq/1000G.EUR.QC."
LDSCORE_PREFIX="${ROOT_PATH}/data/1KG/RegDB_custom_baselineLD_v1.2/serial/customized_baselineLD_cts.RegDB.quantile_normalized.${TISSUE}."
SUMSTATS="${ROOT_PATH}/data/GWAS/${TRAIT}/${TRAIT}.sumstats.gz"

OUTPUT_PREFIX="${ROOT_PATH}/results/h2/${TRAIT}/${TRAIT}/${TISSUE}"

echo "Calculating h2 estimates for ${TRAIT} ${TISSUE}"

python ${LDSC_PATH}/ldsc.py --h2 ${SUMSTATS} --ref-ld-chr ${LDSCORE_PREFIX} --w-ld-chr ${WEIGHTS_PREFIX} --overlap-annot --frqfile-chr ${FREQ_PREFIX} --out ${OUTPUT_PREFIX} --print-coefficients --print-delete-vals > ${OUTPUT_PREFIX}.out

echo "Completed h2 estimates for ${TRAIT} ${TISSUE}"