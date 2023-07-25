#!/bin/bash

conda activate ldsc

ROOT_PATH="$1"
TISSUE="$2"
CHROM="$3"

LDSC_PATH="${ROOT_PATH}/tools/ldsc"
KG_PATH="${ROOT_PATH}/data/1KG/1000G_EUR_Phase3_plink"
HM3_PATH="${ROOT_PATH}/data/1KG/hapmap3_snps"

TISSUE_PATH="${ROOT_PATH}/data/1KG/RegDB_custom_baselineLD_v1.2/serial"

echo "Starting LDscore regression"

echo "Calculating LDscores for ${TISSUE} chr${CHROM}"
PREFIX="customized_baselineLD_cts.RegDB.${TISSUE}.${CHROM}"
ANNOT="${PREFIX}.annot.gz"
HM3="${HM3_PATH}/hm.${CHROM}.snp"
python ${LDSC_PATH}/ldsc.py --l2 --bfile ${KG_PATH}/1000G.EUR.QC.RegDB.${CHROM} --ld-wind-cm 1 --out ${TISSUE_PATH}/${PREFIX} --annot ${TISSUE_PATH}/${ANNOT} --print-snps ${HM3} > ${TISSUE_PATH}/${PREFIX}.out

echo "LDscore complete for ${TISSUE}"