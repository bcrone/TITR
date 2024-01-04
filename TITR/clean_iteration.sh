#!/bin/bash

TISSUES=("ADIPOSE_TISSUE" "ADRENAL_GLAND" "ARTERIAL_BLOOD_VESSEL" "BLOOD" "BLOOD_VESSEL" "BONE_ELEMENT" "BONE_MARROW" "BRAIN" "BREAST" "COLON" "CONNECTIVE_TISSUE" "EAR" "EMBRYO" "ENDOCRINE_GLAND" "EPITHELIUM" "ESOPHAGUS" "EXOCRINE_GLAND" "EXTRAEMBRYONIC_COMPONENT" "EYE" "GENERIC" "GONAD" "HEART" "IMMUNE_ORGAN" "INTESTINE" "KIDNEY" "LARGE_INTESTINE" "LIMB" "LIVER" "LUNG" "LYMPHOID_TISSUE" "LYMPH_NODE" "MAMMARY_GLAND" "MOUTH" "MUSCULATURE_OF_BODY" "NERVE" "OVARY" "PANCREAS" "PENIS" "PLACENTA" "PROSTATE_GLAND" "SKIN_OF_BODY" "SKIN_OF_PREPUCE_OF_PENIS" "SMALL_INTESTINE" "SPINAL_CORD" "SPLEEN" "STOMACH" "TESTIS" "THYMUS" "THYROID_GLAND" "UTERUS" "VAGINA" "VASCULATURE")

TRAIT="$1"
ITERATION="$2"
ROOT_PATH="$3"

if [ -f ${ROOT_PATH}/results/traits/${TRAIT}/${TRAIT}.terminate ]
then
	exit 0
fi

for TISSUE in ${TISSUES[@]}
do
	rm ${ROOT_PATH}/results/${TISSUE}/${TRAIT}.${TISSUE}.*.ITERATION_${ITERATION}.*.sscore
	rm ${ROOT_PATH}/results/${TISSUE}/${TRAIT}.${TISSUE}.*.ITERATION_${ITERATION}.*.log
	rm ${ROOT_PATH}/results/${TISSUE}/${TRAIT}.${TISSUE}.*.ITERATION_${ITERATION}.*.SNPs 
	rm ${ROOT_PATH}/results/${TISSUE}/${TRAIT}.${TISSUE}.*.ITERATION_${ITERATION}.*.snplist

rm ${ROOT_PATH}/repo/RegDB-tissue-heritability/PRS/*out
rm ${ROOT_PATH}/repo/RegDB-tissue-heritability/PRS/logs/${TRAIT}.*.run_ldsc.log
rm ${ROOT_PATH}/repo/RegDB-tissue-heritability/PRS/logs/${TRAIT}.*.run_score.log
done

sbatch --account=apboyle99 --partition=standard --time=0:30:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=1g --job-name=link.run_score.run_optimize_TURF_PRS \
--wrap "bash build_run_ldsc.sh ${TRAIT} $((ITERATION+1)) ${ROOT_PATH}"