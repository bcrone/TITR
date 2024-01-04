#!/bin/bash

TISSUES=("ADIPOSE_TISSUE" "ADRENAL_GLAND" "ARTERIAL_BLOOD_VESSEL" "BLOOD" "BLOOD_VESSEL" "BONE_ELEMENT" "BONE_MARROW" "BRAIN" "BREAST" "COLON" "CONNECTIVE_TISSUE" "EAR" "EMBRYO" "ENDOCRINE_GLAND" "EPITHELIUM" "ESOPHAGUS" "EXOCRINE_GLAND" "EXTRAEMBRYONIC_COMPONENT" "EYE" "GENERIC" "GONAD" "HEART" "IMMUNE_ORGAN" "INTESTINE" "KIDNEY" "LARGE_INTESTINE" "LIMB" "LIVER" "LUNG" "LYMPHOID_TISSUE" "LYMPH_NODE" "MAMMARY_GLAND" "MOUTH" "MUSCULATURE_OF_BODY" "NERVE" "OVARY" "PANCREAS" "PENIS" "PLACENTA" "PROSTATE_GLAND" "SKIN_OF_BODY" "SKIN_OF_PREPUCE_OF_PENIS" "SMALL_INTESTINE" "SPINAL_CORD" "SPLEEN" "STOMACH" "TESTIS" "THYMUS" "THYROID_GLAND" "UTERUS" "VAGINA" "VASCULATURE")

TRAIT="$1"
ITERATION="$2"
ROOT_PATH="$3"

SUMSTATS="/path/to/GWAS/sumstats/${TRAIT}/LDSC/${TRAIT}.sumstats.gz"
OUTSTATS="/path/to/GWAS/sumstats/${TRAIT}/LDSC/${TRAIT}.ITERATION_${ITERATION}.sumstats.gz"
MASTER_SNPS="/path/to/results/${TRAIT}/${TRAIT}.TURF.ITERATION_$((ITERATION-1)).snplist"

rm slurm/${TRAIT}.run_ldsc.jobs
rm /path/to/results/${TRAIT}/${TRAIT}.terminate

gzip -cd ${SUMSTATS} | grep -v -w -f ${MASTER_SNPS} | gzip > ${OUTSTATS}

for TISSUE in ${TISSUES[@]}
do
	bash batch_run_ldsc.sh ${TRAIT} ${TISSUE} ${ITERATION} ${ROOT_PATH} > slurm/${TRAIT}.${TISSUE}.run_ldsc.sh
	RES=$(sbatch slurm/${TRAIT}.${TISSUE}.run_ldsc.sh)
	echo "${RES##* }" >> slurm/${TRAIT}.run_ldsc.jobs
done

DEPENDENCY="afterany"
while read -r line 
do
	DEPENDENCY="${DEPENDENCY}:${line}"
done < slurm/${TRAIT}.run_ldsc.jobs

sbatch --dependency=${DEPENDENCY} --account=apboyle99 --partition=standard --time=0:30:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=1g --job-name=link.run_ldsc.run_summarize_enrichment_tau_star \
--wrap "bash build_run_summarize_enrichment_tau_star.sh ${TRAIT} ${ITERATION} ${ROOT_PATH}"