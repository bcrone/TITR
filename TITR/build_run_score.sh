#!/bin/bash

TISSUES=("ADIPOSE_TISSUE" "ADRENAL_GLAND" "ARTERIAL_BLOOD_VESSEL" "BLOOD" "BLOOD_VESSEL" "BONE_ELEMENT" "BONE_MARROW" "BRAIN" "BREAST" "COLON" "CONNECTIVE_TISSUE" "EAR" "EMBRYO" "ENDOCRINE_GLAND" "EPITHELIUM" "ESOPHAGUS" "EXOCRINE_GLAND" "EXTRAEMBRYONIC_COMPONENT" "EYE" "GENERIC" "GONAD" "HEART" "IMMUNE_ORGAN" "INTESTINE" "KIDNEY" "LARGE_INTESTINE" "LIMB" "LIVER" "LUNG" "LYMPHOID_TISSUE" "LYMPH_NODE" "MAMMARY_GLAND" "MOUTH" "MUSCULATURE_OF_BODY" "NERVE" "OVARY" "PANCREAS" "PENIS" "PLACENTA" "PROSTATE_GLAND" "SKIN_OF_BODY" "SKIN_OF_PREPUCE_OF_PENIS" "SMALL_INTESTINE" "SPINAL_CORD" "SPLEEN" "STOMACH" "TESTIS" "THYMUS" "THYROID_GLAND" "UTERUS" "VAGINA" "VASCULATURE")

TRAIT="$1"
ITERATION="$2"
ROOT_PATH="$3"

rm slurm/${TRAIT}.run_score.jobs

if (($ITERATION==1))
then
	rm ${TRAIT}.partitions
    for TISSUE in ${TISSUES[@]}
    do
		echo "${TISSUE},1" >> ${TRAIT}.partitions
    done
fi

readarray -t TISSUES < ${TRAIT}.lead_tissue

while IFS=, read TISSUE PARTITION
do
	if [[ " ${TISSUES[*]} " =~ " ${TISSUE} " ]]
    then
		bash batch_run_score.sh ${TRAIT} ${TISSUE} ${ITERATION} ${PARTITION} ${ROOT_PATH} > slurm/${TRAIT}.${TISSUE}.run_score.sh
		RES=$(sbatch slurm/${TRAIT}.${TISSUE}.run_score.sh)
		echo "${RES##* }" >> slurm/${TRAIT}.run_score.jobs
	fi
done < ${TRAIT}.partitions

DEPENDENCY="afterany"
while read -r line 
do
	DEPENDENCY="${DEPENDENCY}:${line}"
done < slurm/${TRAIT}.run_score.jobs

sbatch --dependency=${DEPENDENCY} --account=apboyle99 --partition=standard --time=0:30:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=1g --job-name=link.run_score.run_optimize_TURF_PRS \
--wrap "bash build_run_optimize_TURF_PRS.sh ${TRAIT} ${ITERATION} ${ROOT_PATH}"