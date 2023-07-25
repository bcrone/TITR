#!/bin/bash

ROOT_PATH="$1"
TRAIT="$2"

rm slurm/run_calculate_h2.jobs

TISSUES=()
while IFS= read -r line; do
	TISSUE+=("$line")
done < ../data/TURF_tissues.txt

for TISSUE in "${TISSUES[@]}"
do
	bash run_calculate_h2.sh ${ROOT_PATH} ${TRAIT} ${TISSUE} > slurm/${TRAIT}.${TISSUE}.run_calculate_LDscores.sh
	RES=$(sbatch slurm/${TRAIT}.${TISSUE}.run_calculate_LDscores.sh)
	echo "${RES##* }" >> slurm/run_calculate_LDscores.jobs
done