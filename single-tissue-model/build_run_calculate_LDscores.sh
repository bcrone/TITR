#!/bin/bash

ROOT_PATH="$1"

rm slurm/run_calculate_LDscores.jobs

TISSUES=()
while IFS= read -r line; do
	TISSUE+=("$line")
done < ../data/TURF_tissues.txt

for TISSUE in "${TISSUES[@]}"
do
	bash run_calculate_LDscores.sh ${ROOT_PATH} ${TISSUE} > slurm/${TISSUE}.run_calculate_LDscores.sh
	RES=$(sbatch slurm/${TISSUE}.run_calculate_LDscores.sh)
	echo "${RES##* }" >> slurm/run_calculate_LDscores.jobs
done