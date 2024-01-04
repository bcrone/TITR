#!/bin/bash
ANCESTRY="$1"
TRAIT="$2"
ITER=$3
ROOT_PATH="$4"

rm slurm/${TRAIT}.run_validate_score.jobs

for ITERATION in $( seq 1 $ITER )
do
	bash batch_validate_score.sh ${ANCESTRY} ${TRAIT} ${ITERATION} ${ROOT_PATH} > slurm/${TRAIT}.${ITERATION}.run_validate.sh
	RES=$(sbatch slurm/${TRAIT}.${ITERATION}.run_validate.sh)
	echo "${RES##* }" >> slurm/${TRAIT}.run_validate_score.jobs
done
