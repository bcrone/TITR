#!/bin/bash
ANCESTRY="$1"
TRAIT="$2"
ITERATION="$3"
ISQUANT=""
ROOT=""

rm slurm/${TRAIT}.run_validate_TURF.jobs

bash run_validate_TURF_PRS.sh ${ANCESTRY} ${TRAIT} ${ITERATION} ${ISQUANT} ${ROOT} > slurm/${TRAIT}.ITERATION_${ITERATION}.validate_TURF_PRS.sh
RES=$(sbatch slurm/${TRAIT}.ITERATION_${ITERATION}.validate_TURF_PRS.sh)
echo "${RES##* }" >> slurm/${TRAIT}.run_validate_TURF.jobs