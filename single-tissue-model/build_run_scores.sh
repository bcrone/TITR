#!/bin/bash

ROOT_PATH="$1"
ANCESTRY="$2"
TRAIT="$3"
TISSUE="$4"

rm slurm/run_calculate_PLINK_scores.jobs
rm slurm/run_calculate_standard_scores.jobs

bash run_standard_scores.sh ${ROOT_PATH} ${ANCESTRY} ${TRAIT} > slurm/${TRAIT}.${ANCESTRY}.run_standard_scores.sh
bash run_PLINK_scores.sh ${ROOT_PATH} ${ANCESTRY} ${TRAIT} ${TISSUE} > slurm/${TRAIT}.${ANCESTRY}.run_PLINK_scores.sh

RES=$(sbatch slurm/${TRAIT}.${ANCESTRY}.run_standard_scores.sh)
echo "${RES##* }" >> slurm/run_standard_scores.jobs
RES=$(sbatch slurm/${TRAIT}.${ANCESTRY}.run_PLINK_scores.sh)
echo "${RES##* }" >> slurm/run_PLINK_scores.jobs