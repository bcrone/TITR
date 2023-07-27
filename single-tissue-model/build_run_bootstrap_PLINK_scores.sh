#!/bin/bash

ROOT_PATH="$1"
ANCESTRY="$2"
TRAIT="$3"
ISQUANT="$4"

rm slurm/${TRAIT}.${ANCESTRY}.run_bootstrap_PLINK_scores.jobs

bash run_bootstrap_PLINK_scores.sh ${ROOT_PATH} ${ANCESTRY} ${TRAIT} ${ISQUANT} > slurm/${TRAIT}.${ANCESTRY}.run_bootstrap_PLINK_scores.sh

RES=$(sbatch slurm/${TRAIT}.${ANCESTRY}.run_bootstrap_PLINK_scores.sh)
echo "${RES##* }" >> slurm/${TRAIT}.${ANCESTRY}.run_bootstrap_PLINK_scores.jobs