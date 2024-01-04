#!/bin/bash
ANCESTRY="EUR"
TRAIT="$1"
ITERATION="$2"
ISQUANT=""
ROOT_PATH="$3"

rm slurm/${TRAIT}.optimize_TURF_PRS.jobs

bash run_optimize_TURF_PRS.sh ${ANCESTRY} ${TRAIT} ${ITERATION} ${ISQUANT} ${ROOT_PATH} > slurm/${TRAIT}.ITERATION_${ITERATION}.optimize_TURF_PRS.sh
RES=$(sbatch slurm/${TRAIT}.ITERATION_${ITERATION}.optimize_TURF_PRS.sh)
echo "${RES##* }" >> slurm/${TRAIT}.optimize_TURF_PRS.jobs

DEPENDENCY="afterany"
while read -r line 
do
	DEPENDENCY="${DEPENDENCY}:${line}"
done < slurm/${TRAIT}.optimize_TURF_PRS.jobs

sbatch --dependency=${DEPENDENCY} --account=apboyle99 --partition=standard --time=0:30:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=1g --job-name=link.run_optimize_TURF_PRS.clean_iteration.sh \
--wrap "bash clean_iteration.sh ${TRAIT} $((ITERATION)) ${ROOT_PATH}"