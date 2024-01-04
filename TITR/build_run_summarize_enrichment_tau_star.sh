#!/bin/bash
TRAIT="$1"
ITERATION="$2"
ROOT_PATH="$3"

rm slurm/${TRAIT}.summarize_enrichment_tau_star.jobs

RESULTS_PATH="${ROOT_PATH}/results"

bash run_summarize_enrichment_tau_star.sh ${TRAIT} ${ITERATION} ${RESULTS_PATH} > slurm/${TRAIT}.ITERATION_${ITERATION}.summarize_enrichment_tau_star.sh
RES=$(sbatch slurm/${TRAIT}.ITERATION_${ITERATION}.summarize_enrichment_tau_star.sh)
echo "${RES##* }" >> slurm/${TRAIT}.summarize_enrichment_tau_star.jobs

DEPENDENCY="afterany"
while read -r line 
do
	DEPENDENCY="${DEPENDENCY}:${line}"
done < slurm/${TRAIT}.summarize_enrichment_tau_star.jobs

sbatch --dependency=${DEPENDENCY} --account=apboyle99 --partition=standard --time=0:30:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=1g --job-name=link.run_summarzie_enrichment_tau_star \
--wrap "bash build_run_score.sh ${TRAIT} ${ITERATION} ${ROOT_PATH}"