#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$2.$1.ITERATION_$3.run_score
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=65g
#SBATCH --time=5:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --output=logs/%x.log

ANCESTRY="$1"
TRAIT="$2"
ITERATION="$3"
ROOT_PATH="$4"

SUMSTATS="\${ROOT_PATH}/data/GWAS/\${TRAIT}/\${TRAIT}.PLINK.TITR"
OUTPATH="\${ROOT_PATH}/results/traits/\${TRAIT}"

PLINK2_PATH="/nfs/turbo/boylelab/plink2/plink2"

BED_FILE="\${ROOT_PATH}/data/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bed"
BIM_FILE="\${ROOT_PATH}/data/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bim"
FAM_FILE="\${ROOT_PATH}/data/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.fam"

\${PLINK2_PATH} --pgen \${BED_FILE} --pvar \${BIM_FILE} --psam \${FAM_FILE} --score \${SUMSTATS} 1 2 3 header no-mean-imputation \
	--extract \${OUTPATH}/\${TRAIT}.TURF.ITERATION_\${ITERATION}.snplist --remove \${ROOT_PATH}/data/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.NA \
	--out \${OUTPATH}/\${TRAIT}.\${ANCESTRY}.ITERATION_\${ITERATION} 

exit
EOF
