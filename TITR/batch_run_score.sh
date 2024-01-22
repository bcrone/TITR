#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$1.$2.run_score
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=65g
#SBATCH --time=5:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --output=logs/%x.log
#SBATCH --array=$4

ANCESTRY="EUR"
TRAIT="$1"
TISSUE="$2"
ITERATION="$3"
PARTITION="$4"
ROOT_PATH="$5"

EUR_1KG_BFILE="\${ROOT_PATH}/data/1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC"
SUMSTATS="\${ROOT_PATH}/data/GWAS/\${TRAIT}/\${TRAIT}.PLINK.TITR"
CLUMP="\${ROOT_PATH}/data/GWAS/\${TRAIT}/\${TRAIT}.PLINK.dedup"
TISSUE_PATH="\${ROOT_PATH}/data/RegDB/\${TISSUE}/1000G_phase3_master_scores_1KG-pruned_quantile_normalized_chrALL.\${TISSUE}.\${SLURM_ARRAY_TASK_ID}.SNPs"
OUTPATH="\${ROOT_PATH}/results/\${TISSUE}"

PLINK_PATH="/nfs/turbo/boylelab/plink/plink"
PLINK2_PATH="/nfs/turbo/boylelab/plink2/plink2"

BED_FILE="\${ROOT_PATH}/data/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bed"
BIM_FILE="\${ROOT_PATH}/data/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bim"
FAM_FILE="\${ROOT_PATH}/data/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.fam"

RANGE_LIST="\${ROOT_PATH}/data/UKB/phenos/range_list.expanded"
SAMPLE_PATH="\${ROOT_PATH}/data/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.sample.IDs"

module load python3.9-anaconda

python sample_LD_blocks.py --trait \${TRAIT} --tissue \${TISSUE} --iteration \${ITERATION} --partition \${PARTITION}

\${PLINK2_PATH} --bed \${BED_FILE} --bim \${BIM_FILE} --fam \${FAM_FILE} --score \${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range \${RANGE_LIST} \${SUMSTATS}.SNP.pvalue \
        --extract \${OUTPATH}/\${TRAIT}.\${TISSUE}.\${SLURM_ARRAY_TASK_ID}.ITERATION_\${ITERATION}.SNPs --keep \${SAMPLE_PATH} --out \${OUTPATH}/\${TRAIT}.\${TISSUE}.\${SLURM_ARRAY_TASK_ID}.ITERATION_\${ITERATION} \
        --write-snplist 

exit
EOF