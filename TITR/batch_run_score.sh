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

EUR_1KG_BFILE="/path/to/1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC"
SUMSTATS="/path/to/GWAS/\${TRAIT}/\${TRAIT}.PLINK.TITR"
CLUMP="/path/to/GWAS/\${TRAIT}/\${TRAIT}.PLINK.dedup"
TISSUE_PATH="/path/to/RegDB/\${TISSUE}/1000G_phase3_master_scores_1KG-pruned_quantile_normalized_chrALL.\${TISSUE}.\${SLURM_ARRAY_TASK_ID}.SNPs"
OUTPATH="/path/to/results/\${TISSUE}"

PLINK_PATH="/path/to/plink"
PLINK2_PATH="/path/to/plink2"

BED_FILE="/path/to/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bed"
BIM_FILE="/path/to/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bim"
FAM_FILE="/path/to/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.fam"

RANGE_LIST="/path/to/UKB/phenos/range_list.expanded"
SAMPLE_PATH="/path/to/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.sample.IDs"

module load python3.9-anaconda

python sample_LD_blocks.py --trait \${TRAIT} --tissue \${TISSUE} --iteration \${ITERATION} --partition \${PARTITION}

\${PLINK2_PATH} --bed \${BED_FILE} --bim \${BIM_FILE} --fam \${FAM_FILE} --score \${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range \${RANGE_LIST} \${SUMSTATS}.SNP.pvalue \
        --extract \${OUTPATH}/\${TRAIT}.\${TISSUE}.\${SLURM_ARRAY_TASK_ID}.ITERATION_\${ITERATION}.SNPs --keep \${SAMPLE_PATH} --out \${OUTPATH}/\${TRAIT}.\${TISSUE}.\${SLURM_ARRAY_TASK_ID}.ITERATION_\${ITERATION} \
        --write-snplist 

exit
EOF