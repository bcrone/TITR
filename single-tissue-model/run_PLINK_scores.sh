#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$2.$3.run_PLINK_scores
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=65g
#SBATCH --time=5:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --output=logs/%x.log
#SBATCH --array=10,50,100,200,500

ROOT_PATH="$1"
ANCESTRY="$2"
TRAIT="$3"
TISSUE="$4"

EUR_1KG_BFILE="\${ROOT_PATH}/data/1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC"
RANGE_LIST="\${ROOT_PATH}/data/UKB/phenos/range_list.expanded"
SUMSTATS="\${ROOT_PATH}/data/GWAS/\${TRAIT}/\${TRAIT}.PLINK.TITR"

BED_PATH="\${ROOT_PATH}/data/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bed"
BIM_PATH="\${ROOT_PATH}/data/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bim"
FAM_PATH="\${ROOT_PATH}/data/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.fam"
SAMPLES_PATH="\${ROOT_PATH}/data/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.sample.IDs"

IMPACT_PATH="\${ROOT_PATH}/data/IMPACT/\${TRAIT}"
SURF_PATH=\${ROOT_PATH}/data/RegDB/GENERIC
TURF_PATH="\${ROOT_PATH}/data/RegDB/\${TISSUE}"

RESULTS_PATH="\${ROOT_PATH}/results/single_tissue/\${TRAIT}"

PLINK_PATH="/nfs/turbo/boylelab/plink/plink"
PLINK2_PATH="/nfs/turbo/boylelab/plink2/plink2"

#IMPACT
\${PLINK_PATH} --bfile \${EUR_1KG_BFILE} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump \${SUMSTATS} --clump-snp-field SNP --clump-field P \
        --extract \${IMPACT_PATH}/IMPACT707_EUR_chrALL.\${SLURM_ARRAY_TASK_ID}.\${TRAIT}.SNPs --out \${RESULTS_PATH}/IMPACT/\${TRAIT}.\${ANCESTRY}.IMPACT.\${SLURM_ARRAY_TASK_ID}.clump

\${PLINK2_PATH} --pgen \${BED_PATH} --pvar \${BIM_PATH} --psam \${FAM_PATH} --score \${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range \${RANGE_LIST} \${SUMSTATS}.SNP.pvalue \
        --extract \${RESULTS_PATH}/IMPACT/\${TRAIT}.\${ANCESTRY}.IMPACT.\${SLURM_ARRAY_TASK_ID}.clump.clumped --keep \${SAMPLES_PATH} --out \${RESULTS_PATH}/IMPACT/\${TRAIT}.\${ANCESTRY}.IMPACT.\${SLURM_ARRAY_TASK_ID}

#SURF
\${PLINK_PATH} --bfile \${EUR_1KG_BFILE} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump \${SUMSTATS} --clump-snp-field SNP --clump-field P \
        --extract \${SURF_PATH}/1000G_phase3_master_scores_1KG-pruned_quantile_normalized_chrALL.\${SLURM_ARRAY_TASK_ID}.GENERIC.SNPs --out \${RESULTS_PATH}/SURF/\${TRAIT}.\${ANCESTRY}.SURF.\${SLURM_ARRAY_TASK_ID}.clump

\${PLINK2_PATH} --pgen \${BED_PATH} --pvar \${BIM_PATH} --psam \${FAM_PATH} --score \${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range \${RANGE_LIST} \${SUMSTATS}.SNP.pvalue \
        --extract \${RESULTS_PATH}/SURF/\${TRAIT}.\${ANCESTRY}.SURF.\${SLURM_ARRAY_TASK_ID}.clump.clumped --keep \${SAMPLES_PATH} --out \${RESULTS_PATH}/SURF/\${TRAIT}.\${ANCESTRY}.SURF.\${SLURM_ARRAY_TASK_ID}

#TURF
\${PLINK_PATH} --bfile \${EUR_1KG_BFILE} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump \${SUMSTATS} --clump-snp-field SNP --clump-field P \
        --extract \${TURF_PATH}/1000G_phase3_master_scores_1KG-pruned_quantile_normalized_chrALL.\${SLURM_ARRAY_TASK_ID}.\${TISSUE}.SNPs --out \${RESULTS_PATH}/TURF/\${TRAIT}.\${ANCESTRY}.TURF.\${SLURM_ARRAY_TASK_ID}.clump

\${PLINK2_PATH} --pgen \${BED_PATH} --pvar \${BIM_PATH} --psam \${FAM_PATH} --score \${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range \${RANGE_LIST} \${SUMSTATS}.SNP.pvalue \
        --extract \${RESULTS_PATH}/TURF/\${TRAIT}.\${ANCESTRY}.TURF.\${SLURM_ARRAY_TASK_ID}.clump.clumped --keep \${SAMPLES_PATH} --out \${RESULTS_PATH}/TURF/\${TRAIT}.\${ANCESTRY}.TURF.\${SLURM_ARRAY_TASK_ID}

exit
EOF