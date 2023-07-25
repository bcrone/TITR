#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$2.$3.run_standard_scores
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=65g
#SBATCH --time=5:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --output=logs/%x.log

ROOT_PATH="$1"
ANCESTRY="$2"
TRAIT="$3"

EUR_1KG_BFILE="\${ROOT_PATH}/data/1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC"
RANGE_LIST="\${ROOT_PATH}/data/UKB/phenos/range_list.expanded"
SUMSTATS="\${ROOT_PATH}/data/GWAS/\${TRAIT}/\${TRAIT}.PLINK.TITR"

BED_PATH="\${ROOT_PATH}/data/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bed"
BIM_PATH="\${ROOT_PATH}/data/UKB/\${ANCESTRY}/ukb_imp_chrALL_v3.bim"
FAM_PATH="\${ROOT_PATH}/data/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.fam"
SAMPLES_PATH="\${ROOT_PATH}/data/UKB/phenos/\${ANCESTRY}/\${ANCESTRY}.\${TRAIT}.sample.IDs"

RESULTS_PATH="\${ROOT_PATH}/results/single_tissue/\${TRAIT}"

PLINK_PATH="/nfs/turbo/boylelab/plink/plink"
PLINK2_PATH="/nfs/turbo/boylelab/plink2/plink2"

\${PLINK_PATH} --bfile \${EUR_1KG_BFILE} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump \${SUMSTATS} --clump-snp-field SNP --clump-field P \
		--out \${RESULTS_PATH}/standard/\${TRAIT}.\${ANCESTRY}.standard.clump

\${PLINK2_PATH} --bed \${BED_PATH} --bim \${BIM_FILE} --fam \${FAM_FILE} --score \${SUMSTATS} 1 2 3 header no-mean-imputation --q-score-range \${RANGE_LIST} \${SUMSTATS}.SNP.pvalue \
        --extract \${RESULTS_PATH}/standard/\${TRAIT}.\${ANCESTRY}.standard.clump.clumped --keep \${SAMPLES_PATH} --out \${RESULTS_PATH}/standard/\${TRAIT}.\${ANCESTRY}.standard

exit
EOF