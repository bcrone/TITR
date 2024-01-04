#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$1.$2.run_ldsc
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=65g
#SBATCH --time=5:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --output=logs/%x.log

ANCESTRY="EUR"
TRAIT="$1"
TISSUE="$2"
ITERATION="$3"
ROOT_PATH="$4"

SUMSTATS="/path/to/GWAS/sumstats/\${TRAIT}.ITERATION_\${ITERATION}.sumstats.gz"
REF_LD_CHR="/path/to/annotation/ldscore/customized_baselineLD_cts.RegDB.quantile_normalized.\${TISSUE}."
FRQFILE_CHR="/path/to/1KG/1000G_Phase3_frq/1000G.EUR.QC."
W_LD_CHR="/path/to/1KG/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
OUTPATH="/path/to/results/\${TISSUE}/\${TRAIT}.quantile_normalized.v1_2.\${TISSUE}"

module load python3.9-anaconda
source activate ldsc

/nfs/turbo/boylelab/ldsc/ldsc.py --h2 \${SUMSTATS} --ref-ld-chr \${REF_LD_CHR} --out \${OUTPATH} --overlap-annot --frqfile-chr \${FRQFILE_CHR} --w-ld-chr \${W_LD_CHR} --print-coefficients --print-delete-vals

exit
EOF