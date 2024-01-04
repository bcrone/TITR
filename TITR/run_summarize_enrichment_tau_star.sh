#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$1.$2.summarize_enrichment_tau_star
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=65g
#SBATCH --time=5:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --output=logs/%x.log

TRAIT="$1"
ITERATION="$2"
RESULTS_PATH="$3"

PREFIX="$1.quantile_normalized.v1_2"

module load python3.9-anaconda

python summarize_enrichment_tau_star.py --path \${RESULTS_PATH} --prefix \${PREFIX} --trait \${TRAIT} --iteration \${ITERATION}

EOF
exit