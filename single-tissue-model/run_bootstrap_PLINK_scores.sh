#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$2.$3.run_bootstrap_PLINK_scores
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
ISQUANT="$4"

module load python3.9-anaconda

python bootstrap_PLINK_scores.py --ancestry \${ANCESTRY} --trait \${TRAIT} --isQuant \${ISQUANT} --root \${ROOT_PATH}

exit
EOF