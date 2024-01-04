#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$2.$3.optimize_TURF_PRS
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
ISQUANT="$4"
ROOT_PATH="$5"

module load python3.9-anaconda

python optimize_TURF_PRS.py --ancestry \${ANCESTRY} --trait \${TRAIT} --iteration \${ITERATION} --isQuant \${ISQUANT} --root \${ROOT_PATH}

exit
EOF