#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$2.$1.$3.validate_TURF_PRS
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
ROOT="$5"

module load python3.9-anaconda

python validate_TURF_PRS.py --ancestry \${ANCESTRY} --trait \${TRAIT} --iteration \${ITERATION} --isQuant \${ISQUANT} --root \${ROOT}

exit
EOF
