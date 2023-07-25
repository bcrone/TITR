#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$2.run_calculate_LDscores
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=65g
#SBATCH --time=5:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --output=logs/%x.log
#SBATCH --array=1:22

ROOT_PATH="$1"
TISSUE="$2"

bash calculate_LDscores.sh \${ROOT_PATH} \${TISSUE} \${SLURM_ARRAY_TASK_ID}

exit
EOF