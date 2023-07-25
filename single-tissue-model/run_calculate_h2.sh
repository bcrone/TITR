#!/bin/bash
cat <<EOF
#!/bin/bash
#SBATCH --job-name=$2.$3.run_calculate_h2
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=65g
#SBATCH --time=5:00:00
#SBATCH --account=
#SBATCH --partition=
#SBATCH --output=logs/%x.log

ROOT_PATH="$1"
TRAIT="$2"
TISSUE="$3"

bash calculate_h2.sh \${ROOT_PATH} \${TRAIT} \${TISSUE}

exit
EOF