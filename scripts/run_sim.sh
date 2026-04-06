#!/bin/bash
#SBATCH --job-name=UCN_Sim
#SBATCH --account=chenyliu
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --array=0-19
#SBATCH --mem=16G
#SBATCH --output=scripts/logs/run_%a.out
#SBATCH --time=00:10:00

# 1. Load modules
module load gcc
module load openmpi
module load python

# 2. Variables from submit_all.sh
MODE=${1:-Segmented}
DEFECT=${2:-2e-5}
OUT_DIR=${3:-"${MODE}_def_${DEFECT}"}

# 3. Determine holdtime index (4 is for 1550s)
HT_IDX=$((SLURM_ARRAY_TASK_ID / 4))

# 4. Run simulation, longer time limit for 1550s
if [ $HT_IDX -eq 4 ]; then
    echo "Running Long-term task (1550s) for Mode: $MODE, Defect: $DEFECT"
else
    echo "Running Short-term task for Mode: $MODE, Defect: $DEFECT"
fi

# 5. execution
cd ./python
mpirun --bind-to none -- python3 sim_script.py  \
    --mode "$MODE" \
    --defect "$DEFECT" \
    --out_dir "$OUT_DIR" \
    --ntraj 1000 \
    --tracker "energy"