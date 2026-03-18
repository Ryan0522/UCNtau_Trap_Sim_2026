#!/bin/bash
#SBATCH --job-name=UCN_Sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-20
#SBATCH --output=logs/run_%a.out
#SBATCH --time=24:00:00

# 1. Load modules
module load openmpi
module load python3

# 2. Variables from submit_all.sh
MODE=${1:-Fast}
DEFECT=${2:-0}
OUT_DIR="${MODE}_def_${DEFECT}"

# 3. Determine holdtime index (4 is for 1550s)
HT_IDX=$((SLURM_ARRAY_TASK_ID / 20))

# 4. Run simulation, longer time limit for 1550s
if [ $HT_IDX -eq 4 ]; then
    echo "Running Long-term task (1550s) for Mode: $MODE, Defect: $DEFECT"
else
    echo "Running Short-term task for Mode: $MODE, Defect: $DEFECT"
fi

# 5. execution
cd ..
mpirun python3 python/sim_script.py \
    --mode "$MODE" \
    --defect "$DEFECT" \
    --out_dir "$OUT_DIR"