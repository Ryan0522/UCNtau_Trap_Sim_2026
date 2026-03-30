#!/bin/bash

MODE="Segmented"
DEFECT="2e-5"
OUT_DIR="${MODE}_def_${DEFECT}"

echo "--------------------------------------------------"
echo " Launching UCN Production Pipeline"
echo " Mode:   $MODE"
echo " Defect: $DEFECT"
echo " Folder: $OUT_DIR"
echo "--------------------------------------------------"

# 1. Submit simulation job
SIM_ID=$(sbatch --parsable scripts/run_sim.sh "$MODE" "$DEFECT")

if [ -z "$SIM_ID" ]; then
    echo "Error: Failed to submit simulation job."
    exit 1
fi

echo "Step 1: Simulation Array Job submitted. ID: $SIM_ID"

# 2. Analysis job submission (with afterok:$SIM_ID)
ANA_ID=$(sbatch --parsable --dependency=afterok:$SIM_ID scripts/run_analysis.sh "$OUT_DIR")

echo "Step 2: Analysis Job queued. ID: $ANA_ID"
echo "--------------------------------------------------"
echo "Pipeline established. Results will be in results/$OUT_DIR"