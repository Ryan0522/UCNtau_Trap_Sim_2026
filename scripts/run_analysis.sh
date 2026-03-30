#!/bin/bash
#SBATCH --job-name=UCN_Analysis
#SBATCH --account=chenyliu
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=scripts/logs/analysis_%j.out

module load python
export PYTHONPATH=/projects/illinois/eng/physics/chenyliu/Ryan_ciyouh2/mypython3:${PYTHONPATH}

OUT_DIR=${1:-test}

cd ./python
python3 analysis.py --out_dir "$OUT_DIR"