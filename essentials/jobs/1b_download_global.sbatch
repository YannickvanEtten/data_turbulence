#!/bin/bash
# Stage 1b: the global year 2000 in 36 day-blocks (index = 3*(month-1) + block).
# Then run 1c to join the blocks into 12 months. Submit bare (see 1a).
#SBATCH --job-name=ess-dl-glob
#SBATCH --partition=defq
#SBATCH --array=0-35%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ess-dl-glob-%A_%a.out
set -uo pipefail
source essentials/jobs/_common.sh
run download-global --index "$SLURM_ARRAY_TASK_ID"
