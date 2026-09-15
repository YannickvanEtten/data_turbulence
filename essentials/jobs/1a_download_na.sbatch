#!/bin/bash
# Stage 1a: the 504 North Atlantic months (index 0 = 1979-01 ... 503 = 2020-12).
# %4 is the CDS limit (8 concurrent requests failed 8 of 12). Submit BARE:
#     mkdir -p logs && sbatch essentials/jobs/1a_download_na.sbatch
# Passing --array on the command line silently drops the %4.
# Months already on disk are verified and skipped in seconds, so re-submitting is free.
#SBATCH --job-name=ess-dl-na
#SBATCH --partition=defq
#SBATCH --array=0-503%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ess-dl-na-%A_%a.out
set -uo pipefail
source essentials/jobs/_common.sh
run download-na --index "$SLURM_ARRAY_TASK_ID"
