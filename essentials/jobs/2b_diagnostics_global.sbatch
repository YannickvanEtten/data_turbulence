#!/bin/bash
# Stage 2b: 21 diagnostics for the 12 global months of 2000 (the calibration year).
#     sbatch --export=ALL,CONVENTION=audit essentials/jobs/2b_diagnostics_global.sbatch
# Output: $BASE/essentials/derived/global_<convention>/  (~200 GB per convention!)
# Measured: 42.7 GB peak, 3-6 h per month (a month once needed more than 8 h).
#SBATCH --job-name=ess-diag-glob
#SBATCH --partition=defq
#SBATCH --array=0-11%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ess-diag-glob-%A_%a.out
set -uo pipefail
source essentials/jobs/_common.sh
run diagnostics --domain global --index "$SLURM_ARRAY_TASK_ID" --convention "${CONVENTION:-audit}"
