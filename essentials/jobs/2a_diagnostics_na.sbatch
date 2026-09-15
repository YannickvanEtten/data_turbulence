#!/bin/bash
# Stage 2a: 21 diagnostics for the 504 North Atlantic months.
#     sbatch --export=ALL,CONVENTION=audit    essentials/jobs/2a_diagnostics_na.sbatch
#     sbatch --export=ALL,CONVENTION=baseline essentials/jobs/2a_diagnostics_na.sbatch
# Output: $BASE/essentials/derived/north_atlantic_<convention>/  (~300 GB per convention!)
# Measured: ~15 GB peak, 5-9 min per month; %8 is the QOS job cap; ~7 h for 504.
# Stores already written under the same convention are skipped.
#SBATCH --job-name=ess-diag-na
#SBATCH --partition=defq
#SBATCH --array=0-503%8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ess-diag-na-%A_%a.out
set -uo pipefail
source essentials/jobs/_common.sh
run diagnostics --domain north_atlantic --index "$SLURM_ARRAY_TASK_ID" --convention "${CONVENTION:-audit}"
