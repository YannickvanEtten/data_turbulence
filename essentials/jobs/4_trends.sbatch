#!/bin/bash
# Stage 4: exceedance frequencies and 1979-2020 trends, all seasons, all 21 + ensemble.
#     sbatch --export=ALL,CONVENTION=audit essentials/jobs/4_trends.sbatch
# Optional: DERIVED=<folder with the 504 NA stores>, THRESHOLDS=<json>.
# Output: $BASE/essentials/results/{frequencies,trends}_<convention>_1979-2020.csv
# Memory: one diagnostic-year of the box at a time (~0.2 GB) plus counters.
#SBATCH --job-name=ess-trends
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/ess-trends-%j.out
set -uo pipefail
source essentials/jobs/_common.sh
args=(trends --convention "${CONVENTION:-audit}")
[ -n "${DERIVED:-}" ] && args+=(--derived "$DERIVED")
[ -n "${THRESHOLDS:-}" ] && args+=(--thresholds "$THRESHOLDS")
run "${args[@]}"
