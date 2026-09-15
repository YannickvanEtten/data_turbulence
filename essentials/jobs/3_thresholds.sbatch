#!/bin/bash
# Stage 3: severity thresholds from the 12 global stores.
#     sbatch --export=ALL,CONVENTION=audit essentials/jobs/3_thresholds.sbatch
# Optional: DERIVED=<folder with the 12 stores>, ARCHIVE=1 (write the ~41 GB tail archive).
# Measured: 55.7 GB peak, ~2 h.
#SBATCH --job-name=ess-thresholds
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ess-thresholds-%j.out
set -uo pipefail
source essentials/jobs/_common.sh
args=(thresholds --convention "${CONVENTION:-audit}")
[ -n "${DERIVED:-}" ] && args+=(--derived "$DERIVED")
[ -n "${ARCHIVE:-}" ] && args+=(--archive)
run "${args[@]}"
