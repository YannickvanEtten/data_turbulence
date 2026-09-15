#!/bin/bash
# Second check: Stages 3 and 4 of essentials on the ORIGINAL stores.
#     sbatch --export=ALL,CONVENTION=audit    essentials/jobs/0b_check_thresholds_and_trends.sbatch
#     sbatch --export=ALL,CONVENTION=baseline essentials/jobs/0b_check_thresholds_and_trends.sbatch
# Expected (audit): thresholds equal to thresholds_2026-09-11.json (relative difference ~0),
# and the trend table of STATUS 18.15: 25/25 significant, 25/25 Prosser inside the CI,
# annual trend ratios 0.97 / 0.98 / 1.01 / 1.02 / 1.04.
# (baseline: thresholds_2026-09-07.json and RESULT_full_trend_42yr.md.)
#SBATCH --job-name=ess-check34
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ess-check34-%j.out
set -uo pipefail
source essentials/jobs/_common.sh
CONV=${CONVENTION:-audit}
if [ "$CONV" = audit ]; then GLOB=global_audit; NA=north_atlantic_audit; OLD=thresholds_2026-09-11.json
else GLOB=global; NA=north_atlantic; OLD=thresholds_2026-09-07.json; fi
run thresholds --convention "$CONV" --derived "$BASE/derived/$GLOB" \
    --out "$BASE/essentials/calibration/check_thresholds_${CONV}.json" \
    --compare-with "$BASE/calibration/$OLD"
run trends --convention "$CONV" --derived "$BASE/derived/$NA" --thresholds "$BASE/calibration/$OLD"
