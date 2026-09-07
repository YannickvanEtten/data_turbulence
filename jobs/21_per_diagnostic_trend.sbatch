#!/bin/bash
# ---------------------------------------------------------------------------
# 21_per_diagnostic_trend.sbatch
#
# The fitted per-diagnostic trend over all 42 years, with intervals.
# Re-examines STATUS.md §12.6 -- which is a NINE-SEASON, DJF-ONLY, 1979-vs-2020
# ENDPOINT table, and the sole basis for §7's top open item ("ubf is last on
# trend at +11% where every sibling is +38% to +226%").
#
# WHY THIS IS WORTH A JOB SLOT
#     §12.5 was built on the same nine seasons and the full record overturned
#     it outright: the "trend excess" it described (ratios 1.03 -> 1.43 rising
#     with severity) is absent at n=42, where the DJF ratios run 0.84 -> 0.98.
#     §12.4 had already explained why -- the nine-season fit leaned on two
#     strongly positive-NAO winters at the end of the sample. §12.6 shares that
#     sample and that construction, so it is under the same suspicion, and
#     `ubf` may not be anomalous at all.
#
#     NO NEW DATA IS REQUIRED. The 504 North Atlantic stores already exist and
#     are ALL CLEAR (§14.7). This is the cheapest open question in the project.
#
# WHY NOT full_trend_check.py --per-diagnostic
#     That flag reports the first and last year and the percentage between
#     them. Two points, no interval. §12.8 retired exactly that construction:
#     "a comparison of point estimates from a small sample is not a test ...
#     should carry an interval, or it is decoration." This fits all 42 years
#     and reports t, R^2 and a 95% interval per diagnostic.
#
# THRESHOLDS
#     Defaults to the newest calibration/thresholds_*.json, which after
#     2026-09-07 is the FULL reference year. Pass --thresholds explicitly to
#     pin a different one; the script prints which it used, as jobs/16 does.
#
# SIZING
#     Loads one year of North Atlantic diagnostics at a time and subsets to
#     Prosser's box, the same access pattern as jobs/16 (which ran 4 h 06 m at
#     48G for five seasons x five severities). This does one season and one
#     severity but 21 separate exceedance fields, so the loading dominates and
#     the total should be lower. 64G rather than 48G because memory limits are
#     NOT enforced on this cluster -- an under-request fails silently and can
#     take out a neighbour's job instead of your own.
#
#   sbatch jobs/21_per_diagnostic_trend.sbatch                 # annual, MOG
#   sbatch --export=ALL,SEASON=djf jobs/21_per_diagnostic_trend.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-perdiag
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/perdiag-%j.out
#SBATCH --error=logs/perdiag-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
SEASON=${SEASON:-annual}
SEVERITY=${SEVERITY:-moderate}

module load 2025
module load pixi

export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

mkdir -p cat_outputs

echo "=== per-diagnostic trend, season=$SEASON severity=$SEVERITY on $(hostname) at $(date -Is) ==="

start=$(date +%s)
pixi run python ada/per_diagnostic_trend.py \
    --base "$BASE" \
    --season "$SEASON" \
    --severity "$SEVERITY" \
    --csv "cat_outputs/per_diagnostic_${SEASON}_${SEVERITY}.csv"
status=$?
end=$(date +%s)

echo "=== exit=$status in $(( (end - start) / 60 )) min at $(date -Is) ==="
exit $status
