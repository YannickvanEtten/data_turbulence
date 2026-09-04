#!/bin/bash
# ===========================================================================
# 16_full_trend_check.sbatch -- the DEFINITIVE trend check: all 42 years,
# all four seasons plus annual, against the FULL Prosser (2023) Table 1.
#
#     sbatch jobs/16_full_trend_check.sbatch
#     sbatch jobs/16_full_trend_check.sbatch --season djf --per-diagnostic
#
# REQUIRES jobs/15_check_production_complete.sbatch to have reported ALL
# CLEAR first. This script does not re-verify provenance itself.
#
# --mem=48G   ada/trend_check.py's DJF-only (3 months) run used 16G. This
#             script loads a full calendar year (12 months) of the North
#             Atlantic box at once before slicing into seasons, so budget
#             roughly 4x -- and one year's data is discarded before the next
#             loads, so this does not grow with the number of years. If it
#             still runs out of memory on a busy node, split by season:
#             `--season djf`, `--season mam`, etc. across separate jobs.
# --time=04:00:00   42 years x up to 5 seasons of box-subset reduction. Sized
#             generously first time; tighten once a real wall time is measured
#             (print it in the log, per the project's own convention).
# ===========================================================================

#SBATCH --job-name=cat-full-trend
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
# 04:00:00 was NOT enough. Job 1098191 was still running at 3:20 on the first
# real 504-month attempt -- 42 calendar years opened one at a time, sliced into
# five seasons, against 9 DJF seasons in ada/trend_check.py. Raised to 12 h on
# 2026-09-04. Short honest walltimes get backfilled (STATUS 11.9), so this is
# not free, but a job killed at 3:59 with nothing written costs far more.
# Tighten once a real elapsed time has actually been measured and logged.
#SBATCH --time=12:00:00
#SBATCH --output=logs/full-trend-%j.out
#SBATCH --error=logs/full-trend-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== full 42-year trend check on $(hostname) at $(date -Is) ==="

start=$(date +%s)
pixi run python ada/full_trend_check.py --base "$BASE" "$@"
status=$?
end=$(date +%s)

echo
echo "=== exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
exit $status
