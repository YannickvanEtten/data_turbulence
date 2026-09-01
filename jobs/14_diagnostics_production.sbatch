#!/bin/bash
# ===========================================================================
# 14_diagnostics_production.sbatch -- the 21 diagnostics over ALL 504 months,
# 1979-01 to 2020-12. Phase 4's derived half.
#
#     mkdir -p logs && sbatch jobs/14_diagnostics_production.sbatch
#
# THIS IS DESIGNED TO RUN WHILE THE DOWNLOAD IS STILL GOING
# --------------------------------------------------------
# The raw fetch is CDS-bound and takes days; this is CPU-bound and takes hours.
# There is no reason to serialise them. Each task handles three states:
#
#   raw GRIB not downloaded yet   -> prints NOT-YET-DOWNLOADED, exits 0
#   output already current        -> skips in about a second
#   anything else                 -> computes
#
# So you submit it now, it works through whatever has landed, and you submit
# the SAME array again when the download finishes to pick up the rest. That is
# the resumable pattern STATUS.md 11.9 measured on the download side: recovery
# is "resubmit the same array", with no hand-maintained list of what is left.
#
# WHY "output exists" IS NOT ENOUGH, AND WHY jobs/07's CHECK WOULD BE WRONG HERE
# -----------------------------------------------------------------------------
# jobs/07 skips whenever the output directory exists. That was right when every
# output had been built the same way. It is WRONG now, because two conventions
# changed:
#
#   before 2026-08-29  `deformation` holds DEF^2, not DEF   (FORMULA_AUDIT 5)
#   before 2026-08-30  `f2d` holds variant A, not C         (FORMULA_AUDIT 10.4)
#
# The 27 DJF months already in derived/north_atlantic/ are BOTH. A directory
# check would preserve them and produce a series that is 27 months of one
# convention and 477 of another -- a real, plausible-looking dataset quietly
# built from two different definitions. `--skip-if-matching` reads the output's
# recorded provenance instead and recomputes anything stale.
#
# ORDER OF OPERATIONS
# -------------------
#   1. sbatch jobs/12_tests.sbatch                 (the gate -- must be green)
#   2. sbatch --array=0-11%4 jobs/14_...sbatch     (12-month canary)
#      then check one output:
#      pixi run python ada/check_composition.py \
#          $BASE/derived/north_atlantic/diagnostics_na_1979-01.zarr
#   3. sbatch jobs/14_diagnostics_production.sbatch    (the full 504)
#   4. resubmit (3) once the download finishes
#
# Step 2 costs half an hour and catches a configuration mistake before 500
# months of it.
#
# CONCURRENCY -- READ THIS BEFORE SUBMITTING
# ------------------------------------------
# The QOS caps you at 8 RUNNING JOBS (STATUS.md 11.9) and the download array
# holds 4. This asks for 4, which fills the cap exactly and leaves nothing for
# jobs/12 or the audit checks. That is deliberate -- it is the most work you
# can have in flight -- but if you need a slot back, lower it:
#
#     scontrol update JobId=<thisjob> ArrayTaskThrottle=3
#
# --mem=24G   The same request jobs/07 used for these exact months, measured at
#             5-8 min each (STATUS.md 13.3), and the request 11.9 measured
#             starting in 13 seconds where 120 GB waited 22 hours. Do not raise
#             it. Expect some stretching from shared-filesystem contention with
#             the download: 13.2 measured six concurrent tasks turning a
#             30-minute job into 36, mostly in GRIB decode.
# --time=02:00  Generous against a 5-8 min job. Short honest walltimes backfill.
# ===========================================================================

#SBATCH --job-name=cat-diag-prod
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-503%4
#SBATCH --output=logs/diagprod-%A_%a.out
#SBATCH --error=logs/diagprod-%A_%a.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

# Month index -> YYYY-MM. Identical arithmetic to 1_download_hpc.py's --index,
# so task N here is the same month as task N there:
#     index = (year - 1979) * 12 + (month - 1)
# Confirmed against STATUS.md 11.9's worked example, where 1985-01 is 72 and
# 1995-01 is 192.
IDX=${SLURM_ARRAY_TASK_ID}
YEAR=$(( 1979 + IDX / 12 ))
MON=$(( IDX % 12 + 1 ))
MONTH=$(printf "%04d-%02d" "$YEAR" "$MON")

IN="$BASE/raw/north_atlantic/era5_na_${MONTH}.grib"
OUT="$BASE/derived/north_atlantic/diagnostics_na_${MONTH}.zarr"

echo "=== [$IDX] $MONTH on $(hostname) at $(date -Is) ==="
echo "    repo $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
echo "    in   $IN"
echo "    out  $OUT"

if [ ! -f "$IN" ]; then
    # Not an error while the download is still running. Printed as a single
    # greppable token so a whole array's worth can be counted at a glance:
    #     grep -l NOT-YET-DOWNLOADED logs/diagprod-*.out | wc -l
    echo "NOT-YET-DOWNLOADED $MONTH — resubmit this array after jobs/01 finishes"
    exit 0
fi

mkdir -p "$(dirname "$OUT")"

start=$(date +%s)
# --chunk-days 0 matches jobs/07 exactly, so months computed now are identical
# to the ones computed for the trend check. Chunking is bit-identical anyway
# (STATUS.md 11.11), but keeping the flag the same removes the question.
pixi run python ada/diagnostics_global.py "$IN" "$OUT" \
    --chunk-days 0 --skip-if-matching
status=$?
end=$(date +%s)

echo
echo "=== [$IDX] $MONTH exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
exit $status
