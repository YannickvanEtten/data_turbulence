#!/bin/bash
# ===========================================================================
# 07_diagnostics_na_trend.sbatch -- the 21 diagnostics on the North Atlantic
# months of NINE DJF seasons: 1979, 1985, 1990, 1995, 2000, 2005, 2010, 2015,
# 2020.
#
#     sbatch jobs/07_diagnostics_na_trend.sbatch
#
# 27 array tasks. The fifteen months from the five-season run already have
# output and those tasks detect it and skip, so this can be submitted whole
# rather than with a hand-maintained index list. Re-running after a partial
# failure is therefore free and safe. Expect ~12 tasks of real work.
#
# WHY NINE AND NOT FIVE
# ---------------------
# The five-season run fitted a slope with n=5, 3 degrees of freedom. On MOG it
# gave R^2 = 0.59, which is t = sqrt(0.59 * 3 / 0.41) = 2.1, p ~ 0.13 -- not
# significant at 5 %, and a 95 % interval on the +47 % trend running roughly
# -25 % to +165 %. Prosser's +37 % sat inside that interval, so the comparison
# could not have failed and therefore could not have passed either.
#
# At n=9 with the same R^2, t = sqrt(0.59 * 7 / 0.41) = 3.2, p ~ 0.016. That
# is the difference between a number and a result. Twelve extra months of a
# North Atlantic box is a cheap way to buy it.
#
# Prerequisite: the twelve new months downloaded via
#     sbatch --array=72,73,83,192,193,203,312,313,323,432,433,443%4 \
#            jobs/01_download.sbatch
# (index = (year-1979)*12 + month-1 over 01_download.sbatch's 1979-01..2020-12
# span; no new download script needed.)
#
# THE %4 IS NOT OPTIONAL. Giving --array on the command line replaces the
# script's whole array spec including its throttle, and eight simultaneous CDS
# requests return 400 Bad Request. Eight of twelve downloads were lost that
# way once already.
#
# --mem=24G      MEASURED: NA months peaked at 12.9-14.3 GB whole-file.
#                Same request as the global job, and it fits every defq node.
# --time=01:00   MEASURED: 5-8 min each. Short walltimes get backfilled.
# --chunk-days 0 A North Atlantic month fits comfortably in one pass, so F2D
#                loses only the 2 timesteps at each month's true ends.
#
# NOTE ON MONTH BOUNDARIES: each DJF is assembled from three separately
# downloaded months, so F2D loses 6 of ~2160 timesteps per season, 0.3 %.
# Identical treatment in every season, so it cannot bias the TREND -- which is
# the only thing this run is for. Not worth wiring chunk_stitch.py in for.
# ===========================================================================

#SBATCH --job-name=cat-diag-trend
#SBATCH --partition=defq
#SBATCH --array=0-26
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/diag-trend-%A_%a.out
#SBATCH --error=logs/diag-trend-%A_%a.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

MONTHS=(
  1979-01 1979-02 1979-12
  1985-01 1985-02 1985-12
  1990-01 1990-02 1990-12
  1995-01 1995-02 1995-12
  2000-01 2000-02 2000-12
  2005-01 2005-02 2005-12
  2010-01 2010-02 2010-12
  2015-01 2015-02 2015-12
  2020-01 2020-02 2020-12
)
MONTH=${MONTHS[$SLURM_ARRAY_TASK_ID]}

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

IN="$BASE/raw/north_atlantic/era5_na_${MONTH}.grib"
OUT="$BASE/derived/north_atlantic/diagnostics_na_${MONTH}.zarr"

echo "=== NA diagnostics ${MONTH} on $(hostname) at $(date -Is) ==="

if [ -d "$OUT" ] && [ -f "$OUT/.zmetadata" ]; then
  echo "output already exists, skipping: $OUT"
  du -sh "$OUT" | sed 's/^/  /'
  exit 0
fi

if [ ! -f "$IN" ]; then
  echo "!! input not found: $IN"
  echo "   download it first:"
  echo "     sbatch --array=72,73,83,192,193,203,312,313,323,432,433,443%4 jobs/01_download.sbatch"
  exit 1
fi

mkdir -p "$BASE/derived/north_atlantic"

start=$(date +%s)
pixi run python ada/diagnostics_global.py "$IN" "$OUT" --chunk-days 0
status=$?
end=$(date +%s)

echo
echo "=== exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
exit $status
