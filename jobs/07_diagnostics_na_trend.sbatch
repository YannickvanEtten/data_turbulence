#!/bin/bash
# ===========================================================================
# 07_diagnostics_na_trend.sbatch -- the 21 diagnostics on the North Atlantic
# months of five DJF seasons: 1979, 1990, 2000, 2010, 2020.
#
#     sbatch jobs/07_diagnostics_na_trend.sbatch
#
# 15 array tasks. The three 1979 months already have output from
# 05_diagnostics_na_djf.sbatch; those tasks detect it and skip, so this can be
# submitted whole rather than with a hand-maintained index list. Re-running
# after a partial failure is therefore free and safe.
#
# Prerequisite: the twelve new months downloaded via
#     sbatch --array=132,133,143,252,253,263,372,373,383,492,493,503 \
#            jobs/01_download.sbatch
# (index = (year-1979)*12 + month-1 over 01_download.sbatch's 1979-01..2020-12
# span; no new download script needed.)
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
#SBATCH --array=0-14
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
  1990-01 1990-02 1990-12
  2000-01 2000-02 2000-12
  2010-01 2010-02 2010-12
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
  echo "     sbatch --array=132,133,143,252,253,263,372,373,383,492,493,503 jobs/01_download.sbatch"
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
