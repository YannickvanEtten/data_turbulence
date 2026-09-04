#!/bin/bash
# ---------------------------------------------------------------------------
# 17_download_global_fullyear.sbatch
#
# The FULL global year 2000, 3-hourly, for the reference-year calibration.
# Supersedes 02b (days 1/9/17/25) as the calibration input; 02b is kept, since
# comparing the two is what turns "we sub-sampled" from an unexamined choice
# into a measured sensitivity.
#
# WHY 36 TASKS AND NOT 12
#     A full global month is 31 x 8 x 3 x 7 = 5,208 fields and ~10.6 GB in ONE
#     CDS request. The Copernicus per-request limit has never been confirmed
#     (STATUS 9a, "Open sizing question") -- it was sidestepped, not tested.
#     Three ~10-day blocks per month is 1,680 fields and ~3.4 GB, comfortably
#     inside anything plausible, and it parallelises better at the %4 ceiling.
#     jobs/18 concatenates the blocks back into contiguous months.
#
# WHY 3-HOURLY AND NOT FINER
#     Not a storage decision. The thresholds must be computed on the SAME
#     temporal sampling as the series they are applied to, or f2d's threshold
#     comes from an undamped distribution and is applied to a damped one
#     (STATUS 4d: a diurnal signal at 3-hourly spacing is understated by 36%).
#     The North Atlantic production series is 3-hourly, so this is too.
#
# WHY 3 LEVELS AND 7 VARIABLES
#     Both LOCKED in download_plan.py. 175/200/225 is the finest stencil that
#     exists at 200 hPa on the CDS pressure-level product (STATUS 11.6), and
#     all seven variables are load-bearing -- 'd' in particular, because ti2
#     uses ERA5's own divergence rather than the wind-derived one (STATUS 4c).
#
# ⚠ SUBMIT BARE.  `sbatch --array=... jobs/17_...` REPLACES the entire array
#   spec including the `%4`, silently. That has now cost this project three
#   separate incidents (STATUS 11.9, 14.7) -- 8-of-12 lost once, ~40-of-55
#   lost once. If an override is unavoidable, re-append the throttle:
#       sbatch --array=3,4,5%4 jobs/17_download_global_fullyear.sbatch
#
# SIZING          ~128 GB raw GRIB across 36 files (measured 2.00 bytes per
#                 gridpoint, STATUS 10.6; the same constant predicts the
#                 existing sub-sample at 16.7 GB against 16 GB on disk).
# TIMING          ~36 requests at %4; global months measured 17 min - 23 min
#                 at 1/8 this size, so budget 8-12 h wall.
#
#   mkdir -p logs && sbatch jobs/17_download_global_fullyear.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=era5-glob-year
#SBATCH --partition=defq
#SBATCH --array=0-35%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/glob-year-%A_%a.out
#SBATCH --error=logs/glob-year-%A_%a.out

set -euo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
OUT=$BASE/raw/global
DOMAIN=global
YEAR=2000

# 3 day-blocks per month. The last block runs to the month's REAL last day --
# February 2000 has 29, and passing 30/31 would make the integrity check
# expect timesteps the file cannot contain and fail a perfectly good download.
MONTH_IDX=$(( SLURM_ARRAY_TASK_ID / 3 ))
BLOCK=$(( SLURM_ARRAY_TASK_ID % 3 ))
MONTH=$(printf "%02d" $(( MONTH_IDX + 1 )))
NDAYS=$(date -d "${YEAR}-${MONTH}-01 +1 month -1 day" +%d)
case $BLOCK in
  0) LO=1;  HI=10 ;;
  1) LO=11; HI=20 ;;
  2) LO=21; HI=$NDAYS ;;
esac
DAYS=$(seq -w "$LO" "$HI" | paste -sd, -)

module load 2025
module load pixi          # lowercase -- 'module load Pixi' fails (STATUS 10.2)

export PYTHONUNBUFFERED=1 # else a working job is indistinguishable from a hung
                          # one until it exits (STATUS 11.9)

mkdir -p "$OUT" "$REPO/logs"
cd "$REPO"                # pixi.toml and .pixi/ live here, not in $BASE

echo "=== task ${SLURM_ARRAY_TASK_ID}: ${YEAR}-${MONTH} days ${LO}-${HI} on $(hostname) at $(date -Is) ==="
echo "SLURM_NTASKS=${SLURM_NTASKS:-unset}   (must be 1)"
if [ "${SLURM_NTASKS:-1}" != "1" ]; then
    # Two tasks on one node means two processes requesting the same block from
    # CDS and racing on the same .tmp. The atomic-rename guard protects against
    # a half-written file, NOT against two concurrent writers (STATUS 10.7).
    echo "FATAL: SLURM_NTASKS is not 1 -- refusing to run." >&2
    exit 1
fi
df -h "$BASE" | tail -1

pixi run python 1_download_hpc.py \
    --start  "${YEAR}-01" \
    --end    "${YEAR}-12" \
    --domain "$DOMAIN" \
    --days   "$DAYS" \
    --index  "$MONTH_IDX" \
    --out-dir "$OUT"

echo "=== done at $(date -Is) ==="
ls -lh "$OUT" | tail -5
