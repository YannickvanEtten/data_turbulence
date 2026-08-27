#!/bin/bash
# ===========================================================================
# 02_download_global.sbatch -- the year-2000 GLOBAL calibration domain.
#
# This is a DIFFERENT dataset from the North Atlantic series, for a different
# purpose: Prosser calibrates severity thresholds on a global year and applies
# them to the analysis box. 3_pipeline.run_layers_2_to_6 takes the two as
# separate arguments precisely because they are not the same data.
#
# Submit:
#     sbatch jobs/02_download_global.sbatch
#
# RUN ONE MONTH FIRST. Set --array=0 and check it succeeds before letting all
# 12 go. A global month is ~10.6 GB and 5,208 fields in a single CDS request,
# and Copernicus enforces a per-request size limit whose current value we have
# not confirmed. If a whole month is rejected, split into daily requests (168
# fields each) -- download_plan.build_request already accepts a `days` list, so
# that is a small change to this script, not to the library.
#
# Sizing (measured, not estimated):
#     global grid  721 x 1440 = 1,038,240 points = 28.5x the North Atlantic box
#     one day      349 MB
#     one month    10.6 GB
#     year 2000    128 GB   (leap year, 366 days)
# ===========================================================================

#SBATCH --job-name=era5-glob
#SBATCH --partition=defq
#SBATCH --array=0-11%2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/global-%A_%a.out
#SBATCH --error=logs/global-%A_%a.out

# --- why these differ from 01_download.sbatch ------------------------------
# --array=0-11   12 months of the year 2000 only.
#
#        %2      Lower than the North Atlantic job's %4. Each of these requests
#                is ~28x larger, so two at a time is already far more load on
#                CDS -- and on this end, two 10.6 GB transfers at once.
#
# --mem=16G      NOT 4G. The integrity check opens the finished GRIB and counts
#                dimensions; a global month is 28.5x the data. cfgrib is lazy,
#                but the safety margin is cheap and an OOM here wastes the
#                entire download.
#
# --time=24:00:00  A 10.6 GB transfer plus CDS queueing. Still well inside the
#                7-day cap.
# ---------------------------------------------------------------------------

set -euo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
OUT=$BASE/raw/global

DOMAIN=global
START=2000-01
END=2000-12

module load 2025
module load pixi

mkdir -p "$OUT"
cd "$REPO"   # pixi.toml and .pixi/ live here

echo "=== task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -Is) ==="
df -h "$BASE" | tail -1

pixi run python 1_download_hpc.py \
    --start "$START" \
    --end   "$END" \
    --domain "$DOMAIN" \
    --index "$SLURM_ARRAY_TASK_ID" \
    --out-dir "$OUT"

echo "=== done at $(date -Is) ==="
