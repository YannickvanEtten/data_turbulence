#!/bin/bash
# ===========================================================================
# 04_diagnostics_global.sbatch -- the 21 diagnostics on the year-2000 global
# calibration files.
#
# RUN ONE MONTH FIRST, AND MEAN IT THIS TIME.
# ------------------------------------------
#     sbatch --array=0 jobs/04_diagnostics_global.sbatch
#
# The only MaxRSS measurement this project has is 14.8 GB, from a North
# Atlantic month of shape 121 x 301 x 248 (9.0e6 points). A global calibration
# file is 721 x 1440 x 32 -- 3.3e7 points, 3.7x more, but a completely
# different array SHAPE. Linear extrapolation says ~55 GB; that is a guess
# across a shape change, and rojak's gradient operators need not scale with
# point count alone. The first task exists to replace that guess with a
# number.
#
# --mem below is deliberately generous for the FIRST run. After it finishes:
#
#     sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS
#
# then set --mem to roughly 1.5x the real MaxRSS and release the other eleven.
# Memory is scheduled explicitly on ADA, so an honest request is what decides
# how many array tasks fit on a node -- exactly the reasoning that took
# 03_diagnostics.sbatch from 64G down to 24G.
#
# IF IT OOMs
# ----------
# The job will be killed with an out-of-memory message rather than producing
# a wrong answer. Two routes out, in order of preference:
#   1. raise --mem and, if needed, switch to --partition=defq-fat
#   2. fall back to processing each day separately, at the cost of F2D losing
#      2 of every 8 timesteps to its centred material derivative instead of
#      2 of every 32
#
# OUTPUT
#     derived/global/diagnostics_glob_2000-MM.zarr    ~2.8 GB each, ~33 GB total
#
# float32 and zarr only -- NOT 3_pipeline.run(), which writes NetCDF AND zarr
# in float64 (~134 GB for twelve months) and builds a W&J comparison table
# that is not meaningful on a global domain. Layer 2 needs the values and
# their latitudes; that is all this writes.
# ===========================================================================

#SBATCH --job-name=cat-diag-glob
#SBATCH --partition=defq
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/diag-glob-%A_%a.out
#SBATCH --error=logs/diag-glob-%A_%a.out

# --array=0-11   no % throttle: this is CPU-bound local work, not CDS requests.
#                But note --mem=96G means few tasks fit per node, so SLURM will
#                serialise them anyway until the request is right-sized.
# --mem=96G      A DELIBERATE OVERESTIMATE for task 0 only. Replace it with
#                ~1.5x the measured MaxRSS before running the rest.
# --time=06:00   The NA month took 8 min 11 s. This file has 3.7x the points,
#                so expect well under an hour; the rest is slack for a stall.

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
DAYS=d01-09-17-25

module load 2025
module load pixi

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

MONTH=$(printf "%02d" $((SLURM_ARRAY_TASK_ID + 1)))
IN="$BASE/raw/global/era5_glob_2000-${MONTH}_${DAYS}.grib"
OUT="$BASE/derived/global/diagnostics_glob_2000-${MONTH}.zarr"

echo "=== global diagnostics 2000-${MONTH} on $(hostname) at $(date -Is) ==="
if [ ! -f "$IN" ]; then
  echo "!! input not found: $IN"
  ls -lh "$BASE/raw/global/" 2>/dev/null | sed 's/^/     /'
  exit 1
fi

mkdir -p "$BASE/derived/global"

start=$(date +%s)
pixi run python ada/diagnostics_global.py "$IN" "$OUT"
status=$?
end=$(date +%s)

echo
echo "=== exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
echo "Peak memory -- size the remaining eleven from this:"
echo "    sacct -j \$SLURM_JOB_ID --format=JobID,State,Elapsed,MaxRSS"
exit $status
