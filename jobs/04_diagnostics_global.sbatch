#!/bin/bash
# ===========================================================================
# 04_diagnostics_global.sbatch -- the 21 diagnostics on the year-2000 global
# calibration files.
#
# RUN ONE MONTH FIRST.
# --------------------
#     sbatch --array=0 jobs/04_diagnostics_global.sbatch
#
# MEASURED on 2000-01 (job 1092495, node009):
#
#     peak RSS    17.7 GB      <-- against --mem=24G, 74% used
#     wall time   30 min 18 s
#     output      2.28 GB zarr, 21 of 21 diagnostics, no failures
#
# 24G is therefore right and stays: it fits EVERY defq node including the
# 59 GB ones, which is what let this schedule in 13 seconds. Raising it would
# buy nothing and cost scheduling speed -- a 120G request on this cluster was
# queued 22 h out.
#
# The driver prints its own PEAK RSS. Do not chase sacct or sstat: this
# cluster has job-step accounting disabled and both return empty.
#
# --chunk-days 1 COSTS NOTHING
# ----------------------------
# Each chunk is computed with one extra timestep on each side and trimmed
# back, so the output is bit-identical to processing the whole file at once --
# verified across 1-, 2-, 4-day and whole-file chunking. F2D loses exactly the
# 2 steps at the file's true ends either way; NAIVE chunking without the
# overlap would cost it 8. See ada/diagnostics_global.py for the reasoning.
#
# IF IT OOMs
# ----------
# The job is killed with an out-of-memory message rather than producing a
# wrong answer. Raise --mem; the chunk size is already at its minimum useful
# value, so an OOM at 24G would mean the per-chunk footprint is much larger
# than the North Atlantic month's and is worth understanding before the
# 504-month run.
#
# OUTPUT
#     derived/global/diagnostics_glob_2000-MM.zarr    2.28 GB each, ~27 GB total
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
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/diag-glob-%A_%a.out
#SBATCH --error=logs/diag-glob-%A_%a.out

# --array=0-11   no % throttle: this is CPU-bound local work, not CDS requests.
# --mem=24G      With --chunk-days 1 the peak is set by ONE DAY (721x1440x8 =
#                8.3e6 points, essentially the North Atlantic month's footprint
#                at its measured 14.8 GB), not by the whole file. 24G therefore
#                fits EVERY defq node including the 59 GB ones -- which is the
#                difference between starting now and being scheduled tomorrow.
#                A 120G request on this cluster was queued 22 h out while every
#                node sat in 'mix'; the 4-8 GB download jobs started instantly.
# --time=02:00   MEASURED: 2000-01 took 30 min 18 s on node009. 2 h is ~4x
#                that, which covers a slower node (defq spans several
#                generations) without being so long that backfill ignores it.
#                1 h -- my first guess -- would have left only 2x margin.

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
DAYS=d01-09-17-25

module load 2025
module load pixi

# Python block-buffers stdout when it is a file rather than a terminal, so
# without this NOTHING appears in the log until the process exits -- the
# per-chunk progress lines are useless exactly when you want them, and a job
# that is working is indistinguishable from one that is hung. stderr is
# unbuffered either way, so a traceback would still appear immediately.
export PYTHONUNBUFFERED=1

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
pixi run python ada/diagnostics_global.py "$IN" "$OUT" --chunk-days 1
status=$?
end=$(date +%s)

echo
echo "=== exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
echo "Peak memory is the PEAK RSS line above -- sacct on this cluster has not"
echo "been returning MaxRSS for job steps. For a RUNNING job, sstat works:"
echo "    sstat -j \$SLURM_JOB_ID --format=JobID,MaxRSS,AveCPU"
exit $status
