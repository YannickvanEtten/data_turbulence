#!/bin/bash
# ===========================================================================
# 05_diagnostics_na_djf.sbatch -- the 21 diagnostics on North Atlantic
# DJF 1979, the ANALYSIS side of the calibration check.
#
#     sbatch jobs/05_diagnostics_na_djf.sbatch      # all three at once
#
# WHY NOT 03_diagnostics.sbatch?
# ------------------------------
# 03 runs 3_pipeline.run(), which writes the same data as NetCDF *and* Zarr
# in float64 (~2.4 GB per month for output the comparison reads once) and
# builds a W&J comparison table per month. That is the right shape for the
# production analysis run, where those artefacts are the point.
#
# Here the output feeds exactly one consumer -- the threshold comparison --
# and it must be directly comparable to the global calibration fields, which
# ada/diagnostics_global.py wrote as float32 zarr. Running BOTH sides of the
# comparison through the SAME code path removes a whole class of "is this
# difference real or is it the two writers disagreeing" question.
#
# The driver is not global-specific despite its name; it takes any ERA5 GRIB.
#
# --chunk-days 0 (whole file)
# ---------------------------
# A North Atlantic month is 121 x 301 x 248 = 9.0e6 points, and the whole-file
# footprint was measured at 14.8 GB -- comfortably inside 24G, and the same
# ballpark as one global chunk's 17.7 GB. No chunking needed, so F2D loses
# only the 2 timesteps at each month's ends.
#
# NOTE ON THE MONTH BOUNDARIES
# ----------------------------
# DJF 1979 is assembled from three separately-downloaded months, so F2D loses
# its first and last timestep in EACH of them -- 6 of 744, or 0.8%. Closing
# that properly is what chunk_stitch.py exists for at the month level; it is
# not worth wiring in for a check whose acceptance criterion is a frequency
# to two significant figures. Record it, do not fix it here.
#
# DJF 1979 = Jan + Feb + Dec 1979. Prosser's record starts 1 Jan 1979, so
# there is no Dec 1978 to use; 90 days = 2160 hours is the denominator the
# comparison divides by.
# ===========================================================================

#SBATCH --job-name=cat-diag-na
#SBATCH --partition=defq
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/diag-na-%A_%a.out
#SBATCH --error=logs/diag-na-%A_%a.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

MONTHS=(1979-01 1979-02 1979-12)
MONTH=${MONTHS[$SLURM_ARRAY_TASK_ID]}

module load 2025
module load pixi
export PYTHONUNBUFFERED=1        # see 04_diagnostics_global.sbatch

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

IN="$BASE/raw/north_atlantic/era5_na_${MONTH}.grib"
OUT="$BASE/derived/north_atlantic/diagnostics_na_${MONTH}.zarr"

echo "=== NA diagnostics ${MONTH} on $(hostname) at $(date -Is) ==="
if [ ! -f "$IN" ]; then
  echo "!! input not found: $IN"
  ls -lh "$BASE/raw/north_atlantic/" 2>/dev/null | sed 's/^/     /'
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
