#!/bin/bash
# ===========================================================================
# 03_diagnostics.sbatch -- compute the 21 diagnostics for ONE month.
#
# This is the measurement run: the first time the diagnostics touch real ERA5
# through this pipeline on ADA. Its job is as much to tell us what a month
# COSTS -- peak memory, wall time, output size -- as to produce output.
#
# Submit (defaults to the pilot month, 1979-01):
#     sbatch jobs/03_diagnostics.sbatch
#
# Or any other month already downloaded:
#     sbatch jobs/03_diagnostics.sbatch 1979-02
#
# Afterwards, the number that matters most:
#     sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ReqMem
# MaxRSS is the real peak memory. --mem below is a guess; MaxRSS is the truth,
# and it is what the 504-month run should be sized from.
# ===========================================================================

#SBATCH --job-name=cat-diag
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/diag-%j.out
#SBATCH --error=logs/diag-%j.out

# --mem=64G is DELIBERATELY generous for a first run. A month is 248 timesteps
#          x 3 levels x 121 x 301, and 21 diagnostics on top of that, plus
#          rojak's intermediates. Rough arithmetic puts the diagnostics alone
#          near 4 GB at float64 before any intermediate is counted, and an OOM
#          kill would waste the run and tell us nothing. Over-asking once to
#          MEASURE is cheaper than guessing twice. Trim it for production once
#          MaxRSS is known.
# --time=04:00:00 likewise: no idea yet how long 248 timesteps takes.

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

MONTH="${1:-1979-01}"
IN="$BASE/raw/north_atlantic/era5_na_${MONTH}.grib"
OUT="$BASE/derived/north_atlantic/${MONTH}"

module load 2025
module load pixi
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== diagnostics for ${MONTH} on $(hostname) at $(date -Is) ==="
if [ ! -f "$IN" ]; then
  echo "!! input not found: $IN"
  echo "   downloaded months:"
  ls -lh "$BASE/raw/north_atlantic/" 2>/dev/null | sed 's/^/     /'
  exit 1
fi
echo "input : $IN  ($(du -h "$IN" | cut -f1))"
echo "output: $OUT"
mkdir -p "$OUT"

start=$(date +%s)
pixi run python 3_pipeline.py "$IN" "$OUT"
status=$?
end=$(date +%s)

echo
echo "=== RESULT ==="
echo "exit status : $status"
echo "wall time   : $(( (end-start)/60 )) min $(( (end-start)%60 )) s"

if [ $status -eq 0 ]; then
  echo
  echo "--- output size, which is the question this run exists to answer ---"
  du -sh "$OUT" 2>/dev/null | sed 's/^/  total   /'
  du -sh "$OUT"/*.nc "$OUT"/*.zarr 2>/dev/null | sed 's/^/  /'
  echo
  echo "  NOTE: 3_pipeline.save_outputs() writes the SAME data twice, once as"
  echo "  NetCDF and once as Zarr. For the 504-month run pick ONE -- that alone"
  echo "  halves the derived-data footprint. Zarr is the better choice for"
  echo "  chunked parallel reads later; NetCDF is easier to hand to someone."
  echo
  echo "--- extrapolated to the full run ---"
  bytes=$(du -sb "$OUT" 2>/dev/null | cut -f1)
  if [ -n "$bytes" ]; then
    echo "  both formats, x504 months : $(( bytes * 504 / 1000000000 )) GB"
    echo "  one format,   x504 months : $(( bytes * 504 / 2 / 1000000000 )) GB"
    echo "  (against a 2.5 TB share already holding ~183 GB of raw GRIB)"
  fi
  echo
  echo "--- the comparison table (sanity, not proof) ---"
  cat "$OUT/comparison_table.csv" 2>/dev/null
fi

echo
echo "=== done at $(date -Is) ==="
echo "Now check peak memory -- this is what to size the production run from:"
echo "    sacct -j \$SLURM_JOB_ID --format=JobID,State,Elapsed,MaxRSS,ReqMem"
exit $status
