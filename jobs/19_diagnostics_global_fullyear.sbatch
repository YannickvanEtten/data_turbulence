#!/bin/bash
# ---------------------------------------------------------------------------
# 19_diagnostics_global_fullyear.sbatch
#
# All 21 diagnostics on the merged FULL global months from jobs/18.
# Identical driver to jobs/04 (ada/diagnostics_global.py, --chunk-days 1,
# float32, zarr only, 200 hPa). Three things differ, and one of them is the
# one that bites:
#
#   1. INPUT has no _${DAYS} suffix -- these are whole months now.
#   2. --time is 08:00:00, up from jobs/04's 02:00:00.  ** THIS IS THE CHANGE
#      THAT MATTERS. ** A full month is 7.62x the data of a 4-day month, and
#      those took 15-36 min. At 2 h, roughly every task would hit walltime and
#      die most of the way through -- the most expensive possible failure,
#      because a killed task leaves nothing and costs the whole run.
#   3. --mem is 32G, up from 24G. --chunk-days 1 still sets the peak by ONE
#      DAY (721x1440x8), so this is headroom for opening a 31-day GRIB rather
#      than a real increase in the working set. Do NOT round it up further:
#      a 120 GB request was scheduled 22 hours out while 24 GB started in
#      13 seconds (STATUS 11.9).
#
# OUTPUT      derived/global/diagnostics_glob_2000-MM.zarr, ~17 GB each,
#             ~200-255 GB for the year.
#
#             KEEP IT. The full field is the only global, all-season,
#             all-21-diagnostic object this project will have, and it is what
#             answers questions the tails cannot: per-diagnostic spatial
#             pathology, the rank correlations between the 21, and above all
#             ubf's latitude structure -- ubf is the one diagnostic behaving
#             anomalously (STATUS 12.6), its known failure mode scales with
#             tan(phi) (STATUS 4g), and a global field is exactly the picture
#             that would show it. Tails would hide it.
#
# ⚠ MOVE THE 48-DAY DERIVED OUTPUT ASIDE FIRST, or these overwrite it:
#       mv $BASE/derived/global $BASE/derived/global_sub48
#       mkdir -p $BASE/derived/global
#   Keep it: ada/tail_thresholds.py --derived-subdir derived/global_sub48
#   --validate-against calibration/thresholds_2026-09-03.json is the free
#   correctness check on the whole tail path, and it needs those stores.
#
# ⚠ jobs/15_check_production_complete.sbatch expects 12 stores in
#   derived/global and will report NOT CLEAR between the move and the refill.
#   That is expected, not a regression.
#
# ⚠ RUN ONE TASK FIRST:
#       sbatch --array=0 jobs/19_diagnostics_global_fullyear.sbatch
#   Check the peak RSS and elapsed it prints, THEN submit the rest. There is
#   no %N throttle to lose here (this job makes no CDS requests), so an
#   --array override is safe in this script and only in this script.
#
#   sbatch jobs/19_diagnostics_global_fullyear.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-diag-glob-yr
#SBATCH --partition=defq
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/diag-glob-yr-%A_%a.out
#SBATCH --error=logs/diag-glob-yr-%A_%a.out

set -euo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
YEAR=2000

module load 2025
module load pixi

export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

MONTH=$(printf "%02d" $(( SLURM_ARRAY_TASK_ID + 1 )))
IN="$BASE/raw/global/era5_glob_${YEAR}-${MONTH}.grib"
OUT="$BASE/derived/global/diagnostics_glob_${YEAR}-${MONTH}.zarr"

echo "=== global diagnostics ${YEAR}-${MONTH} (full month) on $(hostname) at $(date -Is) ==="
if [ ! -f "$IN" ]; then
  echo "!! input not found: $IN"
  echo "!! did jobs/18_merge_global_months.sbatch run?"
  ls -lh "$BASE/raw/global/" 2>/dev/null | tail -20 | sed 's/^/     /'
  exit 1
fi
ls -lh "$IN" | sed 's/^/     /'

mkdir -p "$BASE/derived/global"

start=$(date +%s)
pixi run python ada/diagnostics_global.py "$IN" "$OUT" --chunk-days 1
status=$?
end=$(date +%s)

echo "=== ${YEAR}-${MONTH} exit=$status in $(( (end - start) / 60 )) min at $(date -Is) ==="
du -sh "$OUT" 2>/dev/null | sed 's/^/     /'
df -h "$BASE" | tail -1
exit $status
