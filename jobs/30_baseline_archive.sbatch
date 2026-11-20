#!/bin/bash
# ---------------------------------------------------------------------------
# 30_baseline_archive.sbatch                     (BASELINE-ARCHIVE-2026-11-20)
#
# Run ONCE, BEFORE derived/north_atlantic/ (the baseline series, ~294 GB) is
# deleted. Writes the baseline's per-diagnostic and ensemble annual series at
# all five severities -- a few hundred kB of CSV -- so the convention
# sensitivity (baseline vs audit) stays reproducible as numbers after the
# stores are gone. Same driver and same output format as the audit tables of
# job 1124312, so the two can be compared column for column.
#
# Nothing is overwritten: output goes to cat_outputs/baseline_archive/.
# The driver refuses a mismatched series/threshold pairing (exit 2).
#
#   mkdir -p logs && sbatch jobs/30_baseline_archive.sbatch
#   afterwards: git add -f cat_outputs/baseline_archive logs/baseline-archive-<jobid>.out
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-base-archive
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:45:00
#SBATCH --output=logs/baseline-archive-%j.out
#SBATCH --error=logs/baseline-archive-%j.out

set -uo pipefail
BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
module load 2025
module load pixi
export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== baseline archive on $(hostname) at $(date -Is), git rev $(git rev-parse --short HEAD 2>/dev/null) ==="
for p in derived/north_atlantic calibration/thresholds_2026-09-07.json; do
    [ -e "$BASE/$p" ] || { echo "!! $BASE/$p missing -- nothing to archive"; exit 2; }
done
echo "    baseline stores: $(ls -d $BASE/derived/north_atlantic/*.zarr | wc -l) (need 504)"

pixi run python ada/per_diagnostic_all_severities.py --base "$BASE" \
    --derived-subdir derived/north_atlantic \
    --thresholds "$BASE/calibration/thresholds_2026-09-07.json" \
    --outdir cat_outputs/baseline_archive --tag _baseline
st=$?
echo "=== exit $st at $(date -Is) ==="
ls -la cat_outputs/baseline_archive/
exit $st
