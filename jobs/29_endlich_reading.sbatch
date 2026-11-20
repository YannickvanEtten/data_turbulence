#!/bin/bash
# ---------------------------------------------------------------------------
# 29_endlich_reading.sbatch                        (ENDLICH-READING-2026-11-20)
#
# Which reading of Sharman A25 should `endlich` use as a separate series?
#   E0 = angle difference of the wind direction across 175/225 hPa (rojak,
#        used by BOTH production series), E1 = closed form
#        (u dv/dz - v du/dz)/|V|^2 (audit F10). Only endlich is computed.
# Rationale and pre-registered Q1-Q4: docstring of ada/endlich_test.py.
#
# WHAT IS ALREADY KNOWN
#   - E0 and E1 reproduce essentials.endlich exactly (checked locally).
#   - On the July-2016 validation subset the readings rank the bulk alike
#     (Spearman 0.96) but swap 69-75 % of the p97 / p99.6 cells, so the tail
#     WILL move. The open question is whether it moves toward Prosser.
#
# READ IN THIS ORDER
#   1. Stage 1: Spearman >= 0.99999, else exit 2.
#   2. Stage 2's "box density, E1 / E0" line against Q3.
#   3. Stage 3's implied fitted ratios against the stencil-group median
#      0.93 / 0.85 / 0.83 and the class-rule band 0.75-1.33.
#
# INPUTS as jobs/28 (no new data). SIZING as jobs/28: two fields, needs all
# three levels per day chunk; 32G, 4 h walltime, expected < 2 h. Can run at the
# same time as jobs/28 (both only read the raw files).
#
#   git pull && grep -l ENDLICH-READING-2026-11-20 ada/endlich_test.py jobs/29_endlich_reading.sbatch
#   mkdir -p logs && sbatch jobs/29_endlich_reading.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-endlich
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/endlich-%j.out
#SBATCH --error=logs/endlich-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

if [ "${SLURM_NTASKS:-1}" != "1" ]; then
    echo "FATAL: SLURM_NTASKS=${SLURM_NTASKS}"; exit 1
fi
if ! grep -q "ENDLICH-READING-2026-11-20" ada/endlich_test.py; then
    echo "FATAL: ada/endlich_test.py is not the 2026-11-20 version -- was it pushed and pulled?"
    exit 2
fi

echo "=== endlich reading test on $(hostname) at $(date -Is) ==="
echo "    git rev $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
sha256sum ada/endlich_test.py jobs/29_endlich_reading.sbatch
echo "    global sub-sample files: $(ls $BASE/raw/global_sub48/*.grib 2>/dev/null | wc -l) (need 12)"
echo "    NA 1979 files:           $(ls $BASE/raw/north_atlantic/era5_na_1979-*.grib 2>/dev/null | wc -l) (need 12)"
echo

s=$(date +%s)
pixi run python ada/endlich_test.py --base "$BASE" --out cat_outputs/endlich_reading
st=$?
echo
echo "=== exit $st after $(( ($(date +%s) - s) / 60 )) min at $(date -Is) ==="
echo "Outputs: $REPO/cat_outputs/endlich_reading/{summary.json,global_tail.csv}"
echo "Bring back: the log, summary.json and global_tail.csv (git add -f if ignored)."
exit $st
