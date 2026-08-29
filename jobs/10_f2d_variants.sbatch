#!/bin/bash
# ===========================================================================
# 10_f2d_variants.sbatch -- which of the four readings of Sharman A9 is #20?
# FORMULA_AUDIT.md §4.
#
#     sbatch jobs/10_f2d_variants.sbatch
#
# This is STEP 1 of two. It measures which variant is DISTRIBUTED like the
# published diagnostic (Williams 2017 p97 / W&J 2013 median = 13.6). It cannot
# separate variant A from variant B, because B = -A and every symmetric
# statistic agrees -- what separates those two is the SIGN OF THE TREND against
# Prosser Figure S5, which needs a full re-run and is step 2:
#
#     # step 2, only once step 1 has narrowed it:
#     sbatch --array=0-11%8 jobs/04_diagnostics_global.sbatch   # + --f2d-variant X
#     sbatch jobs/06_calibration_check.sbatch
#     sbatch --array=0-26%8 jobs/07_diagnostics_na_trend.sbatch # + --f2d-variant X
#     sbatch jobs/08_trend_check.sbatch --per-diagnostic
#
# Run on a GLOBAL calibration month, because Williams' published pair is for
# 50-75N/10-60W and only the global domain reaches past 60N -- the North
# Atlantic download stops there.
#
# --mem=48G   Recomputes f2d twice (variants A and D are genuinely different
#             computations; B and C are exact functions of A) on a global month
#             at three levels. Higher than the 24G that the full 21-diagnostic
#             global run needed, because both variants are held at once before
#             subsetting. Still well inside the 59G nodes.
# --time=02:00  A global month of diagnostics took 15-36 min (STATUS.md §13.3);
#             this computes a small subset of them twice.
# ===========================================================================

#SBATCH --job-name=cat-f2d
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/f2d-%j.out
#SBATCH --error=logs/f2d-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== f2d variant check on $(hostname) at $(date -Is) ==="
start=$(date +%s)

# January of the calibration year: DJF, and the month whose global diagnostics
# every threshold in calibration/thresholds_*.json already descends from.
GLOB=$BASE/raw/global/era5_glob_2000-01_d01-09-17-25.grib

if [ ! -f "$GLOB" ]; then
    echo "!! $GLOB not found"
    echo "   available:"
    ls -1 "$BASE"/raw/global/*.grib 2>/dev/null | head
    exit 1
fi

echo
echo "############################################################"
echo "# GLOBAL month, Williams box 50-75N/10-60W"
echo "############################################################"
pixi run python ada/check_f2d_variants.py "$GLOB" --box williams
rc=$?

# The same question asked in the box the replication actually reports on. The
# shape verdict should not depend on which box it is asked in; if it does, that
# is itself a finding and it belongs in the write-up.
echo
echo "############################################################"
echo "# The same file, Prosser box 36-60N/55-10W"
echo "############################################################"
pixi run python ada/check_f2d_variants.py "$GLOB" --box prosser
rc=$(( rc | $? ))

end=$(date +%s)
echo
echo "=== exit ${rc}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
exit $rc
