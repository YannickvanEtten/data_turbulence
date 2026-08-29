#!/bin/bash
# ===========================================================================
# 10_f2d_variants.sbatch -- which of the four readings of Sharman A9 is #20?
# FORMULA_AUDIT.md §4.
#
#     mkdir -p logs && sbatch jobs/10_f2d_variants.sbatch
#
# This is STEP 1 of two. It measures which variant is DISTRIBUTED like the
# published diagnostic (Williams 2017 p97 / W&J 2013 median = 13.6). It cannot
# separate variant A from variant B, because B = -A and every symmetric
# statistic agrees -- what separates those two is the SIGN OF THE TREND against
# Prosser Figure S5, which needs a full re-run and is step 2:
#
#     # step 2, only once step 1 has narrowed it. Both diagnostics scripts need
#     # --f2d-variant X adding to their ada/diagnostics_global.py line first;
#     # that edit is deliberately NOT made yet, so the scripts that produced
#     # the current results stay untouched until there is a reason.
#     sbatch --array=0-11%8  jobs/04_diagnostics_global.sbatch
#     sbatch                 jobs/06_calibration_check.sbatch
#     sbatch --array=0-26%8  jobs/07_diagnostics_na_trend.sbatch
#     sbatch                 jobs/08_trend_check.sbatch --per-diagnostic
#
# Run on a GLOBAL calibration month, because Williams' published pair is for
# 50-75N/10-60W and only the global domain reaches past 60N -- the North
# Atlantic download stops there.
#
# --mem=24G   check_f2d_variants.py chunks over time at --chunk-days 1, the
#             same device ada/diagnostics_global.py uses, so peak memory is set
#             by one day rather than by the file. 24 GB is the request that
#             STATUS.md §13.3 measured starting in 13 seconds while a 120 GB
#             one sat 22 hours in the queue. Do not raise it "to be safe" --
#             that costs a day, and the chunking is what makes it unnecessary.
# --time=02:00  A global month of all 21 took 15-36 min (STATUS.md §13.3); this
#             computes two variants of one of them.
# ===========================================================================

#SBATCH --job-name=cat-f2d
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
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
echo "    repo $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
start=$(date +%s)

# Discover the calibration month rather than naming it. The sub-sampled global
# files carry a day suffix (era5_glob_2000-01_d01-09-17-25.grib, STATUS.md
# §11.11) and a hard-coded name would break the moment that sampling changes.
# January of the calibration year: DJF, and the month every threshold in
# calibration/thresholds_*.json already descends from.
GLOB=""
for CAND in "$BASE"/raw/global/era5_glob_2000-01*.grib; do
    if [ -f "$CAND" ]; then GLOB="$CAND"; break; fi
done
if [ -z "$GLOB" ]; then
    for CAND in "$BASE"/raw/global/*.grib; do
        if [ -f "$CAND" ]; then GLOB="$CAND"; break; fi
    done
fi

if [ -z "$GLOB" ]; then
    echo "!! no global GRIB under $BASE/raw/global"
    exit 1
fi
echo "using $GLOB"

echo
echo "############################################################"
echo "# GLOBAL month, Williams box 50-75N/10-60W"
echo "############################################################"
pixi run python ada/check_f2d_variants.py "$GLOB" --box williams --chunk-days 1
rc=$?

# The same question asked in the box the replication actually reports on. The
# shape verdict should not depend on which box it is asked in; if it does, that
# is itself a finding and it belongs in the write-up.
echo
echo "############################################################"
echo "# The same file, Prosser box 36-60N/55-10W"
echo "############################################################"
pixi run python ada/check_f2d_variants.py "$GLOB" --box prosser --chunk-days 1
rc=$(( rc | $? ))

end=$(date +%s)
echo
echo "=== exit ${rc}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
echo
echo "Read the p97/median column against 13.6, and the tail-overlap matrix."
echo "A vs B overlapping near 0.000 is expected and is the point: the sign"
echo "choice selects two disjoint populations of grid cells, so it cannot be"
echo "left to a default."
exit $rc
