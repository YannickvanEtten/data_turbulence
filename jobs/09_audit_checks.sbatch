#!/bin/bash
# ===========================================================================
# 09_audit_checks.sbatch -- the two CHEAP checks from FORMULA_AUDIT.md, in one
# job: pointwise composition identities (§7.3) and the p97/median shape test
# across all 21 (§7.1).
#
#     mkdir -p logs && sbatch jobs/09_audit_checks.sbatch
#
# `mkdir -p logs` goes on the COMMAND LINE, not in this script: SLURM opens the
# --output file before the script body runs, so creating the directory inside
# would be too late (STATUS.md §11.9).
#
# Neither check needs a new CDS request. Both run on the 39 months already on
# disk, so this can be submitted while jobs/01_download.sbatch is still
# fetching the production series -- the raw download is diagnostic-independent
# and nothing here can invalidate a byte of it.
#
# CONCURRENCY. The default QOS caps you at 8 RUNNING JOBS (STATUS.md §11.9),
# and a download array at %4 is holding four of them. That leaves four, which
# is exactly enough for 12_tests plus 09, 10 and 11. Nothing here needs to wait
# for the download to finish.
#
# WHAT EACH ONE DECIDES
# ---------------------
# check_composition   Verifies nine of the 21 against each other, cell by cell,
#                     with no external reference. Catches a mis-wired argument
#                     or a chunk-boundary error -- faults that produce a
#                     plausible field and survive every distributional test.
#                     Also settles empirically whether the zarr on disk holds
#                     DEF or DEF^2 (FORMULA_AUDIT.md §5). Exit 2 means an
#                     identity FAILED; that is a stop sign, not a warning.
#
# check_shape_ratios  Turns the two Williams tables from an untransferable
#                     magnitude comparison into a dimensionless shape test
#                     covering all 21 at once. Predicted in advance to flag
#                     f2d; see FORMULA_AUDIT.md §7.1.
#
# --mem=16G   Measured shapes, not guesses. One NA month is 121 x 301 x 248 =
#             9.0e6 cells; eleven fields as float64 plus the derived arrays is
#             under 2 GB. The shape test subsets to Williams' box BEFORE
#             loading, so three DJF months of the global zarr resident at once
#             is a few hundred MB. STATUS.md §11.9 measured a 120 GB request
#             scheduled 22 hours out against 13 seconds for 24 GB -- memory
#             footprint is the biggest lever on time-to-start here, so this is
#             deliberately not rounded up.
# --time=01:00  Both are minutes of real work. Short walltimes backfill.
# ===========================================================================

#SBATCH --job-name=cat-audit
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/audit-%j.out
#SBATCH --error=logs/audit-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== FORMULA_AUDIT checks on $(hostname) at $(date -Is) ==="
echo "    repo $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
start=$(date +%s)
rc=0

# --- pick a month rather than hard-coding one ------------------------------
# The identities hold POINTWISE, so any month tests them exactly as well as any
# other and more months add confidence about nothing. Discovering the file
# instead of naming it also avoids a stale assumption: STATUS.md §13.1 records
# derived/north_atlantic/1979-01/ being deleted on 2026-08-29, and a hard-coded
# path is exactly how a script starts failing for a reason that has nothing to
# do with what it tests.
MONTH=""
for CAND in "$BASE"/derived/north_atlantic/diagnostics_na_*.zarr; do
    if [ -d "$CAND" ]; then MONTH="$CAND"; break; fi
done

echo
echo "############################################################"
echo "# 1/3  POINTWISE COMPOSITION IDENTITIES  (FORMULA_AUDIT 7.3)"
echo "############################################################"
if [ -n "$MONTH" ]; then
    echo "using $MONTH"
    pixi run python ada/check_composition.py "$MONTH"
    rc=$(( rc | $? ))
else
    echo "!! no diagnostics_na_*.zarr under $BASE/derived/north_atlantic"
    echo "   (run jobs/07_diagnostics_na_trend.sbatch first)"
    rc=$(( rc | 1 ))
fi

# --- 2. shape ratios, daily-mean (matches Williams' sampling) ---------------
echo
echo "############################################################"
echo "# 2/3  SHAPE RATIOS, DAILY MEAN  (FORMULA_AUDIT 7.1)"
echo "############################################################"
pixi run python ada/check_shape_ratios.py --base "$BASE" --box williams
rc=$(( rc | $? ))

# --- 3. the same, instantaneous --------------------------------------------
# Williams computes his diagnostics from daily-mean WIND AND TEMPERATURE; step
# 2 averages the DIAGNOSTIC, which is the closer of the two approximations
# available to us but not the same operation. The gap between step 2 and step 3
# bounds how much that approximation can be worth, which is the honest way to
# report it rather than arguing about its size.
echo
echo "############################################################"
echo "# 3/3  SHAPE RATIOS, 3-HOURLY INSTANTANEOUS  (the caveat)"
echo "############################################################"
pixi run python ada/check_shape_ratios.py --base "$BASE" --box williams --no-daily-mean
rc=$(( rc | $? ))

end=$(date +%s)
echo
echo "=== exit ${rc}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
echo
echo "Read in this order:"
echo "  1. every identity PASS at ~1e-7      -> nine of the 21 verified against"
echo "                                          each other, brown2 fully closed"
echo "  2. whether deformation reads DEF^2   -> expected on existing zarr"
echo "  3. f2d's shape ratio vs 13.6         -> hundreds refutes variant A"
exit $rc
