#!/bin/bash
# ===========================================================================
# 09_audit_checks.sbatch -- the two CHEAP checks from FORMULA_AUDIT.md, in one
# job: pointwise composition identities (§7.3) and the p97/median shape test
# across all 21 (§7.1).
#
#     sbatch jobs/09_audit_checks.sbatch
#
# Neither needs a single new CDS request. Both run on the 39 months already on
# disk, so this can be submitted while jobs/01_download.sbatch is still
# fetching the production series -- the raw download is diagnostic-independent
# and nothing here can invalidate a byte of it.
#
# WHAT EACH ONE DECIDES
# ---------------------
# check_composition   Verifies nine of the 21 against each other, cell by cell,
#                     with no external reference. Catches a mis-wired argument
#                     or a chunk-boundary error -- faults that produce a
#                     plausible field and survive every distributional test.
#                     Also settles empirically whether the zarr on disk holds
#                     DEF or DEF^2 (FORMULA_AUDIT.md §5).
#
# check_shape_ratios  Turns the two Williams tables from an untransferable
#                     magnitude comparison into a dimensionless shape test that
#                     covers all 21 at once. Predicted in advance to flag f2d;
#                     see FORMULA_AUDIT.md §7.1.
#
# --mem=32G   The shape test loads three DJF months of the GLOBAL derived zarr
#             but subsets to Williams' box (50-75N/10-60W) before loading, so
#             the resident array is small. 32G is generous; do NOT round it up
#             further -- STATUS.md §11.9 measured a 120G request scheduled 22 h
#             out while 24G started in 13 seconds.
# --time=01:00  Both are minutes of real work. Short walltime backfills.
# ===========================================================================

#SBATCH --job-name=cat-audit
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
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
start=$(date +%s)
rc=0

# --- 1. composition identities, on one NA month -----------------------------
# One month is enough: the identities hold pointwise, so more data adds
# confidence about nothing. DJF 1979 is the month every other result in
# STATUS.md is anchored to, which is why it is the one checked.
MONTH=$BASE/derived/north_atlantic/diagnostics_na_1979-01.zarr

echo
echo "############################################################"
echo "# 1/3  POINTWISE COMPOSITION IDENTITIES  (FORMULA_AUDIT 7.3)"
echo "############################################################"
if [ -d "$MONTH" ]; then
    pixi run python ada/check_composition.py "$MONTH"
    rc=$(( rc | $? ))
else
    echo "!! $MONTH not found -- skipping"
    echo "   (run jobs/07_diagnostics_na_trend.sbatch first)"
    rc=1
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
# report it.
echo
echo "############################################################"
echo "# 3/3  SHAPE RATIOS, 3-HOURLY INSTANTANEOUS  (the caveat)"
echo "############################################################"
pixi run python ada/check_shape_ratios.py --base "$BASE" --box williams --no-daily-mean
rc=$(( rc | $? ))

end=$(date +%s)
echo
echo "=== exit ${rc}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
exit $rc
