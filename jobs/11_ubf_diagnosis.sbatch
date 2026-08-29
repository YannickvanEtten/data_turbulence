#!/bin/bash
# ===========================================================================
# 11_ubf_diagnosis.sbatch -- why UBF behaves unlike its 20 siblings.
# FORMULA_AUDIT.md §7.2.
#
#     sbatch jobs/11_ubf_diagnosis.sbatch
#
# STATUS.md §12.6 has UBF at +11 % where every sibling ran +38 % to +226 %.
# Prosser's Figure S5 says it should be near the TOP of the 21 over the North
# Atlantic, not the bottom. FORMULA_AUDIT.md §2 has already cleared the
# equation against Sharman A30 term by term, so this looks at the three things
# that remain: how hard the four terms cancel, whether float32 arithmetic is
# reordering the tail, and whether ERA5's own vorticity product disagrees with
# the vorticity derived from u and v by enough to matter.
#
# Two months, both DJF, forty-one years apart -- because the anomaly is in the
# TREND, not the level, and a mechanism that explains it has to look different
# at the two ends of the record.
#
# --mem=32G   Computes the four terms twice (float32 and float64) on a North
#             Atlantic month. float64 doubles the arrays, which is why this is
#             above the 24G that jobs/07 uses.
# --time=01:00  A NA month of all 21 took 5-8 min; this computes four terms.
# ===========================================================================

#SBATCH --job-name=cat-ubf
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ubf-%j.out
#SBATCH --error=logs/ubf-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== UBF diagnosis on $(hostname) at $(date -Is) ==="
start=$(date +%s)
rc=0

for MONTH in 1979-01 2020-01; do
    F=$BASE/raw/north_atlantic/era5_na_${MONTH}.grib
    echo
    echo "############################################################"
    echo "# $MONTH"
    echo "############################################################"
    if [ -f "$F" ]; then
        pixi run python ada/check_ubf.py "$F" --box prosser
        rc=$(( rc | $? ))
    else
        echo "!! $F not found -- skipping"
        rc=$(( rc | 1 ))
    fi
done

end=$(date +%s)
echo
echo "=== exit ${rc}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
echo
echo "Compare the two months. A cancellation ratio and a vorticity discrepancy"
echo "that are similar at both ends of the record cannot by themselves explain"
echo "a flat TREND -- they would depress both years equally. A discrepancy that"
echo "GROWS across the record can, and would also be a candidate for STATUS.md"
echo "§5's risk 3 (time-varying bias from ERA5's changing observing system),"
echo "which is a standing referee objection for the whole project and not only"
echo "for this diagnostic."
exit $rc
