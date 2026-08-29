#!/bin/bash
# ===========================================================================
# 11_ubf_diagnosis.sbatch -- why UBF behaves unlike its 20 siblings.
# FORMULA_AUDIT.md §7.2.
#
#     mkdir -p logs && sbatch jobs/11_ubf_diagnosis.sbatch
#
# STATUS.md §12.6 has UBF at +11 % where every sibling ran +38 % to +226 %.
# Prosser's Figure S5 says it should be near the TOP of the 21 over the North
# Atlantic, not the bottom. FORMULA_AUDIT.md §2 has already cleared the
# equation against Sharman A30 term by term, so this looks at the three things
# that remain: how hard the four terms cancel, whether float32 arithmetic is
# reordering the tail, and whether ERA5's own vorticity product disagrees with
# the vorticity derived from u and v by enough to matter.
#
# Two months, both January, forty-one years apart -- because the anomaly is in
# the TREND, not the level, and a mechanism that explains it has to look
# different at the two ends of the record.
#
# SAFE TO RUN DURING THE DOWNLOAD. This reads raw/north_atlantic/*.grib while
# jobs/01_download.sbatch may be writing into the same directory. That is fine:
# the downloader writes <name>.tmp and atomically renames (STATUS.md §10.7), so
# a complete file is never seen half-written, and months already on disk are
# skipped rather than rewritten. cfgrib will drop a .idx sidecar next to each
# file it opens, which is expected and harmless (STATUS.md §13.1).
#
# --mem=16G   One NA month is 121 x 301 x 248 = 9.0e6 cells. This computes four
#             terms twice, once in float32 and once in float64, so roughly
#             twenty arrays of ~72 MB plus the decoded GRIB -- a few GB. 16 GB
#             schedules immediately; STATUS.md §11.9 measured what asking for
#             more costs.
# --time=01:00  A NA month of all 21 took 5-8 min; this computes four terms.
# ===========================================================================

#SBATCH --job-name=cat-ubf
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
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
echo "    repo $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
start=$(date +%s)
rc=0
ran=0

# The two ends of the trend check's record (STATUS.md §12.1). Both are among
# the 27 NA months already on disk; the script reports which it actually used
# rather than assuming, because the production download is adding months to
# this directory as it runs.
for MONTH in 1979-01 2020-01; do
    F=$BASE/raw/north_atlantic/era5_na_${MONTH}.grib
    echo
    echo "############################################################"
    echo "# $MONTH"
    echo "############################################################"
    if [ -f "$F" ]; then
        echo "using $F"
        pixi run python ada/check_ubf.py "$F" --box prosser
        rc=$(( rc | $? ))
        ran=$(( ran + 1 ))
    else
        echo "!! $F not found -- skipping"
    fi
done

if [ "$ran" -eq 0 ]; then
    echo
    echo "!! neither month was available. Available NA months:"
    ls -1 "$BASE"/raw/north_atlantic/era5_na_*.grib 2>/dev/null | head -20
    rc=1
fi

end=$(date +%s)
echo
echo "=== exit ${rc}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="

if [ "$ran" -eq 2 ]; then
    echo
    echo "Compare the two months. A cancellation ratio and a vorticity"
    echo "discrepancy that are SIMILAR at both ends of the record cannot by"
    echo "themselves explain a flat TREND -- they would depress both years"
    echo "equally. A discrepancy that GROWS across the record can, and would"
    echo "also be a candidate for STATUS.md 5's risk 3 (time-varying bias as"
    echo "ERA5's observing system changes), which is a standing referee"
    echo "objection for the whole project and not only for this diagnostic."
fi
exit $rc
