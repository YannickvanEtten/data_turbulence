#!/bin/bash
# ===========================================================================
# 15_check_production_complete.sbatch -- is the 42-year series actually
# complete, and is it ONE consistent convention throughout, including the
# thresholds it will be compared against?
#
#     sbatch jobs/15_check_production_complete.sbatch
#
# Run this FIRST, before jobs/16_full_trend_check.sbatch or any calibration
# comparison. It reads only zarr .zattrs/.zmetadata JSON sidecars and the
# thresholds file -- no diagnostic data is opened -- so it is cheap and fast
# on purpose: it is meant to be run again after every resubmission of
# jobs/14 or jobs/04, not just once.
#
# --mem=4G    Nothing here touches array data.
# --time=00:10:00   516 small JSON reads.
# ===========================================================================

#SBATCH --job-name=cat-prod-check
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/prodcheck-%j.out
#SBATCH --error=logs/prodcheck-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== production completeness check on $(hostname) at $(date -Is) ==="

start=$(date +%s)
pixi run python ada/check_production_complete.py --base "$BASE" "$@"
status=$?
end=$(date +%s)

echo
echo "=== exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
exit $status
