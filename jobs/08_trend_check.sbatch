#!/bin/bash
# ===========================================================================
# 08_trend_check.sbatch -- fit the five DJF seasons and compare the CHANGE
# against Prosser (2023) Table 1.
#
#     sbatch jobs/08_trend_check.sbatch
#     sbatch jobs/08_trend_check.sbatch --per-diagnostic
#
# The calibration check verified a LEVEL. This verifies the TREND, which is
# the paper's actual claim and the half of the pipeline (Layers 5-6) that has
# never seen real data.
#
# --per-diagnostic also breaks DJF 1979 and DJF 2020 down by all 21, which is
# the check the 21-member ensemble mean cannot do: a badly wrong diagnostic
# moves that mean by 1/21 and hides. Costs a few minutes; worth it.
#
# --mem=16G   Small. The Prosser box is 97 x 181 gridpoints; a whole DJF is
#             ~1.6e7 points per diagnostic, an order less than one global
#             calibration month. No sorting either -- the thresholds are read
#             from the JSON, not recomputed.
# --time=01:00  Minutes of real work. Short walltime, easy backfill.
# ===========================================================================

#SBATCH --job-name=cat-trend
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/trend-%j.out
#SBATCH --error=logs/trend-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== trend check on $(hostname) at $(date -Is) ==="

start=$(date +%s)
pixi run python ada/trend_check.py --base "$BASE" "$@"
status=$?
end=$(date +%s)

echo
echo "=== exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
exit $status
