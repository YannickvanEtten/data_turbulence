#!/bin/bash
# THE CHECK TO RUN FIRST. Uses the data already on the share; overwrites nothing.
#     sbatch essentials/jobs/0_check_against_original.sbatch
#
#  1. unit tests of the building blocks                                     (~10 s)
#  2. both conventions, one real month (1979-01), essentials vs original code
#     in memory, bit for bit                                                (~40 min)
#  3. essentials writes 1979-01 under both conventions to $BASE/essentials/
#     and compares with the stores the original pipeline wrote             (~20 min)
# Expected: "ALL 21 BIT-IDENTICAL" twice (x2 with the F10 variant) and "ALL 21 IDENTICAL" twice.
#SBATCH --job-name=ess-check
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=logs/ess-check-%j.out
set -uo pipefail
source essentials/jobs/_common.sh
M=1979-01
pixi run python -m pytest essentials/tests/test_operators.py -q
pixi run python essentials/tests/compare_with_original.py "$BASE/raw/north_atlantic/era5_na_$M.grib"
for CONV in audit baseline; do
    run diagnostics --domain north_atlantic --year 1979 --month 1 --convention $CONV --no-skip
done
run compare-store "$BASE/essentials/derived/north_atlantic_audit/diagnostics_na_$M.zarr" \
                  "$BASE/derived/north_atlantic_audit/diagnostics_na_$M.zarr"
run compare-store "$BASE/essentials/derived/north_atlantic_baseline/diagnostics_na_$M.zarr" \
                  "$BASE/derived/north_atlantic/diagnostics_na_$M.zarr"
