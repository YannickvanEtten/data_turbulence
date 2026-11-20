#!/bin/bash
# ---------------------------------------------------------------------------
# 28_ubf_stencil.sbatch                                (UBF-STENCIL-2026-11-20)
#
# Can `ubf` be made usable as a stand-alone indicator? One test, three
# discretisations of the SAME formula (Sharman A30), side by side. Only ubf is
# computed; nothing on the share is written except cat_outputs/ubf_stencil/.
# Full rationale, versions and pre-registered predictions: the docstring of
# ada/ubf_stencil_test.py. In one paragraph:
#
#   grad^2 Phi is two first differences, f*zeta is one, so even an exactly
#   balanced wave leaves a numerical residual (36 % of f*zeta at 4 grid steps,
#   10 % at 8). If that mismatch is regionally concentrated it can take over a
#   GLOBAL percentile tail -- which would explain ubf's polar reference-year
#   tail and its NA level ratio of 0.39 / 0.19 / 0.08 against Prosser.
#   V0 = current, V1 = compact Laplacian, V2 = 4th-order differences
#   (mismatch 15 % / 1.2 % at 4 / 8 grid steps).
#
# WHAT IS ALREADY KNOWN, so the log is read with the right prior
#   - V0 is bit-identical to essentials.ubf (checked locally) and the script
#     checks it against the production store before anything else (stage 1).
#   - On the 2-step validation subset (July 2016, NA) the residual is ~60 % of
#     its largest term under ALL THREE versions, i.e. the mismatch is not the
#     bulk of ubf. The tail still moves (~15 % of the top cells change). So the
#     expectation going in is MODEST, and P4 (falsification) is a live outcome.
#
# READ IN THIS ORDER
#   1. Stage 1: Spearman >= 0.99999. If not, the job exits 2 -- stop there.
#   2. Stage 2's "box density, V / V0" line and the 60S-90S / 60N-90N columns
#      against P3 and P4 (printed right under the table).
#   3. Stage 3's implied fitted level ratios against the class-rule band
#      0.75-1.33.
#
# NO NEW DATA. Inputs, all already on the share:
#   raw/global_sub48/*.grib              12 files, year 2000, days 1/9/17/25
#   raw/north_atlantic/era5_na_1979-*.grib  12 files
#   derived/global_audit/diagnostics_glob_2000-01.zarr   (control only)
#
# SIZING. Three ubf fields of 721 x 1440 x 384 float32 (~4.8 GB) plus the
# exceedance masks: ~12-15 GB measured on the synthetic run scaled up. 32G
# fits every defq node. Walltime 4 h is generous on purpose (STATUS 15.7 #2:
# a TIMEOUT costs everything) -- expected run time well under 2 h, mostly GRIB
# decode. Peak RSS is printed after every file.
#
#   cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
#   git pull && grep -l UBF-STENCIL-2026-11-20 ada/ubf_stencil_test.py jobs/28_ubf_stencil.sbatch
#   mkdir -p logs && sbatch jobs/28_ubf_stencil.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-ubf-stencil
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/ubfstencil-%j.out
#SBATCH --error=logs/ubfstencil-%j.out

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
if ! grep -q "UBF-STENCIL-2026-11-20" ada/ubf_stencil_test.py; then
    echo "FATAL: ada/ubf_stencil_test.py is not the 2026-11-20 version -- was it pushed and pulled?"
    exit 2
fi

echo "=== ubf stencil test on $(hostname) at $(date -Is) ==="
echo "    git rev $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
sha256sum ada/ubf_stencil_test.py jobs/28_ubf_stencil.sbatch
echo "    global sub-sample files: $(ls $BASE/raw/global_sub48/*.grib 2>/dev/null | wc -l) (need 12)"
echo "    NA 1979 files:           $(ls $BASE/raw/north_atlantic/era5_na_1979-*.grib 2>/dev/null | wc -l) (need 12)"
echo

s=$(date +%s)
pixi run python ada/ubf_stencil_test.py --base "$BASE" --out cat_outputs/ubf_stencil
st=$?
echo
echo "=== exit $st after $(( ($(date +%s) - s) / 60 )) min at $(date -Is) ==="
echo "Outputs: $REPO/cat_outputs/ubf_stencil/{summary.json,global_tail.csv}"
echo "Bring back: the log, summary.json and global_tail.csv (git add -f if ignored)."
exit $st
