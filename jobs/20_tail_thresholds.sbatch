#!/bin/bash
# ---------------------------------------------------------------------------
# 20_tail_thresholds.sbatch
#
# Percentile calibration on the FULL global year, via ada/tail_thresholds.py.
# Replaces the compute_thresholds call inside jobs/06 for the full-year case
# only -- jobs/06 still works and still reproduces the 48-day result, which is
# the point of keeping both.
#
# WHY A SEPARATE JOB AND NOT A CHANGE TO jobs/06
#     Because the full derived field is being KEPT on disk (jobs/19), the
#     threshold computation no longer has to be fused into the diagnostics
#     run. Keeping it separate means it can be re-run at a different cut, or
#     validated against the existing 48-day thresholds, without recomputing
#     21 diagnostics over a year.
#
# THE VALIDATION RUN COMES FIRST, AND IT IS FREE
#     Before jobs/17 has downloaded anything, run:
#
#       sbatch --export=ALL,MODE=validate jobs/20_tail_thresholds.sbatch
#
#     That runs the tail path over the EXISTING 48-day stores in
#     derived/global_sub48 and requires it to reproduce
#     calibration/thresholds_2026-09-03.json. It takes minutes and needs no
#     CDS queue. If it fails, the bug is in the tail code and a day of
#     downloading would have been spent on a broken method.
#
# MEMORY
#     Peak is one diagnostic's retained tail, held while it is sorted and
#     tie-consolidated: at a p88 cut on the full year that is ~3.6e8 values,
#     ~3 GB as float64 plus the argsort index and np.unique's buffer, so
#     roughly 12-15 GB transient. 48G is comfortable; if the queue makes that
#     slow, --cut-percentile 94 halves it and still covers Lee's p95.
#
# TIME
#     Two passes over ~255 GB of zarr per diagnostic-ish (one strided, one
#     full), 21 diagnostics. Budget generously; it is I/O bound on a shared
#     filesystem and six concurrent tasks once stretched a 30-min job to 36
#     (STATUS 13.2).
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-tail-calib
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/tailcalib-%j.out
#SBATCH --error=logs/tailcalib-%j.out

set -euo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
MODE=${MODE:-full}

module load 2025
module load pixi

export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== tail_thresholds MODE=$MODE on $(hostname) at $(date -Is) ==="

if [ "$MODE" = "validate" ]; then
    # Reproduce the existing 48-day thresholds through the new code path.
    # No archive is written: this run exists to answer one question.
    pixi run python ada/tail_thresholds.py \
        --base "$BASE" \
        --derived-subdir derived/global_sub48 \
        --cut-percentile 88 \
        --no-archive \
        --validate-against calibration/thresholds_2026-09-03.json
else
    pixi run python ada/tail_thresholds.py \
        --base "$BASE" \
        --derived-subdir derived/global \
        --cut-percentile 88
fi

echo "=== done at $(date -Is) ==="
ls -lh "$BASE/calibration/" | tail -10
