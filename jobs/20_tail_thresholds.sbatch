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
    #
    # Which directory holds the 48-day stores depends on whether the rename
    # ahead of jobs/19 has happened yet, so find them rather than assume.
    # Override with SUB=... if they are somewhere else again.
    SUB=${SUB:-}
    if [ -z "$SUB" ]; then
        for cand in derived/global_sub48 derived/global; do
            if compgen -G "$BASE/$cand/diagnostics_glob_*.zarr" > /dev/null; then
                SUB=$cand; break
            fi
        done
    fi
    if [ -z "$SUB" ]; then
        echo "FATAL: no diagnostics_glob_*.zarr under $BASE/derived/global_sub48" >&2
        echo "       or $BASE/derived/global. Pass SUB=<subdir> explicitly." >&2
        ls -d "$BASE"/derived/*/ 2>/dev/null | sed 's/^/   /' >&2
        exit 1
    fi
    n=$(compgen -G "$BASE/$SUB/diagnostics_glob_*.zarr" | wc -l)
    echo "validating against the $n stores in $SUB"
    if [ "$n" -ne 12 ]; then
        # 12 is the 48-day calibration set thresholds_2026-09-03.json was built
        # from. A different count means this is not that dataset, and the
        # comparison would be meaningless rather than merely wrong.
        echo "FATAL: expected 12 stores, found $n -- refusing to compare against" >&2
        echo "       thresholds_2026-09-03.json, which was built from 12." >&2
        exit 1
    fi
    pixi run python ada/tail_thresholds.py \
        --base "$BASE" \
        --derived-subdir "$SUB" \
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
