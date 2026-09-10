#!/bin/bash
# ---------------------------------------------------------------------------
# 26_diagnostics_global_audit.sbatch
#
# Step 2 of the re-derive: the year-2000 global reference year, recomputed
# under the SAME audit fix set as jobs/25, written to derived/global_audit/.
#
# WHY THIS IS NOT OPTIONAL, AND WHY IT MUST FINISH BEFORE ANY COMPARISON
#     Severity thresholds are cos(phi)-weighted percentiles of the reference
#     year. Applying thresholds computed under one convention to diagnostics
#     computed under another is not a small error -- it is meaningless. That
#     is exactly what jobs/15 caught on 2026-09-03 (STATUS 14.1 pt 3), when
#     the North Atlantic side had moved to the corrected deformation/f2d
#     conventions and the global side had not. For a rank-preserving change it
#     would merely be wrong; for f2d, whose A->C variant is not monotone
#     (rho(A,D) = -0.94), it is nonsense.
#
#     F1 changes 13 of 21 and F12 changes 10 of 21 (STATUS 18.2, 18.12), so
#     the thresholds MUST be rebuilt. jobs/25's 504 audit stores are not
#     comparable to anything until this has run and jobs/20 has been re-run
#     against its output.
#
# THE ORDER, WRITTEN DOWN SO IT IS NOT GUESSED AT
#     1. jobs/25   NA diagnostics, audit fixes      -> derived/north_atlantic_audit/   ~7 h
#     2. jobs/26   global reference year, same fixes -> derived/global_audit/          ~14 h   <-- THIS
#     3. jobs/20   tail thresholds on THAT directory -> calibration/thresholds_<date>.json  ~2 h
#                  sbatch --export=ALL,DERIVED=derived/global_audit jobs/20_tail_thresholds.sbatch
#     4. jobs/16   the 42-year trend, on the new thresholds and the audit stores  ~6 h
#     5. jobs/21   per-diagnostic trends                                          ~1.5 h
#
#     1 and 2 are independent and could overlap, but the QOS caps concurrency
#     at 8 jobs total (STATUS 11.9) so they would only contend -- and STATUS
#     15.11 measured a 54 % slowdown from concurrent readers on the same
#     stores. Run them in sequence.
#
# SIZING -- ALL MEASURED ON jobs/19, WHICH THIS MIRRORS
#     Peak RSS 42.7 GB at --chunk-days 1 (job 1098825). 64G, unchanged: F12's
#     three substituted fields are computed per chunk, which at 8 timesteps on
#     the 721x1440 grid is ~0.2 GB each -- noise against 42.7.
#
#     WALLTIME 12 h, NOT 8, and this is the expensive lesson in STATUS 15.7
#     item 2. Job 1098825 measured 3:18-5:01 per task with later waves slower
#     as derived/global grew, and task 4 (May) blew through an 8 h cap and
#     wrote NOTHING after burning 8 hours. A TIMEOUT is the worst failure
#     available: full cost, zero output.
#
#     %4 rather than %8, matching jobs/19, because that is the concurrency the
#     3:18-5:01 timings were measured at. Three waves, ~14 h total. This is a
#     DAY run, not an overnight one.
#
# DISK. ~196 GB, beside the existing 196 GB baseline year, which stays.
# Together with jobs/25's ~294 GB and jobs/20's ~41 GB tail archive that is
# ~531 GB on top of 2026-09-08's 985 GB, i.e. ~1.52 TB of 2.5 TB, 61 %.
# Comfortable, and STATUS 16.2's headroom note anticipated exactly this.
#
#   mkdir -p logs
#   sbatch jobs/26_diagnostics_global_audit.sbatch     # BARE -- see STATUS 11.9
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-glob-audit
#SBATCH --partition=defq
#SBATCH --array=0-11%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/diag-glob-audit-%A_%a.out
#SBATCH --error=logs/diag-glob-audit-%A_%a.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
YEAR=2000
OUTDIR=$BASE/derived/global_audit

# MUST match jobs/25 exactly. If these two ever diverge, the thresholds and the
# data they are applied to are built under different conventions and every
# downstream number is meaningless -- the 2026-09-03 trap again.
FIXES=${FIXES:---fix meridional_metric --fix four_variable_provenance}

module load 2025
module load pixi

export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

MONTH=$(printf "%02d" $(( SLURM_ARRAY_TASK_ID + 1 )))
IN="$BASE/raw/global/era5_glob_${YEAR}-${MONTH}.grib"
OUT="$OUTDIR/diagnostics_glob_${YEAR}-${MONTH}.zarr"

echo "=== global audit re-derive ${YEAR}-${MONTH} on $(hostname) at $(date -Is) ==="
echo "    fixes: $FIXES"
if [ ! -f "$IN" ]; then
    echo "!! input not found: $IN"
    echo "!! these are the MERGED contiguous months from jobs/18, not the"
    echo "!! 10-day day-blocks. Check raw/global/ holds 12 era5_glob_2000-MM.grib."
    ls -lh "$BASE/raw/global/" 2>/dev/null | tail -20 | sed 's/^/     /'
    exit 1
fi
ls -lh "$IN" | sed 's/^/     /'

mkdir -p "$OUTDIR"

start=$(date +%s)
pixi run python ada/diagnostics_global.py "$IN" "$OUT" \
    --chunk-days 1 $FIXES --skip-if-matching
status=$?
end=$(date +%s)

echo "=== ${YEAR}-${MONTH} exit=$status in $(( (end - start) / 60 )) min at $(date -Is) ==="
du -sh "$OUT" 2>/dev/null | sed 's/^/     /'
df -h "$BASE" | tail -1
exit $status
