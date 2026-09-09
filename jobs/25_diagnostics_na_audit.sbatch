#!/bin/bash
# ---------------------------------------------------------------------------
# 25_diagnostics_na_audit.sbatch
#
# THE RE-DERIVE. 504 North Atlantic months, 1979-2020, recomputed under the
# audit fix set, written BESIDE the existing production series rather than
# over it.
#
# WHAT CHANGES, AND WHY THESE AND NOT THE OTHERS
#     --fix meridional_metric        F1. dy is the plate-carree map distance
#                                    a*dphi, so the meridional map-scale factor
#                                    is applied exactly ONCE. Measured 4.2e-3
#                                    -> 2e-7 against exact ellipsoidal
#                                    gradients; moves 13 of 21 (STATUS 18.2,
#                                    18.5). NOTE this is a deliberate departure
#                                    from MetPy, whose nominal lat/lon grid
#                                    deltas mix the two coordinates.
#     --fix four_variable_provenance F12. zeta, delta and PV rebuilt from
#                                    u, v, T, z -- Prosser's four fields --
#                                    instead of ERA5's archive. Moves 10 of 21,
#                                    0 inert (STATUS 18.12). SUPERSEDES F2, F3
#                                    and F11, which are retired as separate
#                                    fixes and must NOT also be passed.
#
#     F4 (deformation un-squared), F6 (f2d variant A) and F8 (the grid
#     assertions) need no flag -- they are the code defaults now, and the
#     store records f2d_variant and deformation_convention either way.
#
#     F10 (Endlich closed-form shear) is DELIBERATELY EXCLUDED. Both readings
#     of "centred second-order finite differences" are defensible (AUDIT 4.14)
#     and neither paper resolves it, so it is reported as a sensitivity, not
#     baked into the production series.
#
# THE OUTPUT DIRECTORY IS DIFFERENT ON PURPOSE
#     derived/north_atlantic_audit/, not derived/north_atlantic/. The existing
#     504 stores are the published-convention series behind every number in
#     STATUS 15-17, and they stay. This is a side-by-side, like jobs/24 was at
#     the diagnostic level. Cost: ~294 GB against 1.5 TB free (STATUS 16.2).
#
# THE SKIP GUARD NOW CHECKS THE FIX SET -- THIS IS NEW AND IT MATTERS
#     `existing_output_matches()` compared f2d_variant and
#     deformation_convention only. Every store on disk was written under
#     BASELINE_FIXES with no `audit_fixes` attribute at all, so a resubmission
#     with --fix would have SKIPPED all 504 as "already current" and reported
#     success -- the 2026-09-03 stale-convention trap (STATUS 14.1 pt 3) a
#     second time, one level deeper. Fixed 2026-09-09: the attribute is written
#     and compared, and absence FAILS CLOSED (absent == baseline). That is what
#     makes --skip-if-matching safe to use here, and it is why this job can be
#     resubmitted freely after any interruption.
#
# SIZING -- MEASURED, NOT ROUNDED UP
#     A North Atlantic month is 121 x 301 x 248 x 3. Production ran at 24G with
#     a measured ~15 GB peak (STATUS 10.12, 13.3). F12 adds three substituted
#     fields at that shape, ~0.7 GB. 32G is the measured peak plus real
#     headroom, and STATUS 11.9 is unambiguous that rounding memory up costs a
#     day of queue: a 120 GB request scheduled 22 hours out where 24 GB started
#     in 13 seconds.
#
#     Walltime 2 h against a measured 5-8.5 min per month, because six
#     concurrent tasks stretched a 30-minute job to 36 on shared filesystem I/O
#     (STATUS 13.2) and STATUS 15.11 measured a 54 % slowdown from concurrent
#     readers alone. Short honest walltimes get backfilled.
#
#     %8 is the QOS cap (MaxJobsPU = 8, STATUS 11.9) -- any higher throttle is
#     a no-op. No CDS here, so the %4 download limit does not apply.
#
# EXPECTED WALL CLOCK: ~7 h for 504 months at %8 (STATUS 13.2).
#
#   mkdir -p logs                 # SLURM opens the log BEFORE the script runs
#   sbatch jobs/25_diagnostics_na_audit.sbatch      # BARE -- see below
#
# DO NOT PASS --array ON THE COMMAND LINE. It replaces the ENTIRE spec, %8
# included, and fails silently (STATUS 11.9). If you must subset, re-append it:
#   sbatch --array=0-11%8 jobs/25_diagnostics_na_audit.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-diag-audit
#SBATCH --partition=defq
#SBATCH --array=0-503%8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/diag-audit-%A_%a.out
#SBATCH --error=logs/diag-audit-%A_%a.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
OUTDIR=$BASE/derived/north_atlantic_audit

# The fix set, in one place. Override for a sensitivity run, e.g.
#   sbatch --export=ALL,FIXES="--fix meridional_metric" jobs/25_...
FIXES=${FIXES:---fix meridional_metric --fix four_variable_provenance}

module load 2025
module load pixi

export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

# index 0 = 1979-01, index 503 = 2020-12
YEAR=$(( 1979 + SLURM_ARRAY_TASK_ID / 12 ))
MONTH=$(printf "%02d" $(( SLURM_ARRAY_TASK_ID % 12 + 1 )))

IN="$BASE/raw/north_atlantic/era5_na_${YEAR}-${MONTH}.grib"
OUT="$OUTDIR/diagnostics_na_${YEAR}-${MONTH}.zarr"

echo "=== NA audit re-derive ${YEAR}-${MONTH} on $(hostname) at $(date -Is) ==="
echo "    fixes: $FIXES"
if [ ! -f "$IN" ]; then
    echo "!! input not found: $IN"
    echo "!! all 504 raw months were confirmed present on 2026-09-04 (STATUS 14.7),"
    echo "!! so this is a path or index error, not a missing download."
    exit 1
fi

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
