#!/bin/bash
# ---------------------------------------------------------------------------
# 27_final_checks.sbatch                      FINAL-CHECKS-2026-09-16
#
# THE ONE BATCH that closes the data-gathering phase (FINAL_DATASET.md §4).
# Read-only against the data. Writes only into cat_outputs/final/.
#
#   step 1  ada/final_sweep.py                  ~5-15 min
#           all 504 derived/north_atlantic_audit stores + 12 derived/global_audit
#           stores: attributes, variables, grid, time axis, chunk census, one
#           read per variable; thresholds_2026-09-11.json ladders; the job log
#           that last wrote each store.
#   step 1b ada/final_tail_latitude.py           ~5-10 min
#           zonal profile of every diagnostic's LOG/MOG/SOG exceedance set in
#           the reference year, read from calibration/tails_2026-09-11 (the
#           archive stores each value's latitude row). Tests the one spatial
#           statement the literature makes per diagnostic (W&S 2022 p. 1428).
#   step 2  ada/per_diagnostic_all_severities.py   ~40-75 min (estimate)
#           annual per-diagnostic fits at all five severities in one pass over
#           the 42 years, plus the finite census of every value in the box.
#           Checkpoints every year: if this job times out, resubmitting the
#           same command resumes where it stopped.
#   step 3  ada/final_scorecard.py               seconds
#           joins step 2 with Prosser Figure 4, S4-LOG and S4-SOG; checks the
#           new MOG table against the existing audit MOG table (job 1108148)
#           and the five convention-invariant rows against the baseline.
#
# BUDGET: hard walltime 2 h 45 min (estimate 1-1.7 h), under the 3 h cap. Nothing else should be
# reading derived/north_atlantic_audit while this runs (STATUS 15.11: three
# concurrent readers cost 54 %).
#
# 24G: one year of the Prosser box in memory is 4.3 GB of float32 (21 x 97 x
# 181 x 2920 x 4 B) plus ~2 GB of work arrays. The driver prints its own peak.
#
# BEFORE SUBMITTING (STATUS 18.13: verify the file on ADA, not the commit):
#   cd $BASE/data_turbulence && git pull
#   grep -l FINAL-CHECKS-2026-09-16 jobs/27_final_checks.sbatch ada/final_sweep.py \
#        ada/per_diagnostic_all_severities.py ada/final_scorecard.py \
#        ada/prosser_published_s4.py ada/final_tail_latitude.py
#   -> must list all six files
#   mkdir -p logs && sbatch jobs/27_final_checks.sbatch
#
# AT 60 SECONDS:  head -40 logs/final-<jobid>.out
#   -> step 0 prints the git revision and the six sha256 sums; step 1 prints
#      the expected audit_fixes and starts listing stores.
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-final
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:45:00
#SBATCH --output=logs/final-%j.out
#SBATCH --error=logs/final-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
SERIES=derived/north_atlantic_audit
THRESHOLDS=$BASE/calibration/thresholds_2026-09-11.json

module load 2025
module load pixi
export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }
mkdir -p cat_outputs/final

echo "=== step 0: configuration on $(hostname) at $(date -Is) ==="
echo "    git     $(git rev-parse --short HEAD) $(git log -1 --format=%cd --date=iso)"
echo "    series  $BASE/$SERIES"
echo "    thresh  $THRESHOLDS"
for f in jobs/27_final_checks.sbatch ada/final_sweep.py ada/per_diagnostic_all_severities.py \
         ada/final_scorecard.py ada/prosser_published_s4.py ada/final_tail_latitude.py; do
    if ! grep -q FINAL-CHECKS-2026-09-16 "$f" 2>/dev/null; then
        echo "!! $f is missing or is not the 2026-09-16 version -- git pull first"; exit 3
    fi
    echo "    $(sha256sum "$f")"
done
if [ ! -f cat_outputs/per_diagnostic_annual_moderate_audit.csv ]; then
    echo "   (note) cat_outputs/per_diagnostic_annual_moderate_audit.csv absent:"
    echo "          step 3's regression check will report 'not available'"
fi

echo; echo "=== step 1: store / thresholds / log sweep  $(date -Is) ==="
s1=$(date +%s)
pixi run python ada/final_sweep.py --base "$BASE" --na-subdir "$SERIES" \
    --global-subdir derived/global_audit --thresholds "$THRESHOLDS"
st1=$?
echo "=== step 1 exit=$st1 in $(( ($(date +%s) - s1) / 60 )) min ==="

echo; echo "=== step 1b: tail latitude profile, reference year  $(date -Is) ==="
s1b=$(date +%s)
pixi run python ada/final_tail_latitude.py --base "$BASE" \
    --archive calibration/tails_2026-09-11 --thresholds "$THRESHOLDS"
st1b=$?
echo "=== step 1b exit=$st1b in $(( ($(date +%s) - s1b) / 60 )) min ==="

echo; echo "=== step 2: per-diagnostic fits, all severities  $(date -Is) ==="
s2=$(date +%s)
pixi run python ada/per_diagnostic_all_severities.py --base "$BASE" \
    --derived-subdir "$SERIES" --thresholds "$THRESHOLDS" \
    --outdir cat_outputs/final
st2=$?
echo "=== step 2 exit=$st2 in $(( ($(date +%s) - s2) / 60 )) min ==="

echo; echo "=== step 3: scorecard  $(date -Is) ==="
pixi run python ada/final_scorecard.py
st3=$?
echo "=== step 3 exit=$st3 ==="

echo; echo "=== done at $(date -Is): sweep=$st1 tails=$st1b fits=$st2 scorecard=$st3 ==="
echo "    bring back: git add cat_outputs/final/*.csv cat_outputs/final/*.md logs/final-${SLURM_JOB_ID}.out"
exit $(( st1 | st1b | st2 | st3 ))
