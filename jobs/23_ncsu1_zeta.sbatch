#!/bin/bash
# ---------------------------------------------------------------------------
# 23_ncsu1_zeta.sbatch
#
# The last cheap shot at the biggest gap in the replication, before committing
# seven hours to a re-derive.
#
# WHAT IS BEING ASKED
#     #19 `ncsu1` sits at a 1979 exceedance level ratio of 6.07 against Prosser
#     -- six times too HIGH, and the only one of the 21 whose error has the
#     opposite sign to every other. AUDIT §5.6 listed three candidate causes:
#
#       1. the MAX(Ri, 1e-5) floor        -- MEASURED 2026-09-09, DEAD.
#          Job 1102936 step 7: the floor binds on 0.0356 % of global cells and
#          reaches 0.58 % of NCSU1's p97 tail. It cannot carry a factor of six.
#       2. |grad zeta| from ERA5's ARCHIVED vorticity   <-- THIS JOB
#       3. the MAX(., 0) clip on the advection term.
#
#     Candidate 2: NCSU1 is the only one of the 21 containing |grad zeta|
#     multiplied by two other flow factors, i.e. effectively a SECOND spatial
#     derivative of the wind. Ours uses ERA5's spectrally-derived vorticity,
#     which carries small-scale power a 0.25 deg centred difference cannot
#     represent, concentrated in the midlatitude storm tracks rather than in
#     the tropics that dominate the global threshold. Prosser (2023) p. 2
#     downloaded only u, v, T and z, so his zeta is a finite difference of the
#     wind and is damped by comparison. That asymmetry predicts the right SIGN
#     for a 6.07 -- which candidate 1 never did.
#
# WHY GLOBAL, AND WHY THE LATITUDE TABLE IS THE ACTUAL RESULT
#     The hypothesis is about region-versus-globe CONTRAST, not about
#     magnitude. A flip rate measured inside the North Atlantic box cannot see
#     that. On a global file, ab_compare prints exceedance frequency per
#     latitude band for both versions, and the prediction is specific:
#
#       CONFIRMED   midlatitude bands (+30..+60, -60..-30) LOSE exceedance
#                   frequency relative to the tropics. That is the 6.07
#                   mechanism, and it would mean #19's disagreement is the
#                   archived-vs-computed axis (AUDIT §3.4, §5.5) rather than
#                   anything about A36.
#       FALSIFIED   a large flip with a FLAT latitude profile. The
#                   substitution matters, but not for this reason, and
#                   candidate 3 (the clip) is what is left.
#       DEAD        flip < 0.5 % everywhere. Candidate 2 joins candidate 1 and
#                   #19 stays unexplained -- which is itself worth knowing
#                   before the re-derive rather than after it.
#
# NOTE WHAT THIS IS NOT. A36's zeta is not specified as archived or computed;
# Sharman (2006) p. 284 and Williams & Storer (2022) Eq. (7) p. 1428 both just
# write |grad zeta|. So neither reading is wrong, and this is a provenance
# choice to be made and cited, exactly like F3 and F10 -- not a bug fix.
#
# SIZING. Same shape as jobs/22 step 2: one global month capped to 32
# timesteps, two full result sets. 48G fits every node in defq (STATUS §11.9);
# 2 h against a measured 57 min for the equivalent step.
#
#   mkdir -p logs
#   sbatch jobs/23_ncsu1_zeta.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-ncsu1
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ncsu1-%j.out
#SBATCH --error=logs/ncsu1-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
OUT=data_prosser/ab
GLOBAL_STEPS=${GLOBAL_STEPS:-32}
GLOBAL_FILE=${GLOBAL_FILE:-$BASE/raw/global/era5_glob_2000-01.grib}

module load 2025
module load pixi

export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }
mkdir -p "$OUT"

if [ ! -f "$GLOBAL_FILE" ]; then
    echo "!! GLOBAL_FILE not found: $GLOBAL_FILE"
    exit 2
fi

echo "=== NCSU1 vorticity provenance (audit §5.6 candidate 2) on $(hostname) at $(date -Is) ==="
echo "  global file $GLOBAL_FILE  (first $GLOBAL_STEPS timesteps)"
echo

start=$(date +%s)
pixi run python ada/ab_compare.py "$GLOBAL_FILE" \
    --max-timesteps "$GLOBAL_STEPS" \
    --fix ncsu1_computed_vorticity \
    --out "$OUT"
status=$?
end=$(date +%s)

echo
echo "=== exit=$status in $(( (end - start) / 60 )) min at $(date -Is) ==="
echo "READ THE LATITUDE TILT TABLE, not just the flip rate. Midlatitude bands"
echo "losing exceedance frequency relative to the tropics confirms candidate 2;"
echo "a flat profile falsifies it and leaves candidate 3 (the MAX(.,0) clip)."
exit $status
