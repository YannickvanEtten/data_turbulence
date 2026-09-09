#!/bin/bash
# ---------------------------------------------------------------------------
# 24_four_variable.sbatch
#
# F12. THE SIDE-BY-SIDE, not a replacement: compute all 21 twice on the same
# input, once with ERA5's archived `vo`, `d` and `pv` and once with all three
# rebuilt from u, v, T and z, and report which cells change their exceedance
# status. Nothing on disk is overwritten and no production store is touched.
#
# WHY THIS JOB EXISTS, AND WHY IT REPLACES THREE EARLIER ONES
#     Prosser (2023) p. 2 downloads FOUR fields -- "zonal and meridional wind
#     speed, dry bulb temperature, and geopotential height" -- and then states
#     "The 21 turbulence diagnostics were then calculated from the extracted
#     reanalysis fields." Every vorticity, divergence and potential vorticity
#     in his 21 is therefore a finite difference of that wind on his own grid.
#     This project downloads SEVEN and takes ERA5's archived vo/d/pv because
#     CDS offers them.
#
#     That single convenience is the whole archived-versus-computed axis, and
#     it was chased three separate times -- F2 (ubf), F3 (magnitude_pv), F11
#     (ncsu1) -- as if it were three questions. It is one. F12 substitutes all
#     three fields on the DATASET, before any diagnostic runs, so it also
#     reaches the consumers no per-diagnostic flag could: rojak's
#     horizontal_divergence and ti2 (archived `d`), and brown1,
#     vorticity_squared, rva_magnitude and nva (archived `vo`).
#
# THIS IS NOT A CORRECTION, AND THE LOG SHOULD NOT BE READ AS ONE
#     ERA5's archived zeta and delta are the IFS's own PROGNOSTIC spectral
#     variables -- IFS Documentation CY49R1, Part III, section 2.2.7: "The
#     model state at time t is defined by the spectral coefficients of zeta,
#     D, T, q and ln ps" -- and u and v are derived FROM them by inverse
#     spectral transform. Recomputing zeta from archived u, v on a 0.25 deg
#     grid is a lossy round trip, not a refinement. The archived field is
#     closer to the model; the computed field is closer to Prosser.
#
#     Potential vorticity runs the same way and further: ERA5's `pv` is
#     diagnosed by ECMWF on the model's L137 vertical grid, so its dtheta/dp
#     is better resolved than anything obtainable from a 50 hPa stencil -- or
#     from Prosser's 18 hPa one. Substituting Sharman A18 is WORSE physics and
#     CLOSER to the replication target. Both are legitimate; the paper has to
#     say which it used and why, and this job is what turns that into a
#     measured choice rather than an inherited one.
#
# WHAT TO READ, IN ORDER
#     1. The MATERIAL/inert line and the per-diagnostic flip table. Expect a
#        LOT to move -- this touches nine or ten of the 21, where F2/F3/F11
#        touched one each. Anything that comes back bit-identical is a
#        diagnostic that never read an archived field, and that list is worth
#        recording in its own right.
#     2. THE LATITUDE TILT TABLE, which is the actual finding. F11 measured
#        this substitution for ncsu1 alone and the tilt was TROPICS-NEGATIVE
#        (-0.106/-0.136/-0.122) and EXTRATROPICS-POSITIVE (+0.02..+0.306) --
#        the opposite sign to the storm-track mechanism that was predicted.
#        The prediction on record here is that the same sign appears across
#        the whole archived-field group, because the excess small-scale power
#        in ERA5's spectral zeta relative to a 0.25 deg difference is
#        TROPICAL (convection), not midlatitude. If the tilt is flat, or
#        reverses, that reading is wrong and should be dropped rather than
#        rescued.
#     3. Step 2's NA numbers, which are the ones that feed the replication.
#
# NO NEW DATA. Both inputs are already on the share. No CDS.
#
# SIZING. Same shape as jobs/22: two full result sets in memory at once, so
# ~2x a normal pass. 48G fits every node in defq (STATUS 11.9) -- do NOT round
# up, a 120G request was scheduled 22 hours out where 24G started in 13 s.
# The global step is capped at 32 timesteps for the same reason jobs/22 caps
# it: a full month is 248 steps and ~43 GB of diagnostic arrays.
#
#   mkdir -p logs
#   sbatch jobs/24_four_variable.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-4var
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/fourvar-%j.out
#SBATCH --error=logs/fourvar-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
OUT=data_prosser/ab
GLOBAL_STEPS=${GLOBAL_STEPS:-32}       # >= 3 required by f2d's centred d/dt
GLOBAL_FILE=${GLOBAL_FILE:-$BASE/raw/global/era5_glob_2000-01.grib}
NA_FILE=${NA_FILE:-$BASE/raw/north_atlantic/era5_na_2000-01.grib}

module load 2025
module load pixi

export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }
mkdir -p "$OUT"

for v in GLOBAL_FILE NA_FILE; do
    f=${!v}
    if [ ! -f "$f" ]; then
        echo "!! $v not found (got '$f')"
        exit 2
    fi
done

echo "=== F12 four-variable provenance side-by-side on $(hostname) at $(date -Is) ==="
echo "  global file $GLOBAL_FILE  (first $GLOBAL_STEPS timesteps)"
echo "  NA file     $NA_FILE"
echo "  output      $REPO/$OUT"
echo

declare -a NAMES STATUSES MINUTES

step () {
    local label="$1"; shift
    echo
    echo "############################################################"
    echo "### $label"
    echo "############################################################"
    local s=$(date +%s)
    "$@"
    local st=$?
    local e=$(date +%s)
    NAMES+=("$label"); STATUSES+=("$st"); MINUTES+=("$(( (e - s) / 60 ))")
    echo "--- $label: exit=$st in $(( (e - s) / 60 )) min ---"
    return 0
}

# 1. The suite first. F12 added a dataset-level substitution that every
#    diagnostic sees; if it broke a pinned analytic value, every number below
#    is meaningless and you want that at the TOP of the log.
step "1/3  test suite" \
    pixi run python -m pytest tests/ -q

# 2. Global -- the latitude tilt table is the result. See note 2 above.
step "2/3  F12 four-variable provenance (global)" \
    pixi run python ada/ab_compare.py "$GLOBAL_FILE" \
        --max-timesteps "$GLOBAL_STEPS" \
        --fix four_variable_provenance --out "$OUT"

# 3. North Atlantic -- the domain the replication is actually scored on.
step "3/3  F12 four-variable provenance (NA)" \
    pixi run python ada/ab_compare.py "$NA_FILE" \
        --fix four_variable_provenance --out "$OUT"

echo
echo "############################################################"
echo "### SUMMARY  $(date -Is)"
echo "############################################################"
printf "%-46s %8s %8s\n" "step" "exit" "minutes"
fail=0
for i in "${!NAMES[@]}"; do
    printf "%-46s %8s %8s\n" "${NAMES[$i]}" "${STATUSES[$i]}" "${MINUTES[$i]}"
    [ "${STATUSES[$i]}" -ne 0 ] && fail=1
done
echo
echo "CSV and JSON written to $REPO/$OUT"
echo
echo "READING ORDER: (a) which of the 21 came back bit-identical -- that is the"
echo "list of diagnostics that never touch an archived field; (b) the LATITUDE"
echo "TILT table under step 2, against the tropics-negative prediction on"
echo "record above; (c) step 3's NA flip rates, which are the ones the"
echo "replication is scored on."
echo
echo "NEITHER SIDE IS A BUG FIX. Archived = the IFS prognostic state."
echo "Computed = Prosser's four fields. The output decides which to REPORT,"
echo "not which is right."
exit $fail
