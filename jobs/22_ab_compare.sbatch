#!/bin/bash
# ---------------------------------------------------------------------------
# 22_ab_compare.sbatch
#
# Stage 1 of the post-audit plan: measure how much each proposed fix actually
# moves the exceedance field, BEFORE committing anything to a 7-hour re-derive.
# Seven checks in one overnight job, one log to read in the morning.
#
# WHAT IT ANSWERS, AND WHY IT IS ONE JOB RATHER THAN SIX
#     Every check is a same-input A/B: compute the 21 twice, once with the fix
#     and once without, and report the SYMMETRIC DIFFERENCE OF THE EXCEEDANCE
#     SET. Thresholds are percentiles of our own data, so a constant scale
#     error cancels exactly (AUDIT §5.2) -- rho and relative differences cannot
#     tell you whether a fix matters and the flip rate can. The runs are
#     independent and could go as seven array tasks, but sequential in one job
#     costs ~3 h of an idle night and produces a single ordered log.
#
#     PRE-REGISTERED DECISION RULE, in ada/ab_compare.py and fixed before the
#     numbers: flip < 0.5 % at every severity = inert, merge on the citation
#     alone; flip >= 0.5 % anywhere = material, goes into ONE batched re-run.
#
# NO NEW DATA. Everything here reads GRIB already on the share. No CDS.
#
# WHY THE GLOBAL FILE FOR CHECK 2 AND THE NA FILE FOR THE REST
#     F1 (the meridional metric) is a LATITUDE-STRUCTURED error: +0.42 % at
#     30N falling to -0.27 % at 75N. Its signature is a monotone ramp across
#     the latitude-band table, which a 30-60N box cannot show. It needs
#     latitudes, so it gets a global month -- CAPPED with --max-timesteps.
#
#     THE CAP IS NOT OPTIONAL. raw/global/ now holds FULL months (the 4-day
#     sub-sampled files from STATUS 13.1 were superseded by the full-year
#     download of 15). era5_glob_2000-01.grib is 31 days x 8 = 248 timesteps at
#     721x1440, and this tool holds TWO complete result sets: 21 diagnostics x
#     248 x 721 x 1440 x 4 bytes x 2 = ~43 GB of diagnostic arrays before a
#     single intermediate. 32 timesteps is 33 million cells per diagnostic --
#     the p99.9 set alone is ~33,000 cells, an order of magnitude above
#     MIN_EXCEED_CELLS -- and a flip rate needs cells, not a long record.
#     F5 needs the global file for a different reason: NCSU1 is calibrated
#     against GLOBAL percentiles, so whether its Ri floor reaches the tail is a
#     global question. The rest are single-diagnostic substitutions whose flip
#     rate is domain-independent, so the cheaper NA month does.
#
# SIZING -- 48G, AND THE REASON IS SCHEDULING, NOT ARITHMETIC
#     ab_compare holds BOTH result sets at once, so peak is ~2x a normal pass.
#     The global calibration month is the binding case: one pass measured 24G
#     (job 1092495), the second result set adds ~6 GB, and the float64
#     flattening inside compare_one adds ~1 GB transient. Call it ~32 GB.
#
#     48G rather than 32G for headroom, and NOT MORE, because STATUS 11.9 is
#     unambiguous: memory is the biggest lever on time-to-start, a 120 GB
#     request was scheduled 22 HOURS out where 24 GB started in 13 seconds.
#     defq is heterogeneous and five of its seventeen nodes hold only 59 GB --
#     48G fits every node in the partition, 96G fits thirteen of them. Same
#     request as jobs/16, which ran 4 h 06 m without trouble.
#
#     Walltime 4 h for the same reason: "short, honest walltimes get
#     backfilled; 6 h waits, 1 h runs". The battery is ~3 h.
#
#   mkdir -p logs                      # SLURM opens the log BEFORE the script runs
#   sbatch jobs/22_ab_compare.sbatch
#   sbatch --export=ALL,NA_FILE=/path/era5_na_1995-01.grib jobs/22_ab_compare.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-ab
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/ab-%j.out
#SBATCH --error=logs/ab-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
OUT=data_prosser/ab

module load 2025
module load pixi

export PYTHONUNBUFFERED=1
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }
# logs/ must already exist -- SLURM opens the output file before this script
# runs, so creating it here is too late (STATUS 11.9). Kept anyway for the
# output directory, which is only touched by the steps below.
mkdir -p logs "$OUT"

GLOBAL_STEPS=${GLOBAL_STEPS:-32}   # see the cap note above; >= 3 required by f2d

# --- locate the inputs, and say which were picked --------------------------
# Overridable with --export=ALL,GLOBAL_FILE=...,NA_FILE=...
if [ -z "${GLOBAL_FILE:-}" ]; then
    # Name it exactly. A glob here would match the full month and the three
    # 10-day chunks alike, and `head -1` would silently pick whichever sorts
    # first -- a 248-timestep file where a capped one was intended.
    GLOBAL_FILE=$BASE/raw/global/era5_glob_2000-01.grib
fi
if [ -z "${NA_FILE:-}" ]; then
    NA_FILE=$BASE/raw/north_atlantic/era5_na_2000-01.grib
fi

for v in GLOBAL_FILE NA_FILE; do
    f=${!v}
    if [ -z "$f" ] || [ ! -f "$f" ]; then
        echo "!! $v not found (got '${f:-<empty>}')."
        echo "   Set it explicitly, e.g.:"
        echo "   sbatch --export=ALL,GLOBAL_FILE=\$BASE/raw/global/<file>.grib \\"
        echo "                   ,NA_FILE=\$BASE/raw/north_atlantic/<file>.grib \\"
        echo "          jobs/22_ab_compare.sbatch"
        exit 2
    fi
done

echo "=== Stage-1 A/B battery on $(hostname) at $(date -Is) ==="
echo "  repo        $REPO"
echo "  global file $GLOBAL_FILE  (first $GLOBAL_STEPS timesteps)"
echo "  NA file     $NA_FILE"
echo "  output      $REPO/$OUT"
echo

declare -a NAMES STATUSES MINUTES

step () {                       # step <label> <command...>
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
    return 0                    # never abort the battery on one failure
}

# 1. The suite must pass first. F4 changed compute_rojak_diagnostics' signature
#    and F8 added two guards in prepare_for_rojak; if either broke a pinned
#    test, every number below is suspect and you want to know that at the top
#    of the log rather than at the bottom.
step "1/7  test suite" \
    pixi run python -m pytest tests/ -q

# 2. F1, on the global month -- the latitude-tilt table is the actual test of
#    AUDIT §3.1's derivation. Expect ~13 of 21 to move and 8 to come back
#    bit-identical (the 8 with no horizontal d/dy). A monotone ramp across the
#    bands confirms the derivation; a flat row falsifies it.
step "2/7  F1 meridional metric (global)" \
    pixi run python ada/ab_compare.py "$GLOBAL_FILE" \
        --max-timesteps "$GLOBAL_STEPS" --fix meridional_metric --out "$OUT"

# 3. F2 -- UBF's vorticity provenance. Prediction on record (AUDIT §6 F2): if
#    the mixed-operator diagnosis is right this is MATERIAL and large; if it
#    comes back inert, the cause is the spherical Jacobian convention instead.
step "3/7  F2 UBF computed vorticity (NA)" \
    pixi run python ada/ab_compare.py "$NA_FILE" \
        --fix ubf_computed_vorticity --out "$OUT"

# 4. F3 -- |PV| from Sharman A18 instead of ERA5's archived full Ertel PV.
#    Note Storer et al. (2017) dropped |PV| from this basket entirely, "because
#    it was found to give unrealistic results", so treat a large flip here as
#    informative about the diagnostic rather than about the pipeline.
step "4/7  F3 magnitude_pv from A18 (NA)" \
    pixi run python ada/ab_compare.py "$NA_FILE" \
        --fix pv_sharman_a18 --out "$OUT"

# 5. F10 -- Endlich's directional shear: rojak's angle-secant against the exact
#    closed form |u dv/dz - v du/dz|/|v|. Measured 67 % flip on the validation
#    subset, so this is expected to be MATERIAL; the question the real month
#    answers is whether it moves the level ratio off 0.51.
step "5/7  F10 Endlich closed-form shear (NA)" \
    pixi run python ada/ab_compare.py "$NA_FILE" \
        --fix endlich_component_shear --out "$OUT"

# 6. F6 -- the frontogenesis reading. C is what the archive holds; A is what
#    Williams & Storer (2022) Eq. (3) p. 1427 actually writes. Both sides are
#    computed here, so this measures the cost of the default change rather than
#    assuming it. rho ~ 0.12 on the validation subset: these are barely the
#    same ranking.
step "6/7  F6 frontogenesis C vs A (NA)" \
    pixi run python ada/ab_compare.py "$NA_FILE" \
        --f2d-variant-baseline C --f2d-variant-variant A --out "$OUT"

# 7. F5 -- NCSU1's Richardson floor. Not an A/B: a measurement of whether
#    MAX(Ri,1e-5) reaches the p97 tail. Global, because NCSU1 is calibrated on
#    global percentiles. The script prints its own verdict.
step "7/7  F5 NCSU1 Ri-floor probe (global)" \
    pixi run python ada/ncsu1_floor_probe.py "$GLOBAL_FILE" \
        --max-timesteps "$GLOBAL_STEPS" --out "$OUT"

# --- one-screen summary ----------------------------------------------------
echo
echo "############################################################"
echo "### SUMMARY  $(date -Is)"
echo "############################################################"
printf "%-42s %8s %8s\n" "step" "exit" "minutes"
fail=0
for i in "${!NAMES[@]}"; do
    printf "%-42s %8s %8s\n" "${NAMES[$i]}" "${STATUSES[$i]}" "${MINUTES[$i]}"
    [ "${STATUSES[$i]}" -ne 0 ] && fail=1
done
echo
echo "CSV and JSON written to $REPO/$OUT"
echo "Read the 'MATERIAL / inert / untouched' line of each A/B block, then the"
echo "LATITUDE TILT table under step 2, then step 7's DECISIVE line."
echo "Anything MATERIAL goes into ONE batched re-derive; anything inert merges"
echo "on its citation without re-deriving 42 years."
exit $fail
