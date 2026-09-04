#!/bin/bash
# ---------------------------------------------------------------------------
# 18_merge_global_months.sbatch
#
# Concatenate the three ~10-day blocks jobs/17 downloaded into one contiguous
# GRIB per month, so ada/diagnostics_global.py sees a real month.
#
# WHY `cat` IS THE RIGHT TOOL AND NOT A HACK
#     GRIB is a sequence of self-describing messages. Byte-concatenating two
#     valid GRIB files produces a valid GRIB file, and cfgrib builds its
#     hypercube from the message keys, not from file order. No decode, no
#     re-encode, no precision loss, no new code path to verify.
#
# WHY CONTIGUITY MATTERS AT ALL -- IT IS ONE DIAGNOSTIC
#     frontogenesis_isentropic (f2d) is the only one of the 21 taking a time
#     derivative, via .differentiate("time", datetime_unit="s"). On the
#     non-contiguous 48-day sub-sample its centred stencil straddles 8-day
#     gaps at 6 of every 32 timesteps. np.gradient's non-uniform formula
#     degrades gracefully there (it becomes a one-sided 3-hourly difference
#     plus a sub-percent contamination term), so this was never the
#     order-of-magnitude corruption it looks like -- but a contiguous year
#     removes it for free, and f2d is the diagnostic with the least margin.
#
# ⚠ RUN THIS AS A JOB, NOT ON A LOGIN NODE. Login nodes auto-log-off after
#   ~5 minutes of anything heavy (STATUS 10.1), and this moves ~128 GB.
#
# ⚠ THE BLOCKS ARE NOT DELETED HERE, DELIBERATELY. Keeping them costs a
#   transient ~128 GB and means a bad merge is recoverable without a second
#   day of CDS queue. Delete them by hand only after jobs/19 has produced
#   derived output from the merged files:
#       rm $BASE/raw/global/era5_glob_2000-??_d*.grib*
#
# ⚠ MOVE THE OLD SUB-SAMPLE OUT FIRST, or its filename matches the d01-* glob:
#       mkdir -p $BASE/raw/global_sub48
#       mv $BASE/raw/global/era5_glob_2000-??_d01-09-17-25.grib* $BASE/raw/global_sub48/
#   The script refuses to run if it finds it still in place.
#
#   sbatch jobs/18_merge_global_months.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=cat-merge-glob
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/merge-glob-%j.out
#SBATCH --error=logs/merge-glob-%j.out

set -euo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
RAW=$BASE/raw/global
YEAR=2000

echo "=== merge global ${YEAR} blocks on $(hostname) at $(date -Is) ==="

shopt -s nullglob
stale=( "$RAW"/era5_glob_${YEAR}-??_d01-09-17-25.grib )
if [ ${#stale[@]} -gt 0 ]; then
    echo "FATAL: the 48-day sub-sample is still in $RAW:" >&2
    printf '   %s\n' "${stale[@]}" >&2
    echo "Move it to $BASE/raw/global_sub48/ first -- otherwise it matches the" >&2
    echo "d01-* glob below and would be concatenated into January." >&2
    exit 1
fi

fail=0
for M in $(seq -w 1 12); do
    OUT="$RAW/era5_glob_${YEAR}-${M}.grib"
    b0=( "$RAW"/era5_glob_${YEAR}-${M}_d01-*.grib )
    b1=( "$RAW"/era5_glob_${YEAR}-${M}_d11-*.grib )
    b2=( "$RAW"/era5_glob_${YEAR}-${M}_d21-*.grib )

    if [ ${#b0[@]} -ne 1 ] || [ ${#b1[@]} -ne 1 ] || [ ${#b2[@]} -ne 1 ]; then
        echo "  ${YEAR}-${M}  MISSING BLOCKS (${#b0[@]}/${#b1[@]}/${#b2[@]}) -- skipped"
        fail=1
        continue
    fi
    if [ -f "$OUT" ]; then
        echo "  ${YEAR}-${M}  already merged, skipping ($(du -h "$OUT" | cut -f1))"
        continue
    fi

    want=0
    for f in "${b0[0]}" "${b1[0]}" "${b2[0]}"; do
        want=$(( want + $(stat -c%s "$f") ))
    done

    # Write to .tmp and rename, so an interrupted merge can never leave a
    # truncated file sitting under the name the diagnostics job will trust.
    cat "${b0[0]}" "${b1[0]}" "${b2[0]}" > "$OUT.tmp"
    got=$(stat -c%s "$OUT.tmp")
    if [ "$got" -ne "$want" ]; then
        echo "  ${YEAR}-${M}  SIZE MISMATCH: got $got want $want -- not renaming"
        rm -f "$OUT.tmp"
        fail=1
        continue
    fi
    mv "$OUT.tmp" "$OUT"
    echo "  ${YEAR}-${M}  merged -> $(du -h "$OUT" | cut -f1)"
done

echo
echo "merged months present: $(ls -1 "$RAW"/era5_glob_${YEAR}-??.grib 2>/dev/null | wc -l) of 12"
df -h "$BASE" | tail -1
echo "=== done at $(date -Is), fail=$fail ==="

# A byte count proves the concatenation, not that the GRIB decodes. Before
# jobs/19, confirm the timestep and grid counts with the project's own checker:
#     pixi run python ada/check_calib_downloads.py --help
exit $fail
