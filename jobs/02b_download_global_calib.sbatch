#!/bin/bash
# ===========================================================================
# 02b_download_global_calib.sbatch -- year-2000 GLOBAL calibration domain,
# SUB-SAMPLED to 4 days per month.
#
# WHY A SUB-SAMPLE, AND WHY IT IS NOT A COMPROMISE
# ------------------------------------------------
# The calibration thresholds are PERCENTILES of a climatological distribution
# (Prosser 2023 s2, following Williams 2017: the 97th, 99.1st, 99.6th, 99.8th
# and 99.9th percentiles of year 2000, latitudinally weighted). A percentile
# needs coverage of the seasonal and diurnal cycles and enough tail samples --
# it does not need consecutive days.
#
#     days 1, 9, 17, 25 of every month  =  48 days spread evenly across 2000
#     all 8 three-hourly steps kept     =  the full diurnal cycle
#     every gridpoint kept              =  the full global field
#
#     points per diagnostic   3.99e8   (vs 3.04e9 for the whole year)
#     points above p99.6      1.6e6    -- sampling error on the threshold is
#                                         far below the 0.2% gradient noise
#                                         floor of the diagnostics themselves
#
# TWO PROBLEMS THIS ALSO SOLVES
# -----------------------------
# 1. calibration.compute_thresholds ravels the whole calibration field into
#    memory and weighted_percentile then sorts it and calls np.unique. A full
#    global year is 721 x 1440 x 2928 = 3.04e9 points PER DIAGNOSTIC: 24 GB of
#    float64 values, 24 GB of broadcast weights, plus a sort buffer. That does
#    not fit on any partition here. The sub-sample is 4.0e8 points -> ~3.2 GB
#    of values, which does. Without this, stage 2 cannot run at all.
#
# 2. A full global month is 5,208 fields in one CDS request and Copernicus
#    enforces a per-request limit whose current value we have never confirmed
#    (see 02_download_global.sbatch). This request is 672 fields -- comfortably
#    under any plausible limit, so the open question stops being a blocker.
#
# IS THE SUB-SAMPLE BIG ENOUGH? THE DATA ITSELF WILL TELL YOU.
# ------------------------------------------------------------
# Raw point counts overstate the effective sample size, because gridpoints
# within one timestep are spatially correlated. So do not assume -- split the
# 48 days in half (days 1+17 vs days 9+25) and compute thresholds on each half
# separately. If the two halves agree to a few percent, 48 days was enough. If
# they disagree by tens of percent, add days 5, 13, 21, 29 and re-run. That
# check costs nothing extra: it uses only the data this job downloads.
#
# OUTPUT FILENAMES
# ----------------
#     era5_glob_2000-01_d01-09-17-25.grib   ... one per month
#
# Deliberately NOT the era5_glob_YYYY-MM.grib name a full month gets: a partial
# month must never be picked up by anything globbing the complete-month series.
# The integrity check needs no special-casing -- expected_counts() derives the
# timestep count from len(request["day"]), so each file is verified against
# 4 x 8 = 32 timesteps and the global 721 x 1440 grid.
#
# SUBMIT
# ------
#     sbatch --array=0 jobs/02b_download_global_calib.sbatch    # ONE month first
#     # check the log and the file, THEN:
#     sbatch jobs/02b_download_global_calib.sbatch              # the other 11
#
# The downloader is resumable: a month that already has a passing final file is
# skipped, and a month that failed left no final file, so re-submitting the
# whole array retries exactly the months that failed.
#
# SIZING (from the measured 349 MB/day global figure)
#     one month, 4 days    ~1.4 GB
#     all 12 months        ~17 GB      (vs 128 GB for the full year)
# ===========================================================================

#SBATCH --job-name=era5-glob-calib
#SBATCH --partition=defq
#SBATCH --array=0-11%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/global-calib-%A_%a.out
#SBATCH --error=logs/global-calib-%A_%a.out

# --- why these differ from 02_download_global.sbatch -----------------------
# --ntasks=1     MANDATORY, not decorative. srun on this cluster allocated TWO
#                tasks by default and ran the payload twice on one node (see
#                STATUS 10.7). Two processes requesting the same month from CDS
#                would race on the same <name>.tmp, and the atomic-rename guard
#                in _fetch_verified protects against a half-written file, NOT
#                against two concurrent writers.
#
#        %4      Raised from the full-month script's %2. Each request here is
#                ~1.4 GB rather than ~10.6 GB, so four at a time is less load
#                on CDS and on this end than two full months were.
#
# --mem=8G       Down from 16G. The integrity check opens the finished GRIB and
#                counts dimensions; these files are 1/8 the size of a full
#                global month. Memory is scheduled explicitly on ADA, so asking
#                for less is what lets four array tasks land concurrently.
#
# --time=24:00:00  Generous on purpose. The transfer is short but the CDS queue
#                is not predictable, and cdsapi is now bounded to 8 attempts
#                (STATUS 10.10) so a genuinely unreachable CDS fails loudly
#                rather than holding the node to walltime.
# ---------------------------------------------------------------------------

set -euo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
OUT=$BASE/raw/global

DOMAIN=global
START=2000-01
END=2000-12
DAYS=1,9,17,25

module load 2025
module load pixi          # lowercase -- 'module load Pixi' fails (STATUS 10.2)

mkdir -p "$OUT" "$REPO/logs"
cd "$REPO"                # pixi.toml and .pixi/ live here

echo "=== task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -Is) ==="
echo "SLURM_NTASKS=${SLURM_NTASKS:-unset}   (must be 1)"
if [ "${SLURM_NTASKS:-1}" != "1" ]; then
    echo "FATAL: SLURM_NTASKS is not 1 -- two processes would race on the same" >&2
    echo "       .tmp file for this month. Refusing to run." >&2
    exit 1
fi
df -h "$BASE" | tail -1

pixi run python 1_download_hpc.py \
    --start  "$START" \
    --end    "$END" \
    --domain "$DOMAIN" \
    --days   "$DAYS" \
    --index  "$SLURM_ARRAY_TASK_ID" \
    --out-dir "$OUT"

echo "=== done at $(date -Is) ==="
ls -lh "$OUT" | tail -5
