#!/bin/bash
# ===========================================================================
# 01_download.sbatch -- download ERA5 months from CDS, one array task per month.
#
# Submit:
#     cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
#     mkdir -p logs
#     sbatch jobs/01_download.sbatch
#
# Resume after any failure or cancellation: submit the SAME command again.
# 1_download_hpc.py skips a month only if its final file EXISTS AND re-passes
# the integrity check, so completed months cost seconds and failed ones are
# retried. Nothing needs editing between attempts.
#
# Check progress:
#     squeue --me
#     ls $BASE/raw/north_atlantic/era5_na_*.grib | wc -l          # of 504
#     grep -l Error logs/download-*.out                                # failures
# ===========================================================================

#SBATCH --job-name=era5-dl
#SBATCH --partition=defq
#SBATCH --array=0-503%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/download-%A_%a.out
#SBATCH --error=logs/download-%A_%a.out

# --- why these numbers -----------------------------------------------------
# --array=0-503   504 months, 1979-01..2020-12. --index maps straight onto
#                 SLURM_ARRAY_TASK_ID; download_plan.month_span defines the
#                 order, so the mapping is deterministic and re-runnable.
#
#         %4      AT MOST 4 TASKS AT ONCE. This is a CDS courtesy limit, not a
#                 cluster one -- Copernicus throttles per user, and firing 504
#                 concurrent requests gets you queued or rate-limited, not
#                 served faster. Raise cautiously after watching the first
#                 batch; the wall-clock is dominated by CDS queueing either way.
#
# --ntasks=1      NOT OPTIONAL. Without it srun/sbatch can allocate more than
#                 one task and run this script twice on the same node -- which
#                 was visible in the egress test, where every line printed
#                 twice. Two tasks downloading the same month would issue two
#                 CDS requests and race on the same <name>.tmp file, and the
#                 atomic-rename guard in _fetch_verified protects against a
#                 half-written file, NOT against two writers.
#
# --time=12:00:00 One month is ~373 MB and downloads in minutes, but a CDS
#                 request can sit in Copernicus's queue for hours. The job is
#                 idle-waiting on a remote queue, not computing. Well under the
#                 7-day partition cap.
#
# --mem=4G        Only the integrity check needs memory: it opens the finished
#                 GRIB with cfgrib to count variables, levels and timesteps.
# ---------------------------------------------------------------------------

set -euo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
OUT=$BASE/raw/north_atlantic

DOMAIN=north_atlantic
START=1979-01
END=2020-12

module load 2025
module load pixi          # lowercase -- 'module load Pixi' does not exist

mkdir -p "$OUT"
cd "$REPO"   # pixi.toml and .pixi/ live here

echo "=== task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -Is) ==="

pixi run python 1_download_hpc.py \
    --start "$START" \
    --end   "$END" \
    --domain "$DOMAIN" \
    --index "$SLURM_ARRAY_TASK_ID" \
    --out-dir "$OUT"

echo "=== done at $(date -Is) ==="
