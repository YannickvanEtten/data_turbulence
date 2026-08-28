#!/bin/bash
# ===========================================================================
# 06_calibration_check.sbatch -- the go/no-go gate for the 504-month run.
#
#     sbatch jobs/06_calibration_check.sbatch
#     sbatch jobs/06_calibration_check.sbatch --skip-split-half   # ~half the time
#
# Calibrates severity thresholds on the global year 2000, applies them to
# North Atlantic DJF 1979 inside Prosser's own box (36-60N / 55-10W), and
# compares the exceedance frequencies against Prosser (2023) Table 1.
#
# WHY IT IS A JOB AND NOT AN INTERACTIVE RUN
# ------------------------------------------
# ADA's login nodes auto-log-off after ~5 minutes of anything heavy
# (STATUS 10.1). This sorts 21 arrays of 4e8 float64 values.
#
# --mem=32G
#     calibration.compute_thresholds loops over diagnostics and the zarr
#     stores are opened lazily, so the peak is set by ONE diagnostic, not by
#     the 33 GB of calibration data. One diagnostic is 4.0e8 float32 (1.6 GB),
#     which weighted_percentile casts to float64 and sorts: values 3.2 GB,
#     broadcast weights 3.2 GB, plus numpy's sort and np.unique buffers.
#     ~20 GB with headroom to 32.
#
#     32G still fits EVERY defq node including the 59 GB ones, which is what
#     keeps it schedulable. A 120G request on this cluster queued 22 hours.
#
# --time=04:00:00
#     Estimated ~50 min for the full-year percentiles, ~50 min for the two
#     split halves, ~20 min for the identity check, minutes for the rest.
#     4 h is roughly 2x that. Pass --skip-split-half to halve it.
#
# OUTPUT
#     calibration/thresholds_<date>.json   (via calibration.save_thresholds,
#                                           schema-versioned, with provenance)
#     everything else to the log
# ===========================================================================

#SBATCH --job-name=cat-calib
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/calib-%j.out
#SBATCH --error=logs/calib-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1        # see 04_diagnostics_global.sbatch

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== calibration check on $(hostname) at $(date -Is) ==="

start=$(date +%s)
pixi run python ada/calibration_check.py --base "$BASE" "$@"
status=$?
end=$(date +%s)

echo
echo "=== exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
if [ $status -eq 0 ]; then
  echo "Threshold file:"
  ls -lh "$BASE/calibration/" 2>/dev/null | sed 's/^/  /'
fi
exit $status
