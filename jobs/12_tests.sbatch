#!/bin/bash
# ===========================================================================
# 12_tests.sbatch -- run the verification suite on a compute node.
#
#     mkdir -p logs && sbatch jobs/12_tests.sbatch
#     sbatch jobs/12_tests.sbatch -k TestF2dVariants      # pytest args pass through
#
# WHY A JOB SCRIPT AND NOT srun OR THE LOGIN NODE
# -----------------------------------------------
# Two ADA facts, both already in STATUS.md, both easy to forget:
#
#   §10.2  Modules do not persist. A bare `srun ... pixi run ...` inherits the
#          login shell's environment but not the module system's, and dies with
#          "execve(): pixi: No such file or directory" on the compute node.
#          Every job script re-loads its own modules; this one is no exception.
#
#   §10.1  Login nodes auto-log-off after ~5 minutes of anything heavy. The
#          suite is only ~30 s, but it imports rojak, numba and metpy, and
#          there is no reason to gamble a session on it.
#
# This is also the cheapest possible place for the two 2026-08-29 code changes
# (FORMULA_AUDIT.md §4 and §5) to fail. Run it BEFORE queuing 09/10/11 --
# if the deformation or f2d edits are wrong, nothing downstream is worth
# looking at.
#
# --mem=8G     The manufactured atmosphere in tests/synthetic.py is small; the
#              memory goes on importing rojak's dependency stack, not on data.
#              STATUS.md §11.9: small requests start in seconds, large ones
#              wait hours. Do not round this up.
# --time=00:20 The suite runs in well under a minute. Short honest walltimes
#              get backfilled.
# ===========================================================================

#SBATCH --job-name=cat-tests
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/tests-%j.out
#SBATCH --error=logs/tests-%j.out

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence

module load 2025
module load pixi
export PYTHONUNBUFFERED=1

cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== test suite on $(hostname) at $(date -Is) ==="
echo "    repo $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
echo

start=$(date +%s)
pixi run python -m pytest tests/ -q "$@"
status=$?
end=$(date +%s)

echo
echo "=== exit ${status}, wall $(( (end-start)/60 )) min $(( (end-start)%60 )) s ==="
if [ $status -ne 0 ]; then
    echo
    echo "FAILURES ABOVE. Do not submit jobs/09-11 until this is green:"
    echo "  TestDeformationUnsquared  -> the 2_diagnostics/3_pipeline pairing"
    echo "                               (FORMULA_AUDIT.md 5)"
    echo "  TestF2dVariants           -> the A9 variant flag (FORMULA_AUDIT.md 4)"
    echo "  anything else             -> the change broke something it should"
    echo "                               not have touched"
fi
exit $status
