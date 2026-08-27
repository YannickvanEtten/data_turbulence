#!/bin/bash
# ===========================================================================
# 00_setup_env.sbatch -- build and verify the environment as a BATCH job.
#
# Use this instead of holding an interactive session when defq has no idle
# nodes (check with `sinfo -s`: the NODES column is Allocated/Idle/Other/Total,
# so "16/0/1/17" means nothing is free and an interactive srun would make you
# wait at the terminal for one).
#
# A batch job queues exactly the same, but it runs whenever a slot opens and
# you can close the laptop meanwhile.
#
# Submit:
#     cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
#     mkdir -p logs
#     sbatch jobs/00_setup_env.sbatch
#     squeue --me                      # PD = pending, R = running
#     cat logs/setup-*.out             # when it's gone from the queue
#
# Safe to re-run: everything it does is idempotent.
# ===========================================================================

#SBATCH --job-name=cat-setup
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/setup-%j.out
#SBATCH --error=logs/setup-%j.out

# --cpus-per-task=4  pixi resolves and links packages in parallel; 4 is a
#                    reasonable ask for a one-off build and does not make the
#                    job noticeably harder to schedule than 1.
# --time=01:00:00    generous. A cold build is minutes; the ceiling is only
#                    there so a hung download cannot occupy a node all week.

set -uo pipefail

REPO=/scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }

echo "=== environment build on $(hostname) at $(date -Is) ==="
echo

bash ada/01_setup.sh
status=$?

echo
echo "=== 01_setup.sh exited with status $status at $(date -Is) ==="
exit $status
