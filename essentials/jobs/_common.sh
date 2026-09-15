# Sourced by every job script. One place for the ADA set-up lessons (STATUS 10-11, 15.7).
BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence            # the pixi environment lives here
module load 2025
module load pixi                      # lowercase; 'Pixi' does not exist
export PYTHONUNBUFFERED=1             # otherwise a running job looks hung in the log
cd "$REPO" || { echo "cannot cd to $REPO"; exit 1; }
if [ "${SLURM_NTASKS:-1}" != "1" ]; then
    echo "FATAL: SLURM_NTASKS=${SLURM_NTASKS}; two copies would race on the same output" >&2
    exit 1
fi
echo "=== ${SLURM_JOB_NAME:-job} ${SLURM_ARRAY_TASK_ID:-} on $(hostname) at $(date -Is) ==="
run() { pixi run python essentials/main.py "$@"; }
