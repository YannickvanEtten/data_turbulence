#!/bin/bash
# ===========================================================================
# ada/01_setup.sh -- build the project on ADA. IDEMPOTENT: safe to re-run.
#
# WHERE TO RUN THIS: an INTERACTIVE node, never the login node.
# `pixi install` takes several minutes and login nodes auto-log-off sessions
# doing heavy work every 5 minutes. Get an interactive shell first:
#
#     srun --partition=defq --ntasks=1 --cpus-per-task=4 --mem=8G \
#          --time=01:00:00 --pty bash
#
#     # (or, if your site allows it:  ssh inter01)
#
# then:
#     cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
#     bash ada/01_setup.sh
#
# WHAT IT DOES
#   1. creates the directory tree (mkdir -p -- never touches existing data)
#   2. checks the code is present
#   3. builds the pixi environment INSIDE the repo, from pixi.toml
#   4. verifies: cfgrib/eccodes, rojak rev, all 21 diagnostics, 25 tests
#
# WHAT IT NEVER DOES
#   deletes anything, overwrites data, or downloads from CDS.
# ===========================================================================

set -uo pipefail

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
REPO=$BASE/data_turbulence
EXPECTED_ROJAK_SUFFIX="g25b8685c6"

fail() { printf '\n!! %s\n' "$1"; exit 1; }
step() { printf '\n=== %s ===\n' "$1"; }

echo "ADA setup -- $(date -Is)  on $(hostname)"
if [ -n "${SLURM_JOB_ID:-}" ]; then
  echo "running as SLURM job ${SLURM_JOB_ID} -- good, this is a compute node"
fi
case "$(hostname)" in
  login*) echo
          echo "!! You are on a LOGIN node. pixi install will likely be killed."
          echo "!! Get an interactive shell first (see the header of this file),"
          echo "!! then re-run. Continuing anyway in 10s -- Ctrl-C to stop."
          sleep 10 ;;
esac

# ---------------------------------------------------------------------------
step "1. directory tree"
# mkdir -p is a no-op on directories that already exist, so re-running this
# script can never disturb data that is already downloaded.
for d in raw/north_atlantic raw/global derived/north_atlantic derived/global \
         calibration results logs; do
  mkdir -p "$BASE/$d" && echo "  ok  $BASE/$d"
done
echo
df -h "$BASE" | sed 's/^/  /'

# ---------------------------------------------------------------------------
step "2. code"
[ -d "$REPO" ] || fail "no code at $REPO.
   Put it there first, either:
     git clone https://github.com/YannickvanEtten/data_turbulence.git $REPO
   (after committing and pushing your local work -- the GitHub copy is stale)
   or by copying the folder across with MobaXterm's SFTP pane."
[ -f "$REPO/pixi.toml" ] || fail "$REPO/pixi.toml missing -- the environment
   definition is what this script builds from. Make sure the copy on ADA is
   current, not the May version from GitHub."
cd "$REPO" || fail "cannot cd to $REPO"
echo "  repo: $REPO"
if [ -d .git ]; then
  echo "  commit: $(git log -1 --format='%h %ad %s' --date=short 2>/dev/null)"
  n=$(git status --porcelain 2>/dev/null | wc -l)
  [ "$n" -gt 0 ] && echo "  note: $n uncommitted change(s) here"
fi
for f in 1_download_hpc.py download_plan.py 2_diagnostics.py 3_pipeline.py calibration.py; do
  [ -f "$f" ] || fail "missing $f -- the copy on ADA is incomplete"
done
echo "  all key files present"

# ---------------------------------------------------------------------------
step "3. modules"
module load 2025 || fail "module load 2025 failed"
module load pixi || fail "module load pixi failed (lowercase 'pixi', not 'Pixi')"
echo "  $(pixi --version)"

# ---------------------------------------------------------------------------
step "4. environment"
# pixi install is incremental: if .pixi/ already matches pixi.lock this returns
# in seconds, so re-running costs nothing.
if [ -d .pixi ]; then
  echo "  .pixi/ exists ($(du -sh .pixi 2>/dev/null | cut -f1)) -- refreshing"
else
  echo "  building from scratch. THIS TAKES SEVERAL MINUTES. Do not interrupt."
fi
pixi install || fail "pixi install failed.
   Most common cause: no outbound internet on this node, or a conda-forge
   hiccup. Retry once; if it fails again, check
     curl -sSI https://conda.anaconda.org/conda-forge/noarch/repodata.json"
echo "  environment ready ($(du -sh .pixi 2>/dev/null | cut -f1))"
[ -f pixi.lock ] && echo "  pixi.lock present -- COMMIT IT, it is what pins the environment"

# ---------------------------------------------------------------------------
step "5. verification"
# NOTE: every check below runs from a FILE, never `pixi run python -c "..."`.
# pixi parses the -c string with its own task-shell parser before python ever
# sees it, and that parser rejects ordinary English words that happen to be
# shell reserved words. The first real run on ADA died on the word "in" inside
# an assertion message, with the environment perfectly healthy. Files have no
# quoting layer to get wrong.

echo
echo "-- 5a. python, rojak rev, GRIB stack, diagnostics, CDS --"
pixi run python ada/verify_env.py || fail "environment checks failed (see above)"

echo
echo "-- 5b. cfgrib selfcheck --"
pixi run python -m cfgrib selfcheck || fail "cfgrib selfcheck failed"

echo
echo "-- 5c. the verification suite, on ADA --"
# Passing here means the numbers under verification/ are reproducible in the
# environment that will actually produce the dataset, not only in the one
# where they were first measured.
pixi run python -m pytest tests -q --no-header -p no:cacheprovider \
  || fail "the test suite did not pass on ADA -- do not proceed to downloads"

# ---------------------------------------------------------------------------
step "SETUP COMPLETE"
cat <<'NEXT'
  Next, in order:

    1. Accept the ERA5 pressure-levels licence (one-time, on the website):
         https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels
       -> "Terms of use" tab -> Accept. The first download fails without it.

    2. Smoke test, including one real trial-day download:
         sbatch jobs/00_smoke_test.sbatch
         cat logs/smoke-*.out

    3. The pilot month:
         sbatch --array=0 jobs/01_download.sbatch
         # index 0 = era5_na_1979-01.grib

  Do NOT launch the full array until 2 and 3 have both come back clean.
NEXT
