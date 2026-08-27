#!/bin/bash
# ===========================================================================
# 00_smoke_test.sbatch -- prove the environment works on a compute node
# BEFORE launching 504 array tasks.
#
# Submit:
#     cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
#     mkdir -p logs
#     sbatch jobs/00_smoke_test.sbatch
#     # then, when it finishes:
#     cat logs/smoke-*.out
#
# Every check below has already failed at least once somewhere in this
# project's history, which is why each is here rather than assumed.
# ===========================================================================

#SBATCH --job-name=cat-smoke
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/smoke-%j.out
#SBATCH --error=logs/smoke-%j.out

set -uo pipefail

REPO=/scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence

module load 2025
module load pixi
cd "$REPO"   # pixi.toml and .pixi/ live here

echo "=== node $(hostname) at $(date -Is) ==="

echo
echo "--- 1. exactly one task? (the egress test printed everything twice) ---"
echo "SLURM_NTASKS=${SLURM_NTASKS:-unset}  SLURM_PROCID=${SLURM_PROCID:-unset}"
echo "This block must appear ONCE. If it appears twice, --ntasks=1 is not"
echo "taking effect and the download array would double-request every month."

echo
echo "--- 2. python and rojak version ---"
pixi run python -c "
import sys, importlib.metadata as md
print('  python :', sys.version.split()[0])
print('  rojak  :', md.version('rojak-cat'))
assert md.version('rojak-cat').endswith('g25b8685c6'), 'WRONG ROJAK REV -- check pixi.lock'
print('  rev matches the pin in pixi.toml')
"

echo
echo "--- 3. cfgrib / eccodes ---"
pixi run python -m cfgrib selfcheck

echo
echo "--- 4. the 21 diagnostics import and register ---"
pixi run python -c "
import importlib.util, sys
from pathlib import Path
s = importlib.util.spec_from_file_location('diagnostics', Path('2_diagnostics.py'))
m = importlib.util.module_from_spec(s); sys.modules['diagnostics'] = m; s.loader.exec_module(m)
print('  REFERENCE_TABLE:', len(m.REFERENCE_TABLE), 'diagnostics')
assert len(m.REFERENCE_TABLE) == 21
"

echo
echo "--- 5. the verification suite, on the cluster ---"
# If this passes here, the numbers in verification/ are reproducible on ADA
# rather than only in the environment they were first produced in.
pixi run python -m pytest tests -q --no-header -p no:cacheprovider

echo
echo "--- 6. CDS credentials and authentication ---"
if [ ! -f ~/.cdsapirc ]; then
  echo "  FAIL: ~/.cdsapirc missing -- every download task will fail"
else
  echo "  ~/.cdsapirc present, mode $(stat -c '%a' ~/.cdsapirc) (should be 600)"
  pixi run python -c "
import cdsapi
c = cdsapi.Client()
c.check_authentication()
print('  CDS authentication OK from a compute node')
" || echo "  FAIL: CDS auth did not succeed from this node"
fi

echo
echo "--- 7. one real month, end to end (Jan 2016: small, recent, non-leap) ---"
# This is the only step that costs a CDS request. It proves the whole download
# path -- request, atomic .tmp write, integrity check, rename -- on real data.
pixi run python 1_download_hpc.py \
    --trial-day 2016-01-15 \
    --out-dir /scistor/SBE-EDS-ClimateKoopman/yen230/era5_smoke \
  && echo "  download + integrity check PASSED" \
  || echo "  FAIL: download or integrity check did not pass"

echo
echo "=== smoke test finished at $(date -Is) ==="
echo "If every block above passed, submit jobs/01_download.sbatch."
