#!/bin/bash
# ada_recon.sh -- read-only reconnaissance of the ADA setup.
#
# Purpose: the ADA configuration for this project is not recorded anywhere in
# the local repo or project folder. This script recovers it from the cluster
# itself and writes ONE file to bring back: ~/ada_recon.txt
#
# It is READ-ONLY. It creates nothing except that output file, changes no
# settings, submits no real work. The single srun is a 30-second test job whose
# only purpose is to find out whether compute nodes can reach the internet --
# which decides whether downloading can happen inside jobs or must be done on
# the interactive nodes.
#
# Usage:
#     ssh ada-login
#     bash ada_recon.sh          # ~1-2 minutes, most of it the test job queueing
#     # then copy ~/ada_recon.txt back
#
# If the srun step hangs because the queue is busy, Ctrl-C it -- everything
# before it has already been written.

OUT=~/ada_recon.txt
SCISTOR=/scistor/SBE-EDS-ClimateKoopman
exec > >(tee "$OUT") 2>&1

section() { printf '\n\n========== %s ==========\n' "$1"; }

echo "ADA recon -- $(date -Iseconds)"
echo "user=$(whoami)  host=$(hostname)"

section "1. PRIOR JOB HISTORY (recovers earlier testing even if scripts are gone)"
# SLURM's accounting DB remembers submissions long after the scripts vanish.
sacct -u "$USER" -S 2025-01-01 \
      --format=JobID%14,JobName%24,Partition%12,State%14,Elapsed%12,ReqMem%10,Submit%20 \
  | head -80 || echo "sacct unavailable"

section "2. EXISTING CONFIG AND SCRIPTS IN \$HOME"
ls -la ~ 2>/dev/null | head -40
echo "--- any job scripts or shell scripts anywhere in \$HOME ---"
find ~ -maxdepth 4 \( -name '*.sbatch' -o -name '*.slurm' -o -name '*.sh' \) \
     -not -path '*/envs/*' -not -path '*/.pixi/*' -not -path '*/miniconda*' \
     2>/dev/null | head -40
echo "--- ssh config (host aliases) ---"
sed -n '1,60p' ~/.ssh/config 2>/dev/null || echo "no ~/.ssh/config"
echo "--- CDS credentials present? (key itself NOT printed) ---"
if [ -f ~/.cdsapirc ]; then
  echo "~/.cdsapirc EXISTS; url line:"; grep -i '^url' ~/.cdsapirc
else
  echo "~/.cdsapirc NOT FOUND -- the downloader will fail until this exists"
fi
echo "--- module/conda lines in shell startup files ---"
grep -nE 'module |conda |pixi|PATH=' ~/.bashrc ~/.bash_profile ~/.profile 2>/dev/null | head -30

section "3. PARTITIONS, LIMITS, HARDWARE"
/ada-software/ada-info.sh 2>/dev/null | head -60 || echo "ada-info.sh not found"
echo "--- sinfo: partition / timelimit / cpus / memory / nodes ---"
sinfo -o "%20P %12l %6c %10m %6D %20f" 2>/dev/null | head -30
echo "--- association limits for this user (max jobs, QoS) ---"
sacctmgr -n show assoc user="$USER" \
         format=Account%24,Partition%16,QOS%20,MaxJobs,MaxSubmit,GrpTRES%30 2>/dev/null | head -20

section "4. STORAGE AND QUOTA"
echo "--- \$HOME ---"; echo "HOME=$HOME"; df -h "$HOME" 2>/dev/null
echo "--- project folder ---"
ls -ld "$SCISTOR" 2>/dev/null || echo "cannot stat $SCISTOR"
ls -la "$SCISTOR" 2>/dev/null | head -25
echo "--- my subfolder ---"
ls -la "$SCISTOR/$USER" 2>/dev/null | head -40 || echo "no $SCISTOR/$USER"
echo "--- free space where the 42-year run would land ---"
df -h "$SCISTOR" 2>/dev/null
echo "--- quota (whichever tool exists) ---"
quota -s 2>/dev/null || echo "(no quota cmd)"
lfs quota -h "$SCISTOR" 2>/dev/null || true
echo "--- current usage of the project folder (may take a moment) ---"
timeout 90 du -sh "$SCISTOR" 2>/dev/null || echo "(du timed out or denied -- fine)"

section "5. SOFTWARE STACK"
module --version 2>&1 | head -3
module load 2025 2>&1 | head -5
echo "--- Python modules available ---"
module -t avail Python 2>&1 | grep -i '^python/' | head -20
echo "--- other relevant modules ---"
for m in Miniconda3 Pixi ecCodes eccodes netCDF cdo NCO git; do
  printf '%-14s: ' "$m"; module -t avail "$m" 2>&1 | grep -iv '^/' | head -3 | tr '\n' ' '; echo
done
echo "--- python currently on PATH ---"
which python python3 2>/dev/null; python3 --version 2>&1

section "6. INTERNET ACCESS -- the question that shapes the whole pipeline"
echo "### from THIS node ($(hostname)):"
curl -sS -m 20 -o /dev/null -w "  cds.climate.copernicus.eu -> HTTP %{http_code} in %{time_total}s\n" \
     https://cds.climate.copernicus.eu/api/v2 2>&1 | head -3
curl -sS -m 20 -o /dev/null -w "  pypi.org                  -> HTTP %{http_code}\n" \
     https://pypi.org/simple/ 2>&1 | head -3
echo "  proxy vars: http_proxy='${http_proxy:-unset}' https_proxy='${https_proxy:-unset}'"

echo
echo "### from a COMPUTE node (30s test job; Ctrl-C if the queue is slow):"
srun --partition=defq --time=00:01:00 --mem=1G --cpus-per-task=1 \
     bash -c 'echo "  ran on $(hostname)";
              curl -sS -m 15 -o /dev/null -w "  cds -> HTTP %{http_code}\n" https://cds.climate.copernicus.eu/api/v2 || echo "  cds -> NO INTERNET";
              curl -sS -m 15 -o /dev/null -w "  pypi -> HTTP %{http_code}\n" https://pypi.org/simple/ || echo "  pypi -> NO INTERNET"' \
     2>&1 | head -20

section "7. cfgrib / eccodes -- can we even open a GRIB here?"
python3 -c "import cfgrib, xarray; print('  cfgrib', cfgrib.__version__, '| xarray', xarray.__version__)" 2>&1 | head -5

section "DONE"
echo "Bring back: $OUT"
