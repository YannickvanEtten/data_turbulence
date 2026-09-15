#!/bin/bash
# Stage 1c: join the three day-blocks of each global month (plain concatenation).
#SBATCH --job-name=ess-merge
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ess-merge-%j.out
set -uo pipefail
source essentials/jobs/_common.sh
run merge-global
