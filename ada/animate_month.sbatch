#!/bin/bash
# One month of CAT as an animation + "moving or popping?" statistics.
#     sbatch ada/animate_month.sbatch                                      # Feb 2020, moderate
#     sbatch --export=ALL,YEAR=2020,MONTH=1,SEVERITY=light ada/animate_month.sbatch
# Output: $BASE/explore/animations/cat_<YYYY-MM>_<severity>.{gif,mp4}
#         and cat_<YYYY-MM>_<severity>_moving_or_popping.{png,json}
#SBATCH --job-name=cat-anim
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/cat-anim-%j.out
set -uo pipefail
source essentials/jobs/_common.sh
pixi run python ada/animate_month.py --year "${YEAR:-2020}" --month "${MONTH:-2}" --severity "${SEVERITY:-moderate}"
