#!/bin/bash
# Step 1 of the CAT animations: reduce months of the dataset to small frames files.
#     sbatch ada/cat_month_frames.sbatch                               # Feb 2020
#     sbatch --export=ALL,YEAR=2020,MONTH=1 ada/cat_month_frames.sbatch
#     sbatch --export=ALL,YEAR=2020 ada/cat_month_frames.sbatch        # all 12 months of 2020
# Output: $BASE/explore/frames/cat_frames_YYYY-MM.npz (~3-5 MB per month)
# Step 2 (drawing) needs matplotlib, which the pixi env does not have: run
# ada/cat_month_render.py locally on the .npz (see its docstring).
#SBATCH --job-name=cat-frames
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/cat-frames-%j.out
set -uo pipefail
source essentials/jobs/_common.sh
if [ -n "${MONTH:-}" ] || [ -z "${YEAR:-}" ]; then
    pixi run python ada/cat_month_frames.py --year "${YEAR:-2020}" --month "${MONTH:-2}"
else
    pixi run python ada/cat_month_frames.py --year "$YEAR"
fi
