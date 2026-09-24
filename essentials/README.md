# essentials — the CAT indicator pipeline, short and readable

**Start with [`showcase.html`](showcase.html)** — open it in a browser and read the
whole pipeline end to end: the data, every formula next to the code that runs it,
the maps it produces, and the 42-year trend. `showcase.ipynb` is the same thing
live: it runs on **one day** of ERA5 in about a minute and re-computes nothing
else.

Behind it is a compact re-implementation of the data pipeline in this repository.
It produces **the same 21 clear-air-turbulence diagnostics, thresholds and
trends** as the original code (verified bit for bit, see below), in six plain
Python files with no classes and no rojak dependency. It exists to make the
dataset inspectable, not to rebuild it: the 504 months on the share stay as they
are.

What it does: ERA5 → 21 CAT diagnostics at 200 hPa → severity thresholds from
the global year 2000 → exceedance frequencies over Prosser's North Atlantic
box → 1979–2020 trends, compared with Prosser et al. (2023) Table 1.

**Why every choice was made: [DECISIONS.md](DECISIONS.md).**

## Files

```
showcase.ipynb   THE SHOWCASE: the whole pipeline on one day, formulas + code + figures
showcase.html    the same notebook already executed — nothing to install to read it
config.py        every constant, path and choice (incl. the two CONVENTIONS)
download.py      Stage 1  ERA5 from CDS: request, integrity check, atomic write
diagnostics.py   Stage 2  geometry, derivative operators, the 21 formulas, month driver
thresholds.py    Stage 3  cos(lat)-weighted Hazen percentiles of the global year 2000
trends.py        Stage 4  exceedance -> seasonal frequencies -> OLS trends vs Prosser
main.py          the single command-line entry point
jobs/            SLURM scripts for ADA (one per stage, plus two checks)
tests/           test_operators.py (standalone), compare_*_with_original.py (vs the old code)
DECISIONS.md     what, why, source, and where in the code
```

Read `diagnostics.py` top to bottom: sections 2–3 are the maths every
diagnostic uses (grid spacing, map-scale factors, derivatives of scalars and of
wind components, vertical derivatives), section 5 is one short function per
diagnostic with its formula and equation number.

## The pipeline

```
                     CDS
                      │  download.py  (%4 concurrent requests)
          ┌───────────┴────────────┐
 raw/north_atlantic/          raw/global/
 era5_na_YYYY-MM.grib         era5_glob_2000-MM.grib   (3 day-blocks joined)
 504 months                   12 months
          │ diagnostics.py         │ diagnostics.py
          ▼                        ▼
 essentials/derived/          essentials/derived/
 north_atlantic_<conv>/       global_<conv>/
 21 × (lat, lon, time), 200 hPa, float32 zarr
          │                        │ thresholds.py
          │                        ▼
          │               essentials/calibration/thresholds_<conv>_<date>.json
          │  trends.py  ◄──────────┘
          ▼
 essentials/results/frequencies_<conv>_1979-2020.csv   (per year, season, severity, diagnostic)
 essentials/results/trends_<conv>_1979-2020.csv        (OLS, 95 % CI, Prosser comparison)
```

`<conv>` is `baseline` or `audit` — the two sets of implementation choices that
the existing results were produced under. They differ in three switches (the
north–south grid spacing, whether vorticity/divergence/PV are ERA5's archived
fields or computed from u, v, T, z, and the sign reading of the frontogenesis
function). The choice between them is still open; DECISIONS.md §3 sets out what
each is, the case for each, and what it changes (trends: nothing material;
levels: audit is closer to Prosser).

Nothing is ever written into the original `derived/` or `calibration/`
folders: all output goes under `$BASE/essentials/`.

## The showcase

`showcase.ipynb` imports these files rather than copying them, so it cannot drift
from the pipeline: the formulas it prints are the running functions. It covers

1. the ERA5 input (what is downloaded and why),
2. the grid — degrees to metres, and where the two conventions first differ,
3. the four derivative operators, with the curvature term demonstrated,
4. the 21 diagnostics: formula, source equation, the code, a map and its statistics,
5. the two conventions measured against each other on that day,
6. how the severity thresholds are calibrated (and the streaming trick, checked against a full sort),
7. exceedance → frequency → the 42-year trend, with the real figures read from `data_prosser/`,
8. what was verified and how.

To run it: `jupyter lab showcase.ipynb` from this folder. It needs `numpy`,
`xarray`, `pyproj`, `matplotlib` and `cfgrib` (for the GRIB day; without cfgrib it
falls back to `era5_validation_subset.nc` in the repository root). It looks for the
day file at `../../../Universiteit/Turbulence project/Data/climate_data_01_01_2016.grib`
or in `CAT_DEMO_FILE`; §6 and §7 use the real thresholds if a
`calibration/thresholds_*.json` is copied next to the repo, and say so plainly when
they fall back.

## Running on ADA

From the repository root (`$BASE/data_turbulence`), in the existing pixi
environment (numpy, xarray, zarr, cfgrib, pyproj, cdsapi are all in it):

```bash
module load 2025 && module load pixi
mkdir -p logs
```

**First, check it against what is already on the share** (overwrites nothing):

```bash
sbatch essentials/jobs/0_check_against_original.sbatch
sbatch --export=ALL,CONVENTION=audit essentials/jobs/0b_check_thresholds_and_trends.sbatch
```

The first should print `ALL 21 BIT-IDENTICAL` for each convention and
`ALL 21 IDENTICAL` for the 1979-01 stores it rewrites and compares with the
original ones. The second re-calibrates the original `derived/global_audit`
stores (should match `thresholds_2026-09-11.json`) and recomputes the 42-year
trends from the original `derived/north_atlantic_audit` stores (should match
STATUS §18.15: 25/25 significant, 25/25 inside the interval).

**Only if the dataset ever has to be rebuilt** (per convention; submit each array
bare so `%N` is kept). This is not the intended use — the existing stores are the
dataset, and the checks above show this code reproduces them:

```bash
sbatch essentials/jobs/1a_download_na.sbatch            # skips the 504 months already on disk
sbatch essentials/jobs/1b_download_global.sbatch        # skips the 36 blocks already on disk
sbatch essentials/jobs/1c_merge_global.sbatch
sbatch --export=ALL,CONVENTION=audit essentials/jobs/2a_diagnostics_na.sbatch       # ~7 h, ~300 GB
sbatch --export=ALL,CONVENTION=audit essentials/jobs/2b_diagnostics_global.sbatch   # ~14 h, ~200 GB
sbatch --export=ALL,CONVENTION=audit essentials/jobs/3_thresholds.sbatch            # ~2 h
sbatch --export=ALL,CONVENTION=audit essentials/jobs/4_trends.sbatch
```

Storage: a full re-derive is ~500 GB per convention on top of the ~1.5 TB already
used, so in practice point the later stages at the stores that already exist:

```bash
pixi run python essentials/main.py trends --convention audit \
    --derived $BASE/derived/north_atlantic_audit \
    --thresholds $BASE/calibration/thresholds_2026-09-11.json
```

## Running locally on a small file

```bash
python essentials/main.py demo --input climate_data_01_01_2016.grib --convention audit \
       --thresholds thresholds_2026-09-11.json
```

prints the median and 99th percentile of each diagnostic and, with a thresholds
file, the exceedance frequency at each severity — the same thing the showcase
notebook does, without the figures.

## Tests

```bash
python -m pytest essentials/tests/test_operators.py            # standalone, seconds
python essentials/tests/compare_with_original.py <file.grib>    # needs rojak (pixi env)
python essentials/tests/compare_stages_3_4_with_original.py     # needs the original repo
```

## Regenerating the showcase

```bash
jupyter nbconvert --to notebook --execute --inplace showcase.ipynb
jupyter nbconvert --to html showcase.html showcase.ipynb
```

(The HTML renders its formulas with MathJax from a CDN, so read it online, or
export to PDF from the browser if it has to travel offline.)
