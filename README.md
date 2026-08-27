# Clear-Air Turbulence Pipeline — Prosser et al. (2023) partial replication

Five files, one job each. The number prefixes are the workflow order.

```
1_download.py     Download ERA5 pressure-level data from CDS
2_diagnostics.py  THE diagnostics: data prep + all 21 (14 rojak + 7 hand-written)
3_pipeline.py     Read a GRIB -> compute 21 -> save diagnostics.nc (+ stats, comparison)
4_verify.py       One-time correctness check on the 7 hand-written diagnostics
5_explore.ipynb   Visual exploration of one day's output
archive/          Superseded scripts, kept for reference (not used by anything)
```

## The workflow

```
1_download.py  ──►  climate_data_*.grib
                         │
2_diagnostics.py ◄───────┤  (imported by 3 and 4, never run directly)
                         │
3_pipeline.py  ──────────►  cat_outputs/diagnostics.nc   (21 variables)
                                    │                     cat_outputs/diagnostics.zarr
                                    │                     cat_outputs/summary_stats.csv
                                    │                     cat_outputs/comparison_table.csv
                         ┌──────────┘
5_explore.ipynb  ────────►  plots, animation, correlation matrix
```

`4_verify.py` is run once against a small file to confirm the 7 hand-written
diagnostics are correct; it is not part of the production run.

## How to run

```bash
# 1. one-time: download a trial day (edit the output path inside the file first)
python 1_download.py

# 2. compute all 21 diagnostics for that file
python 3_pipeline.py "path/to/climate_data_02_07_2016.grib" ./cat_outputs

# 3. (optional, one-time) confirm the hand-written diagnostics are correct
python 4_verify.py path/to/era5_validation_subset.nc

# 4. explore visually
jupyter lab 5_explore.ipynb
```

## Why 2_diagnostics.py splits rojak vs hand-written

A Phase-A code review of rojak (@ commit 1a65326) against Sharman et al. (2006)
Appendix A found that 14 of the 21 diagnostics are correct and 7 are not:

| From rojak (14, verified PASS) | Hand-written (7, rojak buggy or absent) |
|---|---|
| PV, Brown1, \|∇T\|, divergence, VWS, Endlich, DEF (√), wind speed, NGM2, NVA, ζ², TI1, NGM1, TI2 | Richardson, Colson–Panofsky, UBF, Frontogenesis2D, NCSU1, RVA, Brown2 |

The 7 hand-written ones live as functions in `2_diagnostics.py`, each citing its
Sharman equation. `4_verify.py` cross-checks every one against a temporarily
corrected rojak and confirms agreement (Spearman ρ ≥ 0.99, two of them exact).

The root bug behind three of the seven was the Richardson number: rojak computes
N²/Sv (units s⁻¹) instead of N²/Sv² (dimensionless), which also silently broke
Colson–Panofsky and NCSU1. Fixed here.

## Editing diagnostics

All diagnostic logic is in `2_diagnostics.py`. To change a formula, edit it there;
both `3_pipeline.py` and `4_verify.py` import from it, so they stay in sync
automatically. `3_pipeline.py` is intentionally thin — only load/compute/save.

## Scaling to the full 1979–2020 run

`1_download.py` already contains a `download_month_range()` scaffold. Once the
trial checks out, loop it over the full period, then run `3_pipeline.py` per
file (or adapt `run()` to `open_mfdataset` for Dask-parallel multi-file input).
