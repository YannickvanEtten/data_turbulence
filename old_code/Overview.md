# CAT Diagnostics Trial Pipeline

Partial replication of **Prosser et al. (2023), GRL** — North Atlantic clear-air turbulence diagnostics on ERA5 data, using the `rojak` package from Imperial College London.

## Install

```bash
# Real rojak — NOT `pip install rojak` (that is a different package)
pip install git+https://github.com/ImperialCollegeLondon/rojak.git

# Dependencies (most pulled in by rojak, but be explicit)
pip install xarray netcdf4 zarr cfgrib pandas numpy
```

## Run

```bash
# Real ERA5 trial data:
python cat_pipeline.py climate_data_01_07_2016_5.grib ./cat_outputs

# Or on a NetCDF you've already downloaded:
python cat_pipeline.py era5.nc ./cat_outputs
```

Output:

```
cat_outputs/
├── diagnostics.nc          # all 21 diagnostics, 200 hPa, 8 time steps
├── diagnostics.zarr        # same, Zarr v2 store
├── summary_stats.csv       # native units
└── comparison_table.csv    # converted to Williams & Joshi (2013) Table 1 units
```

## Files in this directory

| File | Purpose |
|---|---|
| `cat_pipeline.py` | The pipeline. End-to-end loader → rojak diagnostics → RVA → stats → CSV/NC/Zarr. |
| `make_synthetic_era5.py` | Generates a synthetic ERA5-like dataset for offline testing. |
| `cat_outputs/VALIDATION_REPORT.md` | Detailed per-diagnostic verdicts and recommended next steps. |

## What the pipeline does (in one paragraph)

Loads a GRIB or NetCDF file with ERA5 pressure-level fields, renames variables and coordinates to match what rojak expects (`isobaricInhPa` → `pressure_level`, short names → CF names), wraps the dataset in `rojak.core.data.CATData`, runs 20 of the 21 Williams & Joshi (2013) diagnostics through `rojak.turbulence.diagnostic.DiagnosticFactory`, adds the 21st diagnostic (relative vorticity advection magnitude) manually using rojak's geospatial derivative helpers, saves outputs to NetCDF and Zarr, and prints a comparison table flagging any diagnostic whose magnitude or sign looks implausible.

## Bugs from the prior code that this pipeline avoids

- Never passes `pint.Quantity` into `xarray.sel()` — always plain integers.
- Directional shear uses rojak's `Endlich` class (geometric vertical derivative from geopotential), not `differentiate('isobaricInhPa') / δz` (wrong units).
- `NEGATIVE_RICHARDSON` is the explicit enum member, removing sign-name confusion.
- For RVA, derivatives are computed with metres-based `grid_spacing` from `rojak.core.derivatives`, never degree-based `np.gradient`.

## Known limitations of the trial

- **One day of data** → cannot reproduce Prosser's trends. Only the *pipeline* is validated here.
- **175/200/225 hPa is an approximation** of the paper's 188/197/206 hPa (those are model-level only; CDS pressure-level product does not expose them). Expect vertical-derivative diagnostics to be systematically damped by ~2.8×.
- Two diagnostics (#15 Brown EDR, #20 Frontogenesis 2D) need verification regardless of input data — see VALIDATION_REPORT.md §6 for details.
