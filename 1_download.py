"""
Filename: data_download.py

Purpose:
    Download ERA5 pressure-level fields for the North Atlantic CAT
    diagnostics pipeline (partial replication of Prosser et al. 2023, GRL).

    Variables downloaded are the *minimal set* needed for the 21
    Williams & Joshi (2013) / Prosser (2023) clear-air turbulence
    diagnostics computed by the rojak pipeline (cat_pipeline.py).

    Variable justification (see docs section below):
      u, v, t, z, d, vo, pv   — required for one or more W&J diagnostics
      q (specific_humidity)   — stubbed with zeros in the pipeline (rojak
                                schema quirk); no need to download
      r (relative_humidity)   — not used by any W&J diagnostic
      w (vertical_velocity)   — not used by any W&J diagnostic

Author:
    Yannick van Etten
"""
from __future__ import annotations
from pathlib import Path
import cdsapi


# ---------------------------------------------------------------------------
# Variable set (7 fields).
# NOTE: rojak.core.data.CATData.required_variables lists
#     {eastward_wind, northward_wind, temperature, geopotential,
#      divergence_of_wind, vorticity, potential_vorticity, specific_humidity}
# Of those, specific_humidity is unused by the 21 W&J diagnostics, so we
# stub it with zeros in cat_pipeline.prepare_for_rojak() rather than
# download it.
# ---------------------------------------------------------------------------
VARIABLES = [
    "u_component_of_wind",   # u   — 18 of 21 diagnostics
    "v_component_of_wind",   # v   — 18 of 21 diagnostics
    "temperature",           # t   —  6 diagnostics: CP, |∇T|, NGM2, -Ri, F2D, NCSU1
    "geopotential",          # z   — 11 diagnostics (any vertical-derivative one)
    "divergence",            # d   —  2 diagnostics: |div|, TI2
    "vorticity",             # vo  —  7 diagnostics: Brown1/2, RVA, UBF, NVA, vort², NCSU1
    "potential_vorticity",   # pv  —  1 diagnostic: |PV|
]

# ---------------------------------------------------------------------------
# Pressure levels: CDS pressure-level product does not expose the paper's
# 188/197/206 hPa (those are model-level). We use 175/200/225 as the closest
# available substitute, with 200 hPa as the target flight level and 175/225
# as the outer levels for vertical derivatives (δp = 50 hPa vs paper's 18 hPa).
# ---------------------------------------------------------------------------
PRESSURE_LEVELS = ["175", "200", "225"]

# ---------------------------------------------------------------------------
# North Atlantic box — CDS "area" is [N, W, S, E].
# ---------------------------------------------------------------------------
AREA = [60, -75, 30, 0]

# 3-hourly instantaneous timestamps (matches Prosser et al. 2023).
TIMES_3H = ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"]


def build_request(year: str, month: str, days: list[str] | str) -> dict:
    """Build a CDS request dict for a given year/month/day set.

    Same variables and levels are used for the single-day trial and
    for any multi-month scaling — this stops the trial and the full
    replication from silently drifting apart.
    """
    if isinstance(days, str):
        days = [days]
    return {
        "product_type": ["reanalysis"],
        "variable": VARIABLES,
        "year": [year],
        "month": [month],
        "day": days,
        "time": TIMES_3H,
        "area": AREA,
        "pressure_level": PRESSURE_LEVELS,
        "data_format": "grib",
        "download_format": "unarchived",
    }


def download_trial(output_dir: Path) -> Path:
    """Download the single-day trial file: 2016-01-01."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "climate_data_01_01_2016.grib"

    request = build_request(year="2016", month="01", days=["01"])
    client = cdsapi.Client()
    print(f"Requesting ERA5 trial → {out_path}")
    print(f"  variables: {len(VARIABLES)} fields ({', '.join(v.split('_')[0] for v in VARIABLES)})")
    print(f"  levels   : {PRESSURE_LEVELS}")
    print(f"  area     : {AREA}  (N, W, S, E)")

    client.retrieve("reanalysis-era5-pressure-levels", request).download(str(out_path))
    print(f"Saved: {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Multi-month loop scaffold — commented out for the trial run.
# Uncomment and call download_month_range() when scaling to 1979–2020.
# ---------------------------------------------------------------------------
def download_month_range(
    output_dir: Path,
    year_start: int = 1979,
    year_end: int = 2020,
) -> None:
    """Loop the same request over every month in [year_start, year_end].

    One CDS request per month keeps under the CDS per-request size cap
    (roughly 120 000 fields).  For 42 years × 12 months = 504 downloads.
    Each is ~500 MB GRIB on this domain, so plan ~250 GB total on disk.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    client = cdsapi.Client()
    for year in range(year_start, year_end + 1):
        for month in range(1, 13):
            month_str = f"{month:02d}"
            out_path = output_dir / f"era5_{year}-{month_str}.grib"
            if out_path.exists():
                print(f"Skipping existing {out_path.name}")
                continue
            # All days of the month; CDS silently drops missing days for shorter months.
            days = [f"{d:02d}" for d in range(1, 32)]
            request = build_request(year=str(year), month=month_str, days=days)
            print(f"Requesting {year}-{month_str} → {out_path.name}")
            client.retrieve("reanalysis-era5-pressure-levels", request).download(str(out_path))


def main():
    # Point OUTPUT_DIR at wherever you want the .grib files to land.
    # On Windows / OneDrive you might use:
    #     OUTPUT_DIR = Path.home() / "OneDrive" / "Documenten" / "Universiteit" / "Turbulence project" / "Data"
    # On Linux/HPC:
    #     OUTPUT_DIR = Path("/scratch/yannick/era5")
    OUTPUT_DIR = Path.home() / "OneDrive" / "Documenten" / "Universiteit" / "Turbulence project" / "Data"

    # --- Trial run (single day) ------------------------------------------
    download_trial(OUTPUT_DIR)

    # --- Full replication (uncomment when ready) -------------------------
    # download_month_range(OUTPUT_DIR / "monthly", year_start=1979, year_end=2020)


if __name__ == "__main__":
    main()
