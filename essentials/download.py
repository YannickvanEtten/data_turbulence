"""
download.py -- Stage 1: ERA5 pressure-level data from the Copernicus CDS.

One request = one month (North Atlantic) or one ~10-day block (global).
Each file is:
  * downloaded to <name>.tmp,
  * checked: 7 variables, 3 levels, 8 x days time steps, and the grid size of
    the requested domain (the grid is what tells a global file from a North
    Atlantic one -- the other counts are identical),
  * only then renamed to its final name.
A final file that already exists and passes the check is skipped, so the same
command (or SLURM array) can be re-submitted after any failure.

CDS rejects more than ~4 simultaneous requests per user: run arrays with %4.
The three global day-blocks of a month are joined with merge_blocks(); GRIB is
a sequence of self-contained messages, so joining is plain concatenation.
"""
from __future__ import annotations

import calendar
import os
import shutil
from pathlib import Path

import xarray as xr

import config as C


def build_request(year: int, month: int, domain: str, days=None) -> dict:
    """The CDS request for one month (or the given days of it)."""
    if days is None:
        days = range(1, calendar.monthrange(year, month)[1] + 1)
    return {
        "product_type": ["reanalysis"],
        "variable": list(C.VARIABLES),
        "year": [str(year)],
        "month": [f"{month:02d}"],
        "day": [f"{int(d):02d}" for d in days],
        "time": list(C.TIMES_3H),
        "area": list(C.DOMAINS[domain]["area"]),
        "pressure_level": [str(p) for p in C.PRESSURE_LEVELS_HPA],
        "data_format": "grib",
        "download_format": "unarchived",
    }


def block_days(year: int, month: int, block: int) -> list[int]:
    """Days of global day-block 0, 1 or 2 (1-10, 11-20, 21-end of month)."""
    first, last = C.GLOBAL_DAY_BLOCKS[block]
    return list(range(first, min(last, calendar.monthrange(year, month)[1]) + 1))


def block_path(year: int, month: int, block: int) -> Path:
    days = block_days(year, month, block)
    tag = "-".join(f"{d:02d}" for d in days)
    return C.RAW / "global" / f"era5_glob_{year}-{month:02d}_d{tag}.grib"


def expected_grid(domain: str) -> tuple[int, set]:
    north, west, south, east = C.DOMAINS[domain]["area"]
    n_lat = round((north - south) / C.RESOLUTION_DEG) + 1
    span = (east - west) % 360
    if span == 0:
        n = round(360 / C.RESOLUTION_DEG)
        return n_lat, {n, n + 1}
    return n_lat, {round(span / C.RESOLUTION_DEG) + 1}


def check_file(path: Path, request: dict, domain: str) -> list[str]:
    """Problems with a downloaded file (an empty list means it is fine)."""
    ds = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
    try:
        problems = []
        if len(ds.data_vars) != len(request["variable"]):
            problems.append(f"{len(ds.data_vars)} variables, expected {len(request['variable'])}")
        if ds.sizes.get("isobaricInhPa") != len(request["pressure_level"]):
            problems.append(f"levels {ds.sizes.get('isobaricInhPa')}, expected {len(request['pressure_level'])}")
        n_time = len(request["day"]) * C.STEPS_PER_DAY
        if ds.sizes.get("time") != n_time:
            problems.append(f"time steps {ds.sizes.get('time')}, expected {n_time}")
        n_lat, n_lon = expected_grid(domain)
        if ds.sizes.get("latitude") != n_lat or ds.sizes.get("longitude") not in n_lon:
            problems.append(f"grid {ds.sizes.get('latitude')}x{ds.sizes.get('longitude')}, "
                            f"expected {n_lat}x{sorted(n_lon)} for {domain} -- wrong domain?")
        return problems
    finally:
        ds.close()


def fetch(request: dict, final: Path, domain: str, force: bool = False) -> str:
    """Download `request` to `final` atomically; return 'skip' or 'downloaded'."""
    final = Path(final)
    if final.exists() and not force:
        problems = check_file(final, request, domain)
        if problems:
            raise RuntimeError(f"{final.name} exists but fails the check: {problems}. "
                               f"Not overwriting it; delete it or use --force.")
        return "skip"
    import cdsapi                                  # only needed when really downloading
    client = cdsapi.Client(retry_max=8, sleep_max=60, timeout=120)   # fail within minutes, not hours
    final.parent.mkdir(parents=True, exist_ok=True)
    tmp = final.with_name(final.name + ".tmp")
    client.retrieve(C.DATASET, request, str(tmp))
    problems = check_file(tmp, request, domain)
    if problems:
        raise RuntimeError(f"{tmp.name} fails the check: {problems} (left in place)")
    os.replace(tmp, final)
    return "downloaded"


def download_month(year: int, month: int, force: bool = False) -> str:
    """North Atlantic month -> raw/north_atlantic/era5_na_YYYY-MM.grib"""
    request = build_request(year, month, "north_atlantic")
    return fetch(request, C.raw_month_path("north_atlantic", year, month), "north_atlantic", force)


def download_global_block(year: int, month: int, block: int, force: bool = False) -> str:
    """Global day-block -> raw/global/era5_glob_YYYY-MM_dDD-..-DD.grib"""
    request = build_request(year, month, "global", days=block_days(year, month, block))
    return fetch(request, block_path(year, month, block), "global", force)


def merge_blocks(year: int, month: int) -> str:
    """Join the three global day-blocks into raw/global/era5_glob_YYYY-MM.grib."""
    final = C.raw_month_path("global", year, month)
    if final.exists():
        return "skip"
    blocks = [block_path(year, month, b) for b in range(3)]
    missing = [b.name for b in blocks if not b.exists()]
    if missing:
        raise FileNotFoundError(f"missing blocks {missing}")
    tmp = final.with_name(final.name + ".tmp")
    with open(tmp, "wb") as out:
        for b in blocks:
            with open(b, "rb") as src:
                shutil.copyfileobj(src, out)
    if tmp.stat().st_size != sum(b.stat().st_size for b in blocks):
        raise RuntimeError(f"size mismatch while merging {final.name}")
    os.replace(tmp, final)
    return "merged"
