"""
chunk_stitch.py
================
Q-GLOBAL-2 track: the 1-timestep overlap buffer at every month/year chunk
boundary, needed because F2D's isentropic-A9 material derivative (Q-F2D-5)
silently degrades to a one-sided d/dt at every unbuffered chunk edge.

KEY DESIGN POINT: this needs ZERO extra CDS download requests. Each
month's own file already contains its full calendar range (including its
own first and last timesteps) -- month M-1's LAST timestep and month M+1's
FIRST timestep are already sitting on disk once those months are
downloaded via the existing 1_download.py per-month loop. The "buffer" is
a PROCESSING-TIME concatenation, not a download-time duplication:

    1. To compute month M's diagnostics, open month M's file PLUS the
       single bordering timestep borrowed from month M-1 (its last) and
       month M+1 (its first).
    2. Compute F2D (and anything else needing d/dt) across this extended
       range -- every point in month M's OWN calendar range now has a
       true neighbor on both sides for a centered difference.
    3. Trim the two borrowed timesteps back OUT before month M's output
       enters annual aggregation (Q-AGG-5) -- they belong to month M-1
       and M+1's own official ranges, and double-counting them would
       bias the annual normalization divisor.

CAVEAT / fallback: if raw per-month files get deleted immediately after
processing (a real possibility given Q-MODEL-3-CHECK-A's storage
numbers), the natural "borrow from the still-on-disk neighbor" trick
above doesn't work, because the neighbor may already be gone. In that
case, persist a tiny per-month "boundary cache" (just the first + last
timestep, negligible size) before deleting the full raw file, and borrow
from THAT instead. Not implemented here since it's a fallback for a
storage policy not yet decided -- flagged for whoever finalizes the
delete-after-process policy.

TRUE, UNAVOIDABLE edge cases: the very first timestep of the whole run
(Jan 1979 00Z -- no Dec 1978 data exists in this project's scope) and the
very last (Dec 2020 21Z -- no Jan 2021 data planned) have no bordering
data at all, regardless of this scheme. These 2 timesteps out of 122,728
total (Q-MODEL-3-CHECK-A) get a one-sided d/dt; this is a known,
negligible, and irreducible limitation given the project's date range,
not a bug in this stitching logic.
"""
from __future__ import annotations

from pathlib import Path

import xarray as xr


def _open(path: Path) -> xr.Dataset:
    """Open a month file, dispatching on extension.

    The downloader writes GRIB (`download_plan.month_filename` ->
    'era5_YYYY-MM.grib'), and a bare `xr.open_dataset` on a .grib depends on
    engine auto-detection rather than saying so. Dispatch explicitly, and
    suppress the cfgrib .idx sidecar -- these files are opened read-only to
    borrow a single boundary timestep, so a persisted index is pure clutter
    in the data directory.

    NOTE: `filename_pattern` below defaults to '.grib' to match the
    downloader. It previously defaulted to '.nc', which meant this module
    could never have found a single file the pipeline actually produces.
    """
    if path.suffix in (".grib", ".grb", ".grb2"):
        return xr.open_dataset(path, engine="cfgrib",
                               backend_kwargs={"indexpath": ""})
    return xr.open_dataset(path)


def month_bounds(year: int, month: int) -> tuple[int, int]:
    """Return (prev_year, prev_month) and is-there-a-next-month info via
    simple calendar arithmetic (no external date library needed for this)."""
    prev_year, prev_month = (year - 1, 12) if month == 1 else (year, month - 1)
    return prev_year, prev_month


def next_month(year: int, month: int) -> tuple[int, int]:
    return (year + 1, 1) if month == 12 else (year, month + 1)


def load_month_with_buffer(
    month_dir: Path,
    year: int,
    month: int,
    filename_pattern: str = "era5_{year}-{month:02d}.grib",
    time_dim: str = "time",
) -> tuple[xr.Dataset, slice]:
    """Load month (year, month)'s own file, concatenated with 1 bordering
    timestep from the previous and next month's files IF those files
    exist on disk (true edges of the whole 1979-2020 run naturally have
    no border to borrow -- handled gracefully, not an error).

    Returns (buffered_dataset, own_month_time_slice) -- own_month_time_slice
    identifies exactly which part of the returned, buffered dataset
    corresponds to this month's OWN calendar range, for trimming after
    computing any time-derivative diagnostic (see trim_to_own_month below).
    """
    own_path = month_dir / filename_pattern.format(year=year, month=month)
    own_ds = _open(own_path)
    own_start, own_end = own_ds[time_dim].values[0], own_ds[time_dim].values[-1]

    pieces = []
    prev_y, prev_m = month_bounds(year, month)
    prev_path = month_dir / filename_pattern.format(year=prev_y, month=prev_m)
    if prev_path.exists():
        prev_ds = _open(prev_path)
        pieces.append(prev_ds.isel({time_dim: [-1]}))  # borrow ONLY its last timestep

    pieces.append(own_ds)

    next_y, next_m = next_month(year, month)
    next_path = month_dir / filename_pattern.format(year=next_y, month=next_m)
    if next_path.exists():
        next_ds = _open(next_path)
        pieces.append(next_ds.isel({time_dim: [0]}))  # borrow ONLY its first timestep

    buffered = xr.concat(pieces, dim=time_dim) if len(pieces) > 1 else own_ds
    return buffered, slice(own_start, own_end)


def trim_to_own_month(computed: xr.DataArray, own_time_slice: slice, time_dim: str = "time") -> xr.DataArray:
    """After computing a time-derivative diagnostic (e.g. F2D) on a
    buffered dataset, trim back to exactly this month's own calendar
    range -- the borrowed boundary timesteps must NOT enter annual
    aggregation under this month's label (they belong to the neighbor)."""
    return computed.sel({time_dim: own_time_slice})
