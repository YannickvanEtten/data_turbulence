"""
download_plan.py — pure request-planning for the ERA5 CAT download layer.

Role in the project: this is a reusable logic module, the same kind of thing
as aggregate.py / trend.py / chunk_stitch.py -- it is imported by the stage
entrypoint (1_download_hpc.py), not run on its own. It does NO network calls
and NO filesystem I/O. Everything here is deterministic and unit-testable
with no mocking, which is the whole reason it lives apart from the downloader.

What it decides:
  - the LOCKED science constants (variables, levels, times),
  - which calendar days a month actually has,
  - the canonical filename every pipeline layer agrees on,
  - how many timesteps a correct month/day MUST contain (for the integrity
    check in 1_download_hpc.py -- not just "the file exists"),
  - the (year, month) list for an arbitrary span,
  - the CDS request dict for one month (or one explicit set of days).

LOCKED vs PARAMETER (deliberate):
  LOCKED constants -- a wrong value must require a visible code edit, never a
  silent CLI flag: VARIABLES, PRESSURE_LEVELS, TIMES_3H.
  PARAMETER -- varies across runs of the same science: date span, output dir
  (both live in the entrypoint), and area. area defaults to the locked North
  Atlantic box; "worldwide" is a ~28x storage/compute jump, a deliberate
  scope decision, not something a flag should make easy to trip into.
"""

import calendar


# ---------------------------------------------------------------------------
# LOCKED constants -- change these only with a visible code edit, never a flag.
# ---------------------------------------------------------------------------
VARIABLES = [
    "u_component_of_wind",   # u   — 18 of 21 diagnostics
    "v_component_of_wind",   # v   — 18 of 21 diagnostics
    "temperature",           # t   — CP, |grad T|, NGM2, -Ri, F2D, NCSU1
    "geopotential",          # z   — vertical-derivative diagnostics
    "divergence",            # d   — |div|, TI2
    "vorticity",             # vo  — Brown1/2, RVA, UBF, NVA, vort², NCSU1
    "potential_vorticity",   # pv  — |PV|
]

PRESSURE_LEVELS = ["175", "200", "225"]

TIMES_3H = ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"]

DATASET = "reanalysis-era5-pressure-levels"

STEPS_PER_DAY = len(TIMES_3H)  # 8 -- used by the integrity check in the entrypoint

# North Atlantic box -- CDS "area" order is [North, West, South, East].
# The validated superset box every diagnostic has been checked against.
NORTH_ATLANTIC_BOX = [60, -75, 30, 0]

# Global box, for Prosser's year-2000 calibration domain.
GLOBAL_BOX = [90, -180, -90, 180]

RESOLUTION = 0.25  # degrees; ERA5 pressure-level product native grid

# ---------------------------------------------------------------------------
# NAMED DOMAINS
# ---------------------------------------------------------------------------
# The project needs TWO domains at once, for two different purposes:
#   north_atlantic -- where trends are computed (Prosser's analysis region)
#   global         -- where severity thresholds are CALIBRATED (year 2000)
# They are deliberately different datasets; see 3_pipeline.run_layers_2_to_6,
# which takes calibration_fields separately from diagnostic_fields.
#
# Naming them (rather than passing raw --area boxes) makes the two downloads
# symmetric commands instead of "one is the default and the other is a flag
# you must remember to get right".
#
# `code` goes into the filename. That is not cosmetic -- see month_filename.
DOMAINS: dict[str, dict] = {
    "north_atlantic": {"code": "na",   "area": NORTH_ATLANTIC_BOX},
    "global":         {"code": "glob", "area": GLOBAL_BOX},
}
DEFAULT_DOMAIN = "north_atlantic"


def domain_area(domain: str) -> list:
    if domain not in DOMAINS:
        raise KeyError(f"unknown domain {domain!r}; known: {sorted(DOMAINS)}")
    return list(DOMAINS[domain]["area"])


def domain_code(domain: str) -> str:
    if domain not in DOMAINS:
        raise KeyError(f"unknown domain {domain!r}; known: {sorted(DOMAINS)}")
    return DOMAINS[domain]["code"]


def expected_grid(area, resolution: float = RESOLUTION):
    """(n_lat, n_lon) a file covering `area` at `resolution` must have.

    This is what makes a wrong-domain file detectable. The variable, level and
    timestep counts are IDENTICAL for a global and a North Atlantic month --
    only the grid differs -- so counts alone cannot tell them apart.

    Longitude is returned as a set of acceptable values because a full 360-deg
    request may or may not repeat the wrap-around meridian depending on how
    CDS renders it; both 1440 and 1441 are legitimate for a global grid.
    """
    north, west, south, east = area
    n_lat = int(round((north - south) / resolution)) + 1
    span_lon = (east - west) % 360
    if span_lon == 0:  # full circle
        n_lon_full = int(round(360.0 / resolution))
        return n_lat, {n_lon_full, n_lon_full + 1}
    n_lon = int(round(span_lon / resolution)) + 1
    return n_lat, {n_lon}


def month_days(year, month):
    """Exact day strings ('01', '02', ...) for this calendar month.

    Uses calendar.monthrange instead of a fixed range(1, 32) -- we want the
    real day count, not "ask CDS for 31 and hope it silently drops the rest."
    """
    n = calendar.monthrange(year, month)[1]
    return [f"{d:02d}" for d in range(1, n + 1)]


def month_filename(year, month, domain=DEFAULT_DOMAIN):
    """Canonical filename for one month -- the one name every layer of the
    pipeline (download, chunk_stitch, diagnostics) agrees on.

    THE DOMAIN CODE IS PART OF THE NAME, and must stay that way. Without it,
    'era5_2000-01.grib' means either the North Atlantic January 2000 or the
    global January 2000, and downloading one would silently overwrite the
    other. That collision is not caught by the integrity check either: both
    files have exactly 7 variables, 3 levels and 8*31 timesteps. Only the grid
    size differs, which is why expected_grid() exists and why verify_file()
    now checks it.
    """
    return f"era5_{domain_code(domain)}_{year}-{month:02d}.grib"


def expected_timesteps(n_days):
    """Timesteps a correct file MUST contain, given how many days it covers.
    Kept general (takes a day count) so it works for a full month or a single
    trial day, not just full months."""
    return n_days * STEPS_PER_DAY


def month_expected_timesteps(year, month):
    """Convenience wrapper: expected timesteps for a full calendar month."""
    return expected_timesteps(len(month_days(year, month)))


def month_span(start_ym, end_ym):
    """Inclusive list of (year, month) tuples from 'YYYY-MM' to 'YYYY-MM'.

    Fails loud on a malformed string, an out-of-range month, or a reversed
    span, rather than silently returning an empty or backwards range.
    """
    def parse(ym):
        parts = ym.split("-")
        if len(parts) != 2:
            raise ValueError(f"expected 'YYYY-MM', got {ym!r}")
        year, month = int(parts[0]), int(parts[1])
        if month < 1 or month > 12:
            raise ValueError(f"month out of range in {ym!r}")
        return year, month

    start_year, start_month = parse(start_ym)
    end_year, end_month = parse(end_ym)
    if (start_year, start_month) > (end_year, end_month):
        raise ValueError(f"start {start_ym} is after end {end_ym}")

    months = []
    year, month = start_year, start_month
    while (year, month) <= (end_year, end_month):
        months.append((year, month))
        if month == 12:
            year, month = year + 1, 1
        else:
            month = month + 1
    return months


def build_request(year, month, area=None, days=None, domain=DEFAULT_DOMAIN):
    """Build the CDS request dict for one month, or for an explicit set of days.

    days=None means the full calendar month. Passing days (e.g. ["01"]) is how
    the single-day trial is built -- same request shape, one day. The list
    values are copied so a caller that mutates request["variable"] can't
    silently corrupt the module constant for everyone else.

    area defaults to the locked NORTH_ATLANTIC_BOX; a different box has to be
    passed explicitly, never silently defaulted.
    """
    if area is None:
        area = domain_area(domain)
    if days is None:
        days = month_days(year, month)
    return {
        "product_type": ["reanalysis"],
        "variable": list(VARIABLES),
        "year": [str(year)],
        "month": [f"{month:02d}"],
        "day": list(days),
        "time": list(TIMES_3H),
        "area": list(area),
        "pressure_level": list(PRESSURE_LEVELS),
        "data_format": "grib",
        "download_format": "unarchived",
    }