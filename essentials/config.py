"""
config.py -- every fixed choice of the CAT pipeline, in one place.

Nothing in here computes anything. If a number or a path is a decision, it
lives here with a one-line reason, and DECISIONS.md explains it at length.
The four stage files (download.py, diagnostics.py, thresholds.py, trends.py)
import from this file and never hard-code these values themselves.
"""
from pathlib import Path

# ---------------------------------------------------------------------------
# Where things live on ADA
# ---------------------------------------------------------------------------
# BASE is the project share. Raw ERA5 is READ from BASE/raw (it is shared with
# the original pipeline and is never modified). Everything this code WRITES
# goes under OUT_ROOT, so running it can never overwrite the stores and
# threshold files that the results in STATUS.md were computed from.
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
RAW = BASE / "raw"
OUT_ROOT = BASE / "essentials"

# The original pipeline's outputs, used only as READ-ONLY comparison targets
# (main.py compare-store / compare-thresholds).
OLD_DERIVED = {
    ("north_atlantic", "baseline"): BASE / "derived" / "north_atlantic",
    ("north_atlantic", "audit"):    BASE / "derived" / "north_atlantic_audit",
    ("global", "baseline"):         BASE / "derived" / "global",
    ("global", "audit"):            BASE / "derived" / "global_audit",
}
OLD_THRESHOLDS = {
    "baseline": BASE / "calibration" / "thresholds_2026-09-07.json",
    "audit":    BASE / "calibration" / "thresholds_2026-09-11.json",
}

# ---------------------------------------------------------------------------
# ERA5 request (Stage 1). Changing any of these means a new dataset.
# ---------------------------------------------------------------------------
DATASET = "reanalysis-era5-pressure-levels"
VARIABLES = [
    "u_component_of_wind",   # u
    "v_component_of_wind",   # v
    "temperature",           # t
    "geopotential",          # z
    "divergence",            # d   archived; only used by the 'archived' field source
    "vorticity",             # vo  archived; only used by the 'archived' field source
    "potential_vorticity",   # pv  archived; only used by the 'archived' field source
]
PRESSURE_LEVELS_HPA = [175, 200, 225]  # finest CDS stencil around 200 hPa (STATUS 11.6)
TARGET_LEVEL_HPA = 200                 # where every diagnostic is evaluated
TIMES_3H = ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"]
STEPS_PER_DAY = 8
RESOLUTION_DEG = 0.25

# CDS "area" is [North, West, South, East].
DOMAINS = {
    "north_atlantic": {"code": "na",   "area": [60, -75, 30, 0]},    # trend region (superset of Prosser's box)
    "global":         {"code": "glob", "area": [90, -180, -90, 180]}, # calibration region (year 2000)
}
TREND_YEARS = (1979, 2020)
CALIBRATION_YEAR = 2000

# The global months are downloaded in three day-blocks each and then
# concatenated (a full global month in one request was never tested on CDS).
GLOBAL_DAY_BLOCKS = [(1, 10), (11, 20), (21, 31)]

# ---------------------------------------------------------------------------
# Physical and geometric constants (identical to rojak @ 25b8685)
# ---------------------------------------------------------------------------
G = 9.80665                    # gravitational acceleration [m s-2]
OMEGA = 7.292115e-05           # Earth's angular velocity [s-1]
R_EARTH = 6371008.7714         # mean Earth radius [m]  (only used in beta = df/dy)
A_WGS84 = 6378137.0            # WGS84 equatorial radius [m] (map distances)
KAPPA = 0.28571428571428564    # R_d / c_p = 2/7
P0_HPA = 1000                  # reference pressure for potential temperature
RI_CRIT = 0.5                  # Colson-Panofsky critical Richardson number
RI_FLOOR = 1e-5                # NCSU1: Ri is floored at this value (Sharman A36)

# ---------------------------------------------------------------------------
# THE CONVENTIONS -- the three implementation choices that differ between
# the two production series. See DECISIONS.md section 3 for the full story.
# ---------------------------------------------------------------------------
#   dy               How the north-south grid spacing is measured before the
#                    map-scale factor a/M(phi) is applied.
#                      "geodesic": true meridian arc M(phi)*dphi  (rojak/MetPy;
#                                  applies the ellipsoid correction twice)
#                      "map":      map distance a*dphi            (fix F1)
#   fields           Where relative vorticity, divergence and PV come from.
#                      "archived": ERA5's own vo, d, pv fields
#                      "computed": rebuilt from u, v, T, z like Prosser (F12)
#   f2d_variant      Which reading of Sharman (2006) A9 is used for #20.
#                      "C": |0.5 D/Dt Q|   (magnitude)
#                      "A": +0.5 D/Dt Q    (signed; Williams & Storer 2022 Eq. 3)
#   endlich          How the vertical change of wind direction is taken (#14).
#                      "angle":      difference of the direction angle, wrapped
#                                    to the shorter arc (rojak; used in BOTH series)
#                      "components": closed form |u v_z - v u_z| / |V|^2 (F10,
#                                    an open sensitivity, not used in either series)
CONVENTIONS = {
    "baseline": {"dy": "geodesic", "fields": "archived", "f2d_variant": "C", "endlich": "angle"},
    "audit":    {"dy": "map",      "fields": "computed", "f2d_variant": "A", "endlich": "angle"},
}
DEFAULT_CONVENTION = "audit"

# ---------------------------------------------------------------------------
# Severity thresholds (Stage 3) -- Williams (2017) Table 1, used by Prosser
# ---------------------------------------------------------------------------
SEVERITIES = {                 # name -> percentile of the global year-2000 distribution
    "light":              97.0,
    "light_to_moderate":  99.1,
    "moderate":           99.6,
    "moderate_to_severe": 99.8,
    "severe":             99.9,
}
SEVERITY_LABEL = {"light": "LOG", "light_to_moderate": "LMOG", "moderate": "MOG",
                  "moderate_to_severe": "MSOG", "severe": "SOG"}
TAIL_CUT_PERCENTILE = 88.0     # values above p88 are kept; everything below enters via two sums
TAIL_CUT_STRIDE = 7            # every 7th timestep to place the cut (never a multiple of 8)
TAIL_TIME_CHUNK = 64           # timesteps per block in the streaming pass (sets memory)
TAIL_GUARD_MARGIN = 0.005      # the kept tail must start this far below p97

# ---------------------------------------------------------------------------
# Exceedance and trends (Stage 4) -- Prosser (2023)
# ---------------------------------------------------------------------------
PROSSER_BOX = {"lat": (36.0, 60.0), "lon": (-55.0, -10.0)}   # his Figure 2 caption
SEASON_MONTHS = {
    "djf": (12, 1, 2),          # Jan, Feb and Dec of the SAME calendar year
    "mam": (3, 4, 5),
    "jja": (6, 7, 8),
    "son": (9, 10, 11),
    "annual": tuple(range(1, 13)),
}

# Prosser (2023) Table 1: hours per season at an average point in his box,
# 1979 -> 2020, and his relative increase.
PROSSER_TABLE1 = {
    "djf":    {"light": (128.9, 155.6, 0.21), "light_to_moderate": (45.6, 59.3, 0.30),
               "moderate": (22.3, 30.6, 0.37), "moderate_to_severe": (12.1, 17.2, 0.43),
               "severe": (6.4, 9.6, 0.49)},
    "mam":    {"light": (90.4, 113.4, 0.26), "light_to_moderate": (27.2, 38.9, 0.43),
               "moderate": (11.8, 18.6, 0.57), "moderate_to_severe": (5.7, 9.7, 0.71),
               "severe": (2.7, 5.0, 0.85)},
    "jja":    {"light": (114.1, 124.5, 0.09), "light_to_moderate": (36.5, 43.8, 0.20),
               "moderate": (16.1, 21.1, 0.31), "moderate_to_severe": (7.7, 10.9, 0.41),
               "severe": (3.6, 5.5, 0.52)},
    "son":    {"light": (133.1, 153.2, 0.15), "light_to_moderate": (43.4, 53.4, 0.23),
               "moderate": (19.8, 25.8, 0.31), "moderate_to_severe": (10.0, 13.9, 0.39),
               "severe": (5.0, 7.4, 0.47)},
    "annual": {"light": (466.5, 546.8, 0.17), "light_to_moderate": (152.7, 195.4, 0.28),
               "moderate": (70.0, 96.1, 0.37), "moderate_to_severe": (35.5, 51.8, 0.46),
               "severe": (17.7, 27.4, 0.55)},
}

# ---------------------------------------------------------------------------
# The 21 diagnostics, in Williams & Joshi (2013) Table 1 order.
# sign "+" means "more turbulent = larger value" for all 21 (verified against
# Williams 2017 Table 2: every severity ladder increases from light to severe).
# ---------------------------------------------------------------------------
DIAGNOSTICS = [
    # key                     W&J #  name
    ("magnitude_pv",            1,  "Magnitude of potential vorticity |PV|"),
    ("colson_panofsky",         2,  "Colson-Panofsky index"),
    ("brown1",                  3,  "Brown index (Phi)"),
    ("temperature_gradient",    4,  "|Horizontal temperature gradient|"),
    ("horizontal_divergence",   5,  "|Horizontal divergence|"),
    ("vertical_wind_shear",     6,  "Vertical wind shear S_v"),
    ("endlich",                 7,  "Wind speed x directional shear (Endlich)"),
    ("deformation",             8,  "Total deformation DEF"),
    ("wind_speed",              9,  "Wind speed |V|"),
    ("ngm2",                   10,  "|dT/dz| x DEF (NGM2)"),
    ("negative_richardson",    11,  "Negative Richardson number -Ri"),
    ("rva_magnitude",          12,  "|Relative vorticity advection|"),
    ("ubf",                    13,  "|Residual of the nonlinear balance equation| (UBF)"),
    ("nva",                    14,  "Negative absolute vorticity advection"),
    ("brown2",                 15,  "Brown energy dissipation rate"),
    ("vorticity_squared",      16,  "Relative vorticity squared"),
    ("ti1",                    17,  "Ellrod TI1 = S_v x DEF"),
    ("ngm1",                   18,  "|V| x DEF (NGM1)"),
    ("ti2",                    19,  "Ellrod TI2 = S_v x (DEF - div)"),
    ("f2d",                    20,  "Isentropic frontogenesis (Sharman A9)"),
    ("ncsu1",                  21,  "NCSU1 index"),
]
DIAGNOSTIC_KEYS = [k for k, _, _ in DIAGNOSTICS]
SIGN = {k: "+" for k in DIAGNOSTIC_KEYS}


def month_store_name(domain: str, year: int, month: int) -> str:
    """File name of one month of diagnostics. Same pattern as the original."""
    return f"diagnostics_{DOMAINS[domain]['code']}_{year}-{month:02d}.zarr"


def derived_dir(domain: str, convention: str) -> Path:
    """Where THIS code writes diagnostics for a domain and convention."""
    return OUT_ROOT / "derived" / f"{domain}_{convention}"


def raw_month_path(domain: str, year: int, month: int) -> Path:
    """Where the raw ERA5 month lives (shared with the original pipeline)."""
    return RAW / domain / f"era5_{DOMAINS[domain]['code']}_{year}-{month:02d}.grib"
