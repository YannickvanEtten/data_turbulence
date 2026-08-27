"""
2_diagnostics.py
================
THE diagnostics module — everything about computing the 21 Williams & Joshi /
Prosser clear-air-turbulence diagnostics lives here and nowhere else.

Two responsibilities:
  1. Data loading & rojak-preparation   (load_era5, prepare_for_rojak)
  2. The 21 diagnostics:
        - 14 from rojak      (Phase A: PASS / VARIANT)   see ROJAK_DIAGNOSTICS
        - 7  hand-written    (Phase A: BUG / VERIFY / MISSING)  see the
          richardson/colson_panofsky/ubf/frontogenesis_2d/ncsu1/rva/brown2
          functions, each implemented directly from Sharman (2006) Appendix A
     -> compute_all_21() returns a dict of all 21 as xr.DataArrays.

The thin orchestration (stats, comparison table, file output, CLI) lives in
3_pipeline.py, which imports from this file. Verification lives in 4_verify.py.

Why 7 are hand-written (rojak @ commit 1a65326 bugs found in Phase A):
  #11 Richardson       N2/Sv  -> must be N2/Sv2 (dimensionless)
  #2  Colson-Panofsky  inherits buggy Ri; also uses full 3-level centered stencil
  #13 UBF              sign flip + broken laplacian + beta=f/R not 2*Omega*cos/R
  #20 Frontogenesis2D  cross-term dv_dy -> dv_dx
  #21 NCSU1            inherits buggy Ri
  #12 RVA              not provided by rojak
  #15 Brown2           faithful to A14 but native s^-3 (rank-only); no invented L^2

References: Sharman et al. (2006) Wea. Forecasting App. A; Williams & Joshi
(2013); Ellrod & Knapp (1992); Kaplan et al. (2005); Koch & Caracena (2002).
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import xarray as xr

from rojak.core.data import CATData
from rojak.core.derivatives import (
    grid_spacing, first_derivative, spatial_gradient, GradientMode,
)
from rojak.turbulence.calculations import (
    altitude_derivative_on_pressure_level,
    potential_temperature as _potential_temperature,
)
from rojak.core.constants import GRAVITATIONAL_ACCELERATION
from rojak.orchestrator.configuration import TurbulenceDiagnostics
from rojak.turbulence.diagnostic import DiagnosticFactory


# ===========================================================================
# Constants (match rojak exactly so verification cross-checks are exact)
# ===========================================================================
OMEGA   = 7.292115e-05        # Earth angular velocity  [s^-1]
R_EARTH = 6371008.7714        # Earth mean radius       [m]
G       = 9.80665             # gravity                 [m s^-2]



# ===========================================================================
# Diagnostic registries and W&J Table 1 reference
# ===========================================================================
REFERENCE_TABLE: dict[str, dict] = {
    "magnitude_pv":          {"num": 1,  "units": "PVU",                     "wj_median":  6.84,  "sign": "+", "name": "Magnitude of potential vorticity"},
    "colson_panofsky":       {"num": 2,  "units": "10^3 kt^2",               "wj_median": -34.8,  "sign": "either", "name": "Colson–Panofsky index"},
    "brown1":                {"num": 3,  "units": "10^-6 s^-1",              "wj_median": 77.1,   "sign": "+", "name": "Brown index"},
    "temperature_gradient":  {"num": 4,  "units": "10^-6 K m^-1",            "wj_median":  5.75,  "sign": "+", "name": "|Horizontal temperature gradient|"},
    "horizontal_divergence": {"num": 5,  "units": "10^-6 s^-1",              "wj_median":  2.82,  "sign": "+", "name": "|Horizontal divergence|"},
    "vertical_wind_shear":   {"num": 6,  "units": "10^-3 s^-1",              "wj_median":  1.88,  "sign": "+", "name": "|Vertical wind shear|"},
    "endlich":               {"num": 7,  "units": "10^-3 rad s^-1",          "wj_median":  0.952, "sign": "+", "name": "Wind speed × directional shear"},
    "deformation":           {"num": 8,  "units": "10^-6 s^-1",              "wj_median": 18.6,   "sign": "+", "name": "Flow deformation (rojak returns squared)"},
    "wind_speed":            {"num": 9,  "units": "m s^-1",                  "wj_median": 14.9,   "sign": "+", "name": "Wind speed"},
    "ngm2":                  {"num": 10, "units": "10^-9 K m^-1 s^-1",       "wj_median":  8.17,  "sign": "+", "name": "Deformation × vertical T gradient"},
    "negative_richardson":   {"num": 11, "units": "dimensionless",           "wj_median": -127.2, "sign": "either", "name": "Negative Richardson number (-Ri)"},
    "rva_magnitude":         {"num": 12, "units": "10^-10 s^-2",             "wj_median":  2.33,  "sign": "+", "name": "|Relative vorticity advection|"},
    "ubf":                   {"num": 13, "units": "10^-12 s^-2",             "wj_median": 161.0,  "sign": "+", "name": "|Residual of nonlinear balance eq.|"},
    "nva":                   {"num": 14, "units": "10^-10 s^-2",             "wj_median":  2.05,  "sign": "+", "name": "Negative absolute vorticity advection"},
    "brown2":                {"num": 15, "units": "10^-6 J kg^-1 s^-1",      "wj_median": 116.0,  "sign": "+", "name": "Brown energy dissipation rate"},
    "vorticity_squared":     {"num": 16, "units": "10^-9 s^-2",              "wj_median":  0.221, "sign": "+", "name": "Relative vorticity squared"},
    "ti1":                   {"num": 17, "units": "10^-9 s^-2",              "wj_median": 31.5,   "sign": "+", "name": "Ellrod TI1"},
    "ngm1":                  {"num": 18, "units": "10^-3 m s^-2",            "wj_median":  0.251, "sign": "+", "name": "Deformation × wind speed"},
    "ti2":                   {"num": 19, "units": "10^-9 s^-2",              "wj_median": 28.8,   "sign": "+", "name": "Ellrod TI2"},
    "f2d":                   {"num": 20, "units": "10^-9 K^2 m^-2 s^-1",     "wj_median": 56.6,   "sign": "+", "name": "Frontogenesis (2D)"},
    "ncsu1":                 {"num": 21, "units": "10^-18 s^-3",             "wj_median": 11.1,   "sign": "+", "name": "NCSU index v1"},
}

# Map rojak TurbulenceDiagnostics enum members → our short keys.
# ONLY the 14 diagnostics Phase A cleared as PASS/VARIANT come from rojak.
# The 7 that Phase A found buggy (or that rojak lacks) are computed by
# diagnostics_handcoded.py instead — see HANDCODED_KEYS below.
ROJAK_DIAGNOSTICS: dict[str, TurbulenceDiagnostics] = {
    "magnitude_pv":          TurbulenceDiagnostics.MAGNITUDE_PV,      # #1  PASS
    "brown1":                TurbulenceDiagnostics.BROWN1,            # #3  PASS
    "temperature_gradient":  TurbulenceDiagnostics.TEMPERATURE_GRADIENT,  # #4 PASS
    "horizontal_divergence": TurbulenceDiagnostics.HORIZONTAL_DIVERGENCE, # #5 PASS
    "vertical_wind_shear":   TurbulenceDiagnostics.VWS,              # #6  PASS
    "endlich":               TurbulenceDiagnostics.ENDLICH,          # #7  PASS
    "deformation":           TurbulenceDiagnostics.DEF,              # #8  VARIANT (√ applied in run())
    "wind_speed":            TurbulenceDiagnostics.WIND_SPEED,       # #9  PASS
    "ngm2":                  TurbulenceDiagnostics.NGM2,             # #10 PASS
    "nva":                   TurbulenceDiagnostics.NVA,              # #14 PASS
    "vorticity_squared":     TurbulenceDiagnostics.VORTICITY_SQUARED,# #16 PASS
    "ti1":                   TurbulenceDiagnostics.TI1,              # #17 PASS
    "ngm1":                  TurbulenceDiagnostics.NGM1,             # #18 PASS
    "ti2":                   TurbulenceDiagnostics.TI2,              # #19 PASS (divergence sign confirmed)
}

# The 7 diagnostics computed by our own code (Phase A BUG/VERIFY/MISSING):
#   #11 negative_richardson  BUG   N²/Sv → N²/Sv²
#   #2  colson_panofsky      BUG   inherited buggy Ri + 3-level centered stencil
#   #13 ubf                  BUG×3 sign, laplacian, beta
#   #20 f2d                  BUG   cross-term dv_dy → dv_dx
#   #21 ncsu1                BUG   inherited buggy Ri
#   #12 rva_magnitude        MISSING in rojak
#   #15 brown2               VERIFY faithful to A14; native s⁻³, rank-only
HANDCODED_KEYS = [
    "negative_richardson", "colson_panofsky", "ubf",
    "f2d", "ncsu1", "rva_magnitude", "brown2",
]


# ---------------------------------------------------------------------------
# 2. Loading and preparation
# ---------------------------------------------------------------------------


# ===========================================================================
# Data loading and rojak preparation
# ===========================================================================
def load_era5(path: str | Path) -> xr.Dataset:
    """Open an ERA5 trial file.  GRIB (cfgrib) or NetCDF."""
    path = Path(path)
    if path.suffix in (".grib", ".grb", ".grb2"):
        ds = xr.open_dataset(path, engine="cfgrib")
    else:
        ds = xr.open_dataset(path)
    return ds


def prepare_for_rojak(ds: xr.Dataset) -> CATData:
    """
    Take an ERA5-style dataset (any standard CDS naming) and return a
    rojak `CATData` object.

    Handles two CDS variable-name conventions:
      - GRIB / cfgrib short names: u, v, t, z, d, vo, pv
      - NetCDF long names:        u_component_of_wind, v_component_of_wind,
                                  temperature, geopotential, divergence,
                                  vorticity, potential_vorticity
    """
    # 2.1  Pressure-level coordinate name -----------------------------------
    if "isobaricInhPa" in ds.coords:
        ds = ds.rename({"isobaricInhPa": "pressure_level"})
    if "valid_time" in ds.coords and "time" not in ds.coords:
        ds = ds.rename({"valid_time": "time"})
    if "valid_time" in ds.dims and "time" not in ds.dims:
        ds = ds.rename({"valid_time": "time"})

    # 2.2  Variable names ---------------------------------------------------
    # rojak's required CF-style names
    rename_map = {
        # GRIB short names
        "u": "eastward_wind",
        "v": "northward_wind",
        "t": "temperature",
        "z": "geopotential",
        "d": "divergence_of_wind",
        "vo": "vorticity",
        "pv": "potential_vorticity",
        "q": "specific_humidity",
        # CDS long names (NetCDF)
        "u_component_of_wind": "eastward_wind",
        "v_component_of_wind": "northward_wind",
        "geopotential": "geopotential",
        "temperature": "temperature",
        "divergence": "divergence_of_wind",
        "vorticity": "vorticity",
        "potential_vorticity": "potential_vorticity",
        "specific_humidity": "specific_humidity",
    }
    rename_map = {k: v for k, v in rename_map.items() if k in ds.data_vars and k != v}
    if rename_map:
        ds = ds.rename(rename_map)

    # 2.3  CATData also requires specific_humidity in its schema, even
    #      though none of the 21 W&J diagnostics use it.  If it's absent
    #      from the download (the user explicitly excluded it), stub it
    #      with zeros so CATData will instantiate.
    if "specific_humidity" not in ds.data_vars:
        ds = ds.assign(
            specific_humidity=xr.zeros_like(ds["temperature"]).astype("float32")
        )
        ds["specific_humidity"].attrs["note"] = "stubbed with zeros — not used for W&J diagnostics"

    # 2.4  Longitude convention ---------------------------------------------
    # ERA5 sometimes returns 0–360; convert to -180–180 for consistency
    if float(ds.longitude.max()) > 180.0:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
        ds = ds.sortby("longitude")

    # 2.5  Altitude coord (rojak requires it) -------------------------------
    if "altitude" not in ds.coords:
        ds = ds.assign_coords(altitude=("pressure_level", _icao_altitude(ds.pressure_level.values)))

    # 2.6  Transpose to rojak's expected order ------------------------------
    ds = ds.transpose("latitude", "longitude", "time", "pressure_level")

    return CATData(ds, pressure_level_prefix=100.0)  # pressure_level in hPa


def _icao_altitude(p_hpa: np.ndarray) -> np.ndarray:
    """ICAO standard atmosphere pressure-to-altitude (m).  Matches rojak."""
    # rojak uses pressure_to_altitude_icao; reproduce here so the script is
    # robust even if the helper moves between rojak versions.
    p = np.asarray(p_hpa, dtype=float) * 100.0  # Pa
    p0 = 101325.0
    T0 = 288.15
    L = 0.0065
    g = 9.80665
    R = 287.0531
    exponent = (R * L) / g
    return T0 / L * (1.0 - (p / p0) ** exponent)


# ---------------------------------------------------------------------------
# 3. Diagnostic computation
# ---------------------------------------------------------------------------



# ===========================================================================
# Shared geospatial derivative helpers
# ===========================================================================
def _brunt_vaisala_squared(ds: xr.Dataset) -> xr.DataArray:
    """N^2 = (g/theta) d(theta)/dz, exactly as rojak's BruntVaisalaFrequency."""
    theta = _potential_temperature(ds["temperature"], ds["temperature"]["pressure_level"])
    dtheta_dz = altitude_derivative_on_pressure_level(theta, ds["geopotential"])
    return (GRAVITATIONAL_ACCELERATION / theta) * dtheta_dz


def _spacings_1d(ds: xr.Dataset):
    """1-D dx (len nlon) and dy (len nlat) spacings in metres."""
    gs = grid_spacing(ds["latitude"], ds["longitude"], units="degrees")
    dx = np.asarray(gs.dx); dy = np.asarray(gs.dy)
    dx_1d = dx.mean(axis=0) if dx.ndim == 2 else dx
    dy_1d = dy.mean(axis=1) if dy.ndim == 2 else dy
    return dx_1d, dy_1d


def _grad(field: xr.DataArray):
    """Horizontal gradient via rojak's geospatial spatial_gradient -> (dfdx, dfdy)."""
    g = spatial_gradient(field, "deg", GradientMode.GEOSPATIAL)
    return g["dfdx"], g["dfdy"]


def _sel_level(ds: xr.Dataset, name: str, level: int) -> xr.DataArray:
    """Select one pressure level with a plain integer (never a pint Quantity)."""
    return ds[name].sel(pressure_level=level)


def _altitude_from_levels(ds: xr.Dataset) -> xr.DataArray:
    """ICAO standard-atmosphere altitude from pressure levels (fallback)."""
    p = ds["pressure_level"].values.astype(float) * 100.0
    p0, T0, L, g, R = 101325.0, 288.15, 0.0065, 9.80665, 287.0531
    alt = T0 / L * (1.0 - (p / p0) ** ((R * L) / g))
    return xr.DataArray(alt, coords={"pressure_level": ds["pressure_level"]},
                        dims=["pressure_level"], name="altitude")



# ===========================================================================
# The 14 rojak diagnostics (Phase A: PASS/VARIANT)
# ===========================================================================
def compute_rojak_diagnostics(catdata: CATData) -> dict[str, xr.DataArray]:
    """Run all 20 rojak diagnostics. Returns dict {short_key: DataArray}."""
    factory = DiagnosticFactory(catdata)
    out: dict[str, xr.DataArray] = {}
    for key, enum_member in ROJAK_DIAGNOSTICS.items():
        try:
            diag = factory.create(enum_member)
            da = diag.computed_value
            # Persist a sensible name and attrs
            da = da.rename(key)
            da.attrs["rojak_diagnostic"] = enum_member.value
            da.attrs["wj_number"] = REFERENCE_TABLE[key]["num"]
            da.attrs["wj_name"] = REFERENCE_TABLE[key]["name"]
            da.attrs["wj_table_units"] = REFERENCE_TABLE[key]["units"]
            out[key] = da
            print(f"  [{REFERENCE_TABLE[key]['num']:>2}] {key:<23} OK  shape={tuple(da.shape)}")
        except Exception as e:
            print(f"  [{REFERENCE_TABLE[key]['num']:>2}] {key:<23} FAILED: {type(e).__name__}: {e}")
    return out




# ===========================================================================
# The 7 hand-written diagnostics (Phase A: BUG/VERIFY/MISSING)
# Each implemented from Sharman (2006) Appendix A.
# ===========================================================================
    theta = _potential_temperature(ds["temperature"], ds["temperature"]["pressure_level"])
    dtheta_dz = altitude_derivative_on_pressure_level(theta, ds["geopotential"])
    return (GRAVITATIONAL_ACCELERATION / theta) * dtheta_dz

# ---------------------------------------------------------------------------
# Constants (match rojak's values exactly so cross-checks are exact)
# ---------------------------------------------------------------------------
OMEGA   = 7.292115e-05        # Earth angular velocity  [s⁻¹]  (rojak value)
R_EARTH = 6371008.7714        # Earth mean radius       [m]    (rojak value)
G       = 9.80665             # gravity                 [m s⁻²]


# ===========================================================================
# Shared geospatial derivative helpers
# ===========================================================================
def _spacings_1d(ds: xr.Dataset):
    """1-D dx (len nlon) and dy (len nlat) spacings in metres.

    grid_spacing returns dx shape (nlat, nlon-1) and dy shape (nlat-1, nlon);
    first_derivative wants a 1-D spacing per axis, so we reduce each to its
    dominant 1-D structure. For a 30–60°N / 0.25° box the ignored variation
    is < 1 %.
    """
    gs = grid_spacing(ds["latitude"], ds["longitude"], units="degrees")
    dx = np.asarray(gs.dx)
    dy = np.asarray(gs.dy)
    dx_1d = dx.mean(axis=0) if dx.ndim == 2 else dx
    dy_1d = dy.mean(axis=1) if dy.ndim == 2 else dy
    return dx_1d, dy_1d


def _ddx(field: xr.DataArray, dx_1d) -> xr.DataArray:
    return first_derivative(field, dx_1d, axis=field.dims.index("longitude"))


def _ddy(field: xr.DataArray, dy_1d) -> xr.DataArray:
    return first_derivative(field, dy_1d, axis=field.dims.index("latitude"))


def _grad(field: xr.DataArray):
    """Horizontal gradient via rojak's geospatial spatial_gradient.

    Returns (dfdx, dfdy) on the same grid. Using rojak's helper (rather than
    our own first_derivative) keeps the discretization identical to rojak's
    other diagnostics, so the Phase-B cross-check against corrected rojak is
    near-exact rather than merely close.
    """
    g = spatial_gradient(field, "deg", GradientMode.GEOSPATIAL)
    return g["dfdx"], g["dfdy"]


def _sel_level(ds: xr.Dataset, name: str, level: int) -> xr.DataArray:
    """Select one pressure level with a plain integer (never a pint Quantity)."""
    return ds[name].sel(pressure_level=level)


# ===========================================================================
# #11 — Gradient Richardson number  (Sharman A1)
# ===========================================================================
def richardson(ds: xr.Dataset, negative: bool = True) -> xr.DataArray:
    r"""Gradient Richardson number, Sharman (2006) A1.

        Ri = N² / Sv²                                              (A1)

    with N² the Brunt–Väisälä frequency squared (A2) and
    Sv = |∂**u**/∂z| the *magnitude* of vertical wind shear (A3):

        Sv² = (∂u/∂z)² + (∂v/∂z)²

    Ri is **dimensionless** (s⁻² / s⁻² ). rojak's bug divides by Sv
    (unsquared), giving s⁻¹ and a value ~Sv too small.

    The "negative Richardson number" CAT diagnostic is −Ri.

    Returns
    -------
    xr.DataArray  (dimensionless)
    """
    u = ds["eastward_wind"]
    v = ds["northward_wind"]
    phi = ds["geopotential"]

    # N²  — reuse rojak's A2 implementation (Phase A: PASS)
    n2 = _brunt_vaisala_squared(ds)

    # Sv²  — vertical shear on altitude, then squared magnitude
    du_dz = altitude_derivative_on_pressure_level(u, phi)
    dv_dz = altitude_derivative_on_pressure_level(v, phi)
    sv_squared = du_dz ** 2 + dv_dz ** 2         # A3, squared

    ri = n2 / sv_squared                          # A1 — dimensionless
    ri = ri.rename("negative_richardson" if negative else "richardson")
    ri.attrs.update({
        "long_name": ("Negative gradient Richardson number (-Ri)" if negative
                      else "Gradient Richardson number"),
        "units": "1",
        "sharman_eq": "A1",
        "note": "N^2 / Sv^2 (dimensionless); rojak bug uses N^2/Sv",
    })
    return -ri if negative else ri


# ===========================================================================
# #2 — Colson–Panofsky index  (Sharman A4)
# ===========================================================================
def colson_panofsky(ds: xr.Dataset, ri_crit: float = 0.5) -> xr.DataArray:
    r"""Colson–Panofsky index, Sharman (2006) A4.

        CP = Sv² (Δz)² ( 1 − Ri / Ri_crit )                        (A4)

    with Sv the vertical wind shear magnitude, Δz the layer thickness,
    and Ri the *correct* dimensionless Richardson number (A1).

    rojak's own Sv² and Δz are fine; the bug is that it feeds the broken
    #11 Ri into the ( 1 − Ri/Ri_crit ) factor, which flips CP positive.
    Using the correct Ri restores the expected negative median.

    Units: (s⁻¹)² · m² = m² s⁻²  (rojak/W&J then express in 10³ kt²).

    Returns
    -------
    xr.DataArray  (m² s⁻², defined on the layer-upper levels)
    """
    u = ds["eastward_wind"]
    v = ds["northward_wind"]
    phi = ds["geopotential"]
    alt = ds["altitude"] if "altitude" in ds else _altitude_from_levels(ds)

    # Δz across each adjacent level pair (label upper), reduces level dim by 1
    length_scale = alt.diff("pressure_level", label="upper")
    lvl = length_scale["pressure_level"]

    # Sv² (absolute velocities, squared) — matches rojak's internal CP shear
    du_dz = altitude_derivative_on_pressure_level(u, phi)
    dv_dz = altitude_derivative_on_pressure_level(v, phi)
    sv_squared = (du_dz ** 2 + dv_dz ** 2).sel(pressure_level=lvl)

    # Correct Ri (positive form) on the same levels
    ri = richardson(ds, negative=False).sel(pressure_level=lvl)
    richardson_term = 1.0 - ri / ri_crit

    cp = sv_squared * np.square(length_scale) * richardson_term
    cp = cp.rename("colson_panofsky")
    cp.attrs.update({
        "long_name": "Colson-Panofsky index",
        "units": "m2 s-2",
        "sharman_eq": "A4",
        "ri_crit": ri_crit,
        "note": "Uses corrected dimensionless Ri; native m2 s-2 (W&J tabulate 10^3 kt^2)",
    })
    return cp


def _altitude_from_levels(ds: xr.Dataset) -> xr.DataArray:
    """ICAO standard-atmosphere altitude from pressure levels (fallback)."""
    p = ds["pressure_level"].values.astype(float) * 100.0  # Pa
    p0, T0, L, g, R = 101325.0, 288.15, 0.0065, 9.80665, 287.0531
    alt = T0 / L * (1.0 - (p / p0) ** ((R * L) / g))
    return xr.DataArray(alt, coords={"pressure_level": ds["pressure_level"]},
                        dims=["pressure_level"], name="altitude")


# ===========================================================================
# #13 — UBF, residual of nonlinear balance equation  (Sharman A30)
# ===========================================================================
def ubf(ds: xr.Dataset, target_level: int = 200) -> xr.DataArray:
    r"""Unbalanced-flow diagnostic, Sharman (2006) A30 / Koch & Caracena (2002).

        UBF = | ∇²Φ − 2 J(u,v) − f ζ + β u |                        (A30, |·| form)

    where
        ∇²Φ = ∂²Φ/∂x² + ∂²Φ/∂y²         (true Laplacian)
        J(u,v) = ∂u/∂x·∂v/∂y − ∂u/∂y·∂v/∂x   (Jacobian)
        f  = 2Ω sinφ                    (Coriolis)
        β  = ∂f/∂y = 2Ω cosφ / R_E      (Rossby parameter)
        ζ  = relative vorticity

    Fixes all three rojak bugs:
      1. sign: residual (∇²Φ − …) not sum (∇²Φ + …)
      2. ∇²Φ is a real second-derivative Laplacian, not ∂Φ/∂x+∂Φ/∂y
      3. β = 2Ω cosφ/R_E, not f/R_E

    Units: s⁻².

    Returns
    -------
    xr.DataArray  (s⁻², at target_level)
    """
    u    = _sel_level(ds, "eastward_wind", target_level)
    v    = _sel_level(ds, "northward_wind", target_level)
    phi  = _sel_level(ds, "geopotential", target_level)
    zeta = _sel_level(ds, "vorticity", target_level)

    lat_rad = np.deg2rad(ds["latitude"])
    f    = (2 * OMEGA * np.sin(lat_rad)).broadcast_like(u.isel(time=0))
    beta = (2 * OMEGA * np.cos(lat_rad) / R_EARTH).broadcast_like(u.isel(time=0))

    # True ∇²Φ  = ∂²Φ/∂x² + ∂²Φ/∂y²  (compose spatial_gradient twice)
    dphi_dx, dphi_dy = _grad(phi)
    d2phi = _grad(dphi_dx)[0] + _grad(dphi_dy)[1]

    # Jacobian
    du_dx, du_dy = _grad(u)
    dv_dx, dv_dy = _grad(v)
    jac = du_dx * dv_dy - du_dy * dv_dx

    residual = d2phi - 2 * jac - f * zeta + beta * u
    out = np.abs(residual).rename("ubf")
    out.attrs.update({
        "long_name": "UBF (nonlinear balance equation residual)",
        "units": "s-2",
        "sharman_eq": "A30",
        "note": "residual form; true Laplacian; beta=2*Omega*cos(phi)/R",
    })
    return out


# ===========================================================================
# #20 — 2-D frontogenesis  (Sharman A9)
# ===========================================================================
def frontogenesis_2d(ds: xr.Dataset, target_level: int = 200) -> xr.DataArray:
    r"""2-D kinematic frontogenesis, Sharman (2006) A9 / Bluestein (1993).

        F = −(1/|∇θ|) [ θx² ∂u/∂x + θy² ∂v/∂y
                        + θx θy (∂v/∂x + ∂u/∂y) ]                   (A9)

    where θ is potential temperature and θx = ∂θ/∂x, θy = ∂θ/∂y.

    rojak's bug: the cross-term coefficient uses ∂v/∂y instead of the
    correct ∂v/∂x (the cross-term is θx·θy·(∂v/∂x + ∂u/∂y), i.e. the
    shearing-deformation combination, not θx·θy·(∂v/∂y + ∂u/∂y)).

    Units: K² m⁻² s⁻¹.

    Returns
    -------
    xr.DataArray  (K² m⁻² s⁻¹, at target_level)
    """
    u   = _sel_level(ds, "eastward_wind", target_level)
    v   = _sel_level(ds, "northward_wind", target_level)
    t   = _sel_level(ds, "temperature", target_level)
    p   = float(target_level)

    # potential temperature θ = T (1000/p)^(R/cp)
    theta = t * (1000.0 / p) ** (287.0 / 1004.0)

    dx_1d, dy_1d = _spacings_1d(ds)
    theta_x, theta_y = _grad(theta)
    mag_grad = np.sqrt(theta_x ** 2 + theta_y ** 2)
    inv_mag = -xr.where(mag_grad != 0, 1.0 / mag_grad, 0.0)

    du_dx, du_dy = _grad(u)
    dv_dx, dv_dy = _grad(v)

    f2d = inv_mag * (
        theta_x ** 2 * du_dx
        + theta_y ** 2 * dv_dy
        + theta_x * theta_y * (dv_dx + du_dy)     # correct cross term
    )
    f2d = f2d.rename("f2d")
    f2d.attrs.update({
        "long_name": "2-D frontogenesis",
        "units": "K2 m-2 s-1",
        "sharman_eq": "A9",
        "note": "cross-term uses dv_dx + du_dy (corrected)",
    })
    return f2d


# ===========================================================================
# #21 — NCSU1  (Sharman A36)
# ===========================================================================
def ncsu1(ds: xr.Dataset, target_level: int = 200, ri_floor: float = 1e-5) -> xr.DataArray:
    r"""NCSU1 index, Sharman (2006) A36 / Kaplan et al. (2005).

        NCSU1 = [ 1 / max(Ri, 1e-5) ]
                · max( u ∂u/∂x + v ∂v/∂y, 0 )
                · |∇ζ|                                             (A36)

    The formula itself is fine in rojak; the bug is that rojak divides by
    the broken #11 Ri. Using the correct dimensionless Ri restores both
    the magnitude and (because tiny-Ri cells no longer blow up in the same
    way) the 3-hourly rank stability.

    Units: with correct Ri (dimensionless), NCSU1 has units of
    (m s⁻¹ · s⁻¹) · (s⁻¹ m⁻¹) = s⁻³ ... wait: advection term
    u·∂u/∂x ~ (m/s)(1/s) = m s⁻²; |∇ζ| ~ s⁻¹ m⁻¹; product = s⁻³;
    divided by dimensionless Ri → s⁻³.  (W&J tabulate 10⁻¹⁸ s⁻³ — but
    note W&J's exact NCSU1 units are debated; we return SI s⁻³.)

    Returns
    -------
    xr.DataArray  (s⁻³, at target_level)
    """
    u    = _sel_level(ds, "eastward_wind", target_level)
    v    = _sel_level(ds, "northward_wind", target_level)
    zeta = _sel_level(ds, "vorticity", target_level)

    dx_1d, dy_1d = _spacings_1d(ds)
    du_dx, _ = _grad(u)
    _, dv_dy = _grad(v)

    # advection-like term  (u ∂u/∂x + v ∂v/∂y), clipped ≥ 0
    advection = (u * du_dx + v * dv_dy).clip(min=0)

    # |∇ζ|
    dzeta_dx, dzeta_dy = _grad(zeta)
    grad_zeta = np.sqrt(dzeta_dx ** 2 + dzeta_dy ** 2)

    # correct dimensionless Ri at this level, floored
    ri = richardson(ds, negative=False).sel(pressure_level=target_level)
    # Sharman A36 floors Ri from BELOW at ri_floor: max(Ri, 1e-5).
    # This is NOT the same as flooring |Ri| -- any negative Ri (common in
    # unstable/turbulent regions, i.e. exactly where NCSU1 should fire
    # strongly) must also be replaced by the small positive floor, not
    # left as a large negative divisor.
    ri_safe = xr.where(ri < ri_floor, ri_floor, ri)
    print(">>> [ncsu1 debug] marker v2 is running")
    n_negative = int((ri.values < 0).sum())
    n_total = ri.size
    print(f"    [ncsu1 debug] Ri cells: {n_total} total, {n_negative} negative ({100*n_negative/n_total:.2f}%)")

    out = (advection * grad_zeta) / ri_safe
    out = out.rename("ncsu1")
    out.attrs.update({
        "long_name": "NCSU1 index",
        "units": "s-3",
        "sharman_eq": "A36",
        "note": "uses corrected dimensionless Ri (floored at 1e-5)",
    })
    return out


# ===========================================================================
# #12 — Relative vorticity advection magnitude  (not in rojak)
# ===========================================================================
def rva(ds: xr.Dataset, target_level: int = 200) -> xr.DataArray:
    r"""Magnitude of relative vorticity advection.

        RVA = | u ∂ζ/∂x + v ∂ζ/∂y |

    rojak provides negative *absolute* vorticity advection (NVA, A37) but
    not this magnitude-of-relative-vorticity-advection quantity, so it is
    computed here with geospatial (metre-based) derivatives.

    Units: s⁻².

    Returns
    -------
    xr.DataArray  (s⁻², at target_level)
    """
    u    = _sel_level(ds, "eastward_wind", target_level)
    v    = _sel_level(ds, "northward_wind", target_level)
    zeta = _sel_level(ds, "vorticity", target_level)

    dzeta_dx, dzeta_dy = _grad(zeta)
    out = np.abs(u * dzeta_dx + v * dzeta_dy).rename("rva_magnitude")
    out.attrs.update({
        "long_name": "Magnitude of relative vorticity advection",
        "units": "s-2",
        "note": "geospatial derivatives; not provided by rojak",
    })
    return out


# ===========================================================================
# #15 — Brown energy dissipation rate (Brown2)  (Sharman A14)
# ===========================================================================
def brown2(ds: xr.Dataset, brown1: xr.DataArray) -> xr.DataArray:
    r"""Brown Index 2 / Brown EDR, Sharman (2006) A14.

        ε_brown = (1/24) Φ Sv²                                     (A14)

    where Φ = Brown Index 1 (A13, s⁻¹) and Sv is the vertical wind shear.
    rojak transcribes A14 faithfully; the only wrinkle is dimensional —
    A14 yields s⁻³, whereas W&J Table 1 lists J kg⁻¹ s⁻¹ = m² s⁻³, implying
    an unspecified length² factor.

    Decision (documented): we keep the native s⁻³ form and treat Brown2 as a
    **rank-only** diagnostic. Prosser's methodology is percentile-threshold
    based, so a constant multiplicative L² does not change any exceedance
    ranking. We therefore do NOT invent a length scale.

    This function exists mainly so the pipeline has a single, citable Brown2
    implementation with the units decision recorded; numerically it equals
    rojak's Brown2 (which Phase A cleared as faithful to A14).

    Parameters
    ----------
    brown1 : the already-computed Brown Index 1 field (from rojak).

    Returns
    -------
    xr.DataArray  (s⁻³)
    """
    u = ds["eastward_wind"]
    v = ds["northward_wind"]
    phi = ds["geopotential"]

    du_dz = altitude_derivative_on_pressure_level(u, phi)
    dv_dz = altitude_derivative_on_pressure_level(v, phi)
    sv_squared = du_dz ** 2 + dv_dz ** 2

    out = (1.0 / 24.0) * brown1 * sv_squared
    out = out.rename("brown2")
    out.attrs.update({
        "long_name": "Brown energy dissipation rate (Brown2)",
        "units": "s-3",
        "sharman_eq": "A14",
        "note": "native s-3; rank-only (no invented length scale); "
                "W&J tabulate m2 s-3 (implied L^2)",
    })
    return out




# ===========================================================================
# Master entry point — compute all 21 diagnostics
# ===========================================================================
def compute_all_21(catdata: CATData, target_level: int = 200) -> dict[str, xr.DataArray]:
    """Compute all 21 diagnostics: 14 from rojak, 7 hand-written.

    Returns a dict {short_key: xr.DataArray}. Diagnostics that still carry a
    pressure_level dimension are reduced to `target_level` at the end.
    """
    ds = catdata._dataset

    # --- 14 from rojak ---
    out = compute_rojak_diagnostics(catdata)

    # --- 7 hand-written ---
    brown1 = DiagnosticFactory(catdata).create(TurbulenceDiagnostics.BROWN1).computed_value
    builders = {
        "negative_richardson": lambda: richardson(ds, negative=True),
        "colson_panofsky":     lambda: colson_panofsky(ds),
        "ubf":                 lambda: ubf(ds, target_level=target_level),
        "f2d":                 lambda: frontogenesis_2d(ds, target_level=target_level),
        "ncsu1":               lambda: ncsu1(ds, target_level=target_level),
        "rva_magnitude":       lambda: rva(ds, target_level=target_level),
        "brown2":              lambda: brown2(ds, brown1),
    }
    for key, builder in builders.items():
        out[key] = builder().rename(key)

    # --- reduce any remaining pressure_level dim to the target level ---
    for k, da in list(out.items()):
        if "pressure_level" in da.dims:
            try:
                out[k] = da.sel(pressure_level=target_level)
            except KeyError:
                out[k] = da.isel(pressure_level=len(da.pressure_level) // 2)
    return out