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
  #21 NCSU1            inherits buggy Ri; du/dx,dv/dy need projection-corrected
                        vector_derivatives, not plain spatial_gradient
  #12 RVA              not provided by rojak
  #15 Brown2           faithful to A14 but native s^-3 (rank-only); no invented L^2

References: Sharman et al. (2006) Wea. Forecasting App. A; Williams & Joshi
(2013); Ellrod & Knapp (1992); Kaplan et al. (2005); Koch & Caracena (2002).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import xarray as xr

from rojak.core.data import CATData
from rojak.core.derivatives import (
    grid_spacing, first_derivative, spatial_gradient, GradientMode,
    vector_derivatives, VelocityDerivative,
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
    "colson_panofsky":       {"num": 2,  "units": "10^3 kt^2",               "wj_median": -34.8,  "sign": "+", "name": "Colson–Panofsky index"},  # Q-AGG-3: was "either" -- Williams(2017) Table 2 is a single monotonic ladder (light=-29.3 -> severe=-22.2), i.e. one-tailed; negative median != two-tailed criterion
    "brown1":                {"num": 3,  "units": "10^-6 s^-1",              "wj_median": 77.1,   "sign": "+", "name": "Brown index"},
    "temperature_gradient":  {"num": 4,  "units": "10^-6 K m^-1",            "wj_median":  5.75,  "sign": "+", "name": "|Horizontal temperature gradient|"},
    "horizontal_divergence": {"num": 5,  "units": "10^-6 s^-1",              "wj_median":  2.82,  "sign": "+", "name": "|Horizontal divergence|"},
    "vertical_wind_shear":   {"num": 6,  "units": "10^-3 s^-1",              "wj_median":  1.88,  "sign": "+", "name": "|Vertical wind shear|"},
    "endlich":               {"num": 7,  "units": "10^-3 rad s^-1",          "wj_median":  0.952, "sign": "+", "name": "Wind speed × directional shear"},
    "deformation":           {"num": 8,  "units": "10^-6 s^-1",              "wj_median": 18.6,   "sign": "+", "name": "Flow deformation (rojak returns squared)"},
    "wind_speed":            {"num": 9,  "units": "m s^-1",                  "wj_median": 14.9,   "sign": "+", "name": "Wind speed"},
    "ngm2":                  {"num": 10, "units": "10^-9 K m^-1 s^-1",       "wj_median":  8.17,  "sign": "+", "name": "Deformation × vertical T gradient"},
    "negative_richardson":   {"num": 11, "units": "dimensionless",           "wj_median": -127.2, "sign": "+", "name": "Negative Richardson number (-Ri)"},  # Q-AGG-3: was "either" -- same reasoning as colson_panofsky; Williams(2017) Table 2 ladder is monotonic (light=-15.4 -> severe=-5.9), one-tailed
    "rva_magnitude":         {"num": 12, "units": "10^-10 s^-2",             "wj_median":  2.33,  "sign": "+", "name": "|Relative vorticity advection|"},
    "ubf":                   {"num": 13, "units": "10^-12 s^-2",             "wj_median": 161.0,  "sign": "+", "name": "|Residual of nonlinear balance eq.|"},
    "nva":                   {"num": 14, "units": "10^-10 s^-2",             "wj_median":  2.05,  "sign": "+", "name": "Negative absolute vorticity advection"},
    "brown2":                {"num": 15, "units": "10^-6 J kg^-1 s^-1",      "wj_median": 116.0,  "sign": "+", "name": "Brown energy dissipation rate"},
    "vorticity_squared":     {"num": 16, "units": "10^-9 s^-2",              "wj_median":  0.221, "sign": "+", "name": "Relative vorticity squared"},
    "ti1":                   {"num": 17, "units": "10^-9 s^-2",              "wj_median": 31.5,   "sign": "+", "name": "Ellrod TI1"},
    "ngm1":                  {"num": 18, "units": "10^-3 m s^-2",            "wj_median":  0.251, "sign": "+", "name": "Deformation × wind speed"},
    "ti2":                   {"num": 19, "units": "10^-9 s^-2",              "wj_median": 28.8,   "sign": "+", "name": "Ellrod TI2"},
    "f2d":                   {"num": 20, "units": "10^-9 m^2 s^-3 K^-2",     "wj_median": 56.6,   "sign": "+", "name": "Frontogenesis (2D)"},  # Q-UNITS-1: was stale "10^-9 K^2 m^-2 s^-1" (Miller-form units); corrected to A9/isentropic units, agrees with Williams(2017) and W&J(2013)
    "ncsu1":                 {"num": 21, "units": "10^-18 s^-3",             "wj_median": 11.1,   "sign": "+", "name": "NCSU index v1"},
}

# Map rojak TurbulenceDiagnostics enum members → our short keys.
# ONLY the 14 diagnostics Phase A cleared as PASS/VARIANT come from rojak.
# The 7 that Phase A found buggy (or that rojak lacks) are hand-written
# below — see HANDCODED_KEYS.
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
#   #21 ncsu1                BUG   inherited buggy Ri + needs projection-corrected derivatives
#   #12 rva_magnitude        MISSING in rojak
#   #15 brown2               VERIFY faithful to A14; native s⁻³, rank-only
HANDCODED_KEYS = [
    "negative_richardson", "colson_panofsky", "ubf",
    "f2d", "ncsu1", "rva_magnitude", "brown2",
]


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
    # Pressure-level coordinate name
    if "isobaricInhPa" in ds.coords:
        ds = ds.rename({"isobaricInhPa": "pressure_level"})
    if "valid_time" in ds.coords and "time" not in ds.coords:
        ds = ds.rename({"valid_time": "time"})
    if "valid_time" in ds.dims and "time" not in ds.dims:
        ds = ds.rename({"valid_time": "time"})

    # Variable names -> rojak's required CF-style names
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

    # CATData requires specific_humidity in its schema even though none of
    # the 21 W&J diagnostics use it. Stub it with zeros if absent.
    if "specific_humidity" not in ds.data_vars:
        ds = ds.assign(
            specific_humidity=xr.zeros_like(ds["temperature"]).astype("float32")
        )
        ds["specific_humidity"].attrs["note"] = "stubbed with zeros — not used for W&J diagnostics"

    # Longitude convention: ERA5 sometimes returns 0–360; convert to -180–180
    if float(ds.longitude.max()) > 180.0:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
        ds = ds.sortby("longitude")

    # Altitude coord (rojak requires it)
    if "altitude" not in ds.coords:
        ds = ds.assign_coords(altitude=("pressure_level", _icao_altitude(ds.pressure_level.values)))

    # Transpose to rojak's expected order
    ds = ds.transpose("latitude", "longitude", "time", "pressure_level")

    return CATData(ds, pressure_level_prefix=100.0)  # pressure_level in hPa


def _icao_altitude(p_hpa: np.ndarray) -> np.ndarray:
    """ICAO standard atmosphere pressure-to-altitude (m).  Matches rojak."""
    p = np.asarray(p_hpa, dtype=float) * 100.0  # Pa
    p0 = 101325.0
    T0 = 288.15
    L = 0.0065
    g = 9.80665
    R = 287.0531
    exponent = (R * L) / g
    return T0 / L * (1.0 - (p / p0) ** exponent)


# ===========================================================================
# Shared geospatial derivative helpers
# ===========================================================================
def _brunt_vaisala_squared(ds: xr.Dataset) -> xr.DataArray:
    """N^2 = (g/theta) d(theta)/dz, exactly as rojak's BruntVaisalaFrequency."""
    theta = _potential_temperature(ds["temperature"], ds["temperature"]["pressure_level"])
    dtheta_dz = altitude_derivative_on_pressure_level(theta, ds["geopotential"])
    return (GRAVITATIONAL_ACCELERATION / theta) * dtheta_dz


def theta_derivative_on_pressure_level(function: xr.DataArray, theta: xr.DataArray,
                                        level_coord_name: str = "pressure_level") -> xr.DataArray:
    r"""Derivative w.r.t. potential temperature for data on pressure levels
    (Q-F2D-5). Same bulk chain-rule pattern as rojak's
    altitude_derivative_on_pressure_level -- just theta instead of Phi/g:

        df/dtheta = (df/dp) / (dtheta/dp)

    No isentropic remapping -- computed directly on the existing
    175/200/225 hPa levels via xarray's .differentiate (centered in the
    interior, one-sided at the top/bottom level), same mechanism
    altitude_derivative_on_pressure_level already uses for df/dp.
    """
    return function.differentiate(level_coord_name) / theta.differentiate(level_coord_name)


def _spacings_1d(ds: xr.Dataset):
    """1-D dx (len nlon) and dy (len nlat) spacings in metres.

    grid_spacing returns dx shape (nlat, nlon-1) and dy shape (nlat-1, nlon);
    first_derivative wants a 1-D spacing per axis, so we reduce each to its
    dominant 1-D structure. For a 30–60°N / 0.25° box the ignored variation
    is < 1 %.
    """
    gs = grid_spacing(ds["latitude"], ds["longitude"], units="degrees")
    dx = np.asarray(gs.dx); dy = np.asarray(gs.dy)
    dx_1d = dx.mean(axis=0) if dx.ndim == 2 else dx
    dy_1d = dy.mean(axis=1) if dy.ndim == 2 else dy
    return dx_1d, dy_1d


def _grad(field: xr.DataArray):
    """Horizontal gradient via rojak's geospatial spatial_gradient -> (dfdx, dfdy).

    Correct for SCALAR fields (vorticity, geopotential, potential temperature).
    NOT correct for differentiating a VECTOR COMPONENT (u or v itself) --
    those need the extra map-projection correction terms in
    rojak.core.derivatives.vector_derivatives(). See ncsu1() below.
    """
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
@dataclass
class DiagnosticFailure:
    """Q-INTEG-4: one rojak-diagnostic-computation failure, kept as a
    first-class, inspectable record -- not just a print statement that
    scrolls off the terminal on a 42-year, 3-hourly run."""
    key: str
    wj_number: int
    exception_type: str
    message: str


def compute_rojak_diagnostics(
    catdata: CATData,
) -> tuple[dict[str, xr.DataArray], list[DiagnosticFailure]]:
    """Run all 14 cleared rojak diagnostics. Returns
    ({short_key: DataArray}, [DiagnosticFailure, ...]).

    Q-INTEG-4: a bare `except Exception` here used to print "FAILED" and
    silently OMIT the diagnostic from the output dict -- over a 42-year,
    3-hourly run (122,728 timesteps), some diagnostic throwing on some odd
    chunk is a "when", not an "if", and dropping it silently means it gets
    averaged over as if it never existed (the exact same failure SHAPE as
    Bug 2's silent level-substitution and Bug 1's NaN-as-non-exceedance:
    a real, plausible-looking result quietly built from less than it
    claims to be built from).

    Fix: on failure, (1) log loudly (unchanged, still prints), (2) record
    a DiagnosticFailure with the diagnostic name and cause, and (3) if
    ANY other diagnostic in this same call succeeded, use its shape/coords
    as a template to emit an all-NaN placeholder for the failed one under
    its own key -- so it flows into Bug 1's skipna + populated_count
    machinery downstream and gets EXCLUDED from the mean at every
    cell/timestep, rather than silently missing from the dict entirely
    (which would silently change which diagnostics feed
    exceedance_mean_single_severity's dict, with no record of it).

    Deliberately does NOT guess a shape from nothing: if every single
    diagnostic in this call fails, there is no valid sibling template,
    and -- more importantly -- no evidence the input chunk itself is
    sound. That is a structural failure, not an isolated numerical
    hiccup, so this raises RuntimeError instead of returning an empty or
    fabricated-shape result.
    """
    factory = DiagnosticFactory(catdata)
    out: dict[str, xr.DataArray] = {}
    failures: list[DiagnosticFailure] = []
    pending: list[tuple[str, TurbulenceDiagnostics, Exception]] = []

    for key, enum_member in ROJAK_DIAGNOSTICS.items():
        try:
            diag = factory.create(enum_member)
            da = diag.computed_value
            da = da.rename(key)
            da.attrs["rojak_diagnostic"] = enum_member.value
            da.attrs["wj_number"] = REFERENCE_TABLE[key]["num"]
            da.attrs["wj_name"] = REFERENCE_TABLE[key]["name"]
            da.attrs["wj_table_units"] = REFERENCE_TABLE[key]["units"]
            out[key] = da
            print(f"  [{REFERENCE_TABLE[key]['num']:>2}] {key:<23} OK  shape={tuple(da.shape)}")
        except Exception as e:
            print(f"  [{REFERENCE_TABLE[key]['num']:>2}] {key:<23} FAILED: {type(e).__name__}: {e}")
            failures.append(DiagnosticFailure(
                key=key, wj_number=REFERENCE_TABLE[key]["num"],
                exception_type=type(e).__name__, message=str(e),
            ))
            pending.append((key, enum_member, e))

    if pending:
        if not out:
            # Every diagnostic failed -- no sibling template exists, and no
            # evidence this chunk's input is sound. Fail loudly rather than
            # return something built from nothing.
            raise RuntimeError(
                f"Q-INTEG-4: ALL {len(ROJAK_DIAGNOSTICS)} rojak diagnostics failed "
                f"for this chunk -- no successful sibling to template a NaN "
                f"placeholder from, and no evidence the input chunk itself is "
                f"sound. Refusing to fabricate a shape from nothing. Failures: "
                f"{[(f.key, f.exception_type, f.message) for f in failures]}"
            )
        template = next(iter(out.values()))
        for key, enum_member, e in pending:
            placeholder = xr.full_like(template, np.nan).rename(key)
            placeholder.attrs.update({
                "rojak_diagnostic": enum_member.value,
                "wj_number": REFERENCE_TABLE[key]["num"],
                "wj_name": REFERENCE_TABLE[key]["name"],
                "wj_table_units": REFERENCE_TABLE[key]["units"],
                "Q_INTEG_4_status": "FAILED_ALL_NAN",
                "Q_INTEG_4_exception_type": type(e).__name__,
                "Q_INTEG_4_message": str(e),
                "Q_INTEG_4_template_source": str(template.name),
            })
            out[key] = placeholder
            print(f"       {key:<23} -> emitting all-NaN placeholder "
                  f"(shape={tuple(placeholder.shape)}, templated from "
                  f"'{template.name}'); will be excluded via skipna "
                  f"downstream (Q-INTEG-3), not counted as non-exceedance")

    return out, failures


# ===========================================================================
# The 7 hand-written diagnostics (Phase A: BUG/VERIFY/MISSING)
# Each implemented from Sharman (2006) Appendix A.
# ===========================================================================

# ---------------------------------------------------------------------------
# #11 — Gradient Richardson number  (Sharman A1)
# ---------------------------------------------------------------------------
def richardson(ds: xr.Dataset, negative: bool = True) -> xr.DataArray:
    r"""Gradient Richardson number, Sharman (2006) A1.

        Ri = N² / Sv²                                              (A1)

    with N² the Brunt–Väisälä frequency squared (A2) and
    Sv = |∂u/∂z| the magnitude of vertical wind shear (A3):

        Sv² = (∂u/∂z)² + (∂v/∂z)²

    Ri is dimensionless. rojak's bug divides by Sv (unsquared), giving
    units of s⁻¹ and a value ~Sv too small.

    The "negative Richardson number" CAT diagnostic is −Ri.
    """
    u = ds["eastward_wind"]
    v = ds["northward_wind"]

    n2 = _brunt_vaisala_squared(ds)

    du_dz = altitude_derivative_on_pressure_level(u, ds["geopotential"])
    dv_dz = altitude_derivative_on_pressure_level(v, ds["geopotential"])
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


# ---------------------------------------------------------------------------
# #2 — Colson–Panofsky index  (Sharman A4)
# ---------------------------------------------------------------------------
def colson_panofsky(ds: xr.Dataset, ri_crit: float = 0.5) -> xr.DataArray:
    r"""Colson–Panofsky index, Sharman (2006) A4.

        CP = Sv² (Δz)² ( 1 − Ri / Ri_crit )                        (A4)
           = (Δz)² ( Sv² − N² / Ri_crit )        [algebraically-reduced form]

    Q-CP-3: implemented directly in the REDUCED form, not via a Ri
    round-trip. Root cause of the earlier rho=0.24-0.35 vs rojak (Q-CP-1):
    rojak's ColsonPanofsky class combines two DIFFERENT Sv^2 values that
    are supposed to be the same number -- Ri's Sv^2 comes from the VWS
    diagnostic on the FULL unsliced 3-level u/v/geopotential, while its
    own multiplicative Sv^2 is recomputed AFTER pre-slicing to only 2
    levels in __init__. These disagree by up to ~3000x at low-shear grid
    cells (verified: cell 54.75N/48W, 2 Jul 2016 00Z -- Sv^2 from the
    Richardson path was 3.9e-9 vs CP's own re-sliced 1.2e-5), so the
    Sv^2*Ri/Ri_crit -> N^2/Ri_crit cancellation never actually happens and
    the residual (N^2 times the ratio of the two mismatched Sv^2's) blows
    up. Confirmed NOT a float32/near-machine-zero rounding artifact
    (Sv^2=1.2e-5 is nowhere near float32 underflow ~1.2e-7) and NOT a data
    issue (N^2 at that cell sits at the 56th percentile of the domain,
    smooth vs. all 5x5 neighbours).

    Fix: compute Sv^2 and N^2 ONCE, from a single consistent stencil, and
    combine them directly with no Ri intermediate. Both terms use
    altitude_derivative_on_pressure_level on the FULL (unsliced) 3-level
    u/v/geopotential -- the same building blocks that already feed the
    (Phase A/B CLEAN, rho=1.0 vs rojak) #6 VWS and #11 Richardson -- then
    select down to the CP output levels only at the very end. This is the
    same stencil hand's own richardson()/colson_panofsky() always used
    internally (du_dz/dv_dz computed before any level selection), so this
    change does not alter hand's own numbers -- it only removes the Ri
    round-trip, which was never actually the mismatched quantity on the
    hand side.

    Units: m² s⁻² (rojak/W&J then express in 10³ kt²).
    """
    u = ds["eastward_wind"]
    v = ds["northward_wind"]
    alt = ds["altitude"] if "altitude" in ds else _altitude_from_levels(ds)

    length_scale = alt.diff("pressure_level", label="upper")
    lvl = length_scale["pressure_level"]

    # Single consistent (full 3-level) stencil for BOTH terms -- computed
    # before any level selection, exactly like richardson()'s n2/sv_squared.
    du_dz = altitude_derivative_on_pressure_level(u, ds["geopotential"])
    dv_dz = altitude_derivative_on_pressure_level(v, ds["geopotential"])
    sv_squared = (du_dz ** 2 + dv_dz ** 2).sel(pressure_level=lvl)
    n2 = _brunt_vaisala_squared(ds).sel(pressure_level=lvl)

    cp = np.square(length_scale) * (sv_squared - n2 / ri_crit)
    cp = cp.rename("colson_panofsky")
    cp.attrs.update({
        "long_name": "Colson-Panofsky index",
        "units": "m2 s-2",
        "sharman_eq": "A4",
        "ri_crit": ri_crit,
        "note": "Q-CP-3: algebraically-reduced form l^2*(Sv^2-N^2/Ri_crit), "
                "single consistent 3-level stencil for both terms, no Ri "
                "round-trip. native m2 s-2 (W&J tabulate 10^3 kt^2)",
    })
    return cp


# ---------------------------------------------------------------------------
# #13 — UBF, residual of nonlinear balance equation  (Sharman A30)
# ---------------------------------------------------------------------------
def ubf(ds: xr.Dataset, target_level: int = 200) -> xr.DataArray:
    r"""Unbalanced-flow diagnostic, Sharman (2006) A30 / Koch & Caracena (2002).

        UBF = | ∇²Φ − 2 J(u,v) − f ζ + β u |                        (A30)

    where ∇²Φ = ∂²Φ/∂x² + ∂²Φ/∂y² (true Laplacian), J(u,v) is the Jacobian,
    f = 2Ω sinφ, β = ∂f/∂y = 2Ω cosφ/R_E, ζ = relative vorticity.

    Fixes three rojak bugs: (1) sign -- residual not sum, (2) ∇²Φ is a
    real second-derivative Laplacian, not ∂Φ/∂x+∂Φ/∂y, (3) β=2Ω cosφ/R_E,
    not f/R_E.

    Fixes two further GEOMETRY errors found on 2026-08-26 (STATUS.md §4g),
    both of which rojak also has, so the cross-check in 4_verify.py could
    never have caught them -- the two implementations shared the mistake:
      (4) ∇²Φ was ∂²Φ/∂x²+∂²Φ/∂y², missing the spherical -tan(φ)/M ∂Φ/∂y
          term. Verified against the divergence theorem: the old form is
          99% wrong. For a geostrophic Φ the missing term is 1-7x the
          published median UBF, and it grows with latitude, so it does NOT
          cancel out of a percentile calibration.
      (5) J(u,v) used scalar gradients of vector components. Verified
          against Stokes' theorem: that operator gets vorticity 100% wrong.
          Smaller effect on UBF (2J is ~0.5% of the residual) but wrong.

    UBF is a RESIDUAL -- a near-cancellation of large terms -- so a
    consistent geometry matters more here than in any other diagnostic. The
    ζ term uses ERA5's own vorticity, which is the true spherical vorticity;
    computing the other terms in a flat-earth geometry meant the residual
    partly measured that inconsistency rather than atmospheric imbalance.

    Units: s⁻².
    """
    u    = _sel_level(ds, "eastward_wind", target_level)
    v    = _sel_level(ds, "northward_wind", target_level)
    phi  = _sel_level(ds, "geopotential", target_level)
    zeta = _sel_level(ds, "vorticity", target_level)

    lat_rad = np.deg2rad(ds["latitude"])
    f    = (2 * OMEGA * np.sin(lat_rad)).broadcast_like(u.isel(time=0))
    beta = (2 * OMEGA * np.cos(lat_rad) / R_EARTH).broadcast_like(u.isel(time=0))

    # ---- geometry (verified 2026-08-26, see STATUS.md 4g) -----------------
    # Both terms below are DIVERGENCES OF VECTOR FIELDS on a curved earth, so
    # both need the map-projection corrections; the plain scalar _grad() is
    # correct only for the gradient of a scalar. This was established by two
    # independent geometric theorems, not by convention:
    #
    #   Stokes (circulation / area) showed that zeta built from scalar
    #   gradients is 100% wrong -- wrong sign -- because a zonal flow carries
    #   curvature vorticity u tan(phi)/M even with no shear.
    #
    #   The divergence theorem (boundary flux / area) showed that
    #   Phi_xx + Phi_yy is 99% wrong as a spherical Laplacian; the missing
    #   -tan(phi)/M * dPhi/dy term is the bulk of the answer, not a refinement.
    #
    # Both corrections are obtained here by asking rojak for the divergence of
    # a vector field rather than hand-coding tan(phi)/M, so the ellipsoidal
    # scale factors stay consistent with every other diagnostic.

    # Laplacian of geopotential = divergence of grad(Phi).
    # grad(Phi) is correct with _grad (Phi is a scalar); its DIVERGENCE is not.
    dphi_dx, dphi_dy = _grad(phi)
    vd_phi = vector_derivatives(
        dphi_dx, dphi_dy, "deg",
        components=[VelocityDerivative.DU_DX, VelocityDerivative.DV_DY])
    d2phi = vd_phi[VelocityDerivative.DU_DX] + vd_phi[VelocityDerivative.DV_DY]

    # Jacobian of the horizontal velocity: derivatives OF vector components.
    # Matches rojak's own CATData.jacobian_horizontal_velocity().
    vd = vector_derivatives(
        u, v, "deg",
        components=[VelocityDerivative.DU_DX, VelocityDerivative.DU_DY,
                    VelocityDerivative.DV_DX, VelocityDerivative.DV_DY])
    du_dx = vd[VelocityDerivative.DU_DX]
    du_dy = vd[VelocityDerivative.DU_DY]
    dv_dx = vd[VelocityDerivative.DV_DX]
    dv_dy = vd[VelocityDerivative.DV_DY]
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


# ---------------------------------------------------------------------------
# #20 — 2-D frontogenesis  (Sharman A9)
# ---------------------------------------------------------------------------
def frontogenesis_2d(ds: xr.Dataset, target_level: int = 200) -> xr.DataArray:
    r"""2-D kinematic frontogenesis, Sharman (2006) A9 / Bluestein (1993).

        F = −(1/|∇θ|) [ θx² ∂u/∂x + θy² ∂v/∂y
                        + θx θy (∂v/∂x + ∂u/∂y) ]                   (A9)

    rojak's bug: the cross-term coefficient uses ∂v/∂y instead of the
    correct ∂v/∂x. Units: K² m⁻² s⁻¹.
    """
    u   = _sel_level(ds, "eastward_wind", target_level)
    v   = _sel_level(ds, "northward_wind", target_level)
    t   = _sel_level(ds, "temperature", target_level)
    p   = float(target_level)

    theta = t * (1000.0 / p) ** (287.0 / 1004.0)

    theta_x, theta_y = _grad(theta)
    mag_grad = np.sqrt(theta_x ** 2 + theta_y ** 2)
    inv_mag = -xr.where(mag_grad != 0, 1.0 / mag_grad, 0.0)

    du_dx, du_dy = _grad(u)
    dv_dx, dv_dy = _grad(v)

    f2d = inv_mag * (
        theta_x ** 2 * du_dx
        + theta_y ** 2 * dv_dy
        + theta_x * theta_y * (dv_dx + du_dy)     # corrected cross term
    )
    f2d = f2d.rename("f2d")
    f2d.attrs.update({
        "long_name": "2-D frontogenesis",
        "units": "K2 m-2 s-1",
        "sharman_eq": "A9",
        "note": "cross-term uses dv_dx + du_dy (corrected)",
    })
    return f2d


# ---------------------------------------------------------------------------
# #20 (Q-F2D-5 replacement) — Isentropic frontogenesis, true Sharman A9
# ---------------------------------------------------------------------------
# The four candidate readings of Sharman A9. See FORMULA_AUDIT.md §4 for why
# more than one exists; ada/check_f2d_variants.py is what chooses between them.
F2D_VARIANTS = {
    "A": "+0.5 D/Dt[Q]        A10's algebra in theta coordinates (pre-2026-08-30)",
    "B": "-0.5 D/Dt[Q]        A9's leading minus sign applied",
    "C": "|0.5 D/Dt[Q]|       magnitude — matches the published distribution",
    "D": "-D/Dt[sqrt(Q)]      literal A9 including the |dv/dtheta|^-1 normalisation",
}

# DEFAULT CHANGED A -> C ON 2026-08-30. FORMULA_AUDIT.md 10.4 and the note in
# frontogenesis_isentropic below. Three independent lines put it beyond doubt,
# the decisive one being a figure rather than an equation:
#
#   1. Williams (2017) Fig. 1 plots the frontogenesis histogram on 0..300
#      x10^-9 m^2 s^-3 K^-2, anchored at exactly zero and decaying from an
#      18 % first bin. The same figure plots Negative Richardson on -300..0
#      and Colson-Panofsky on -45..-25, so Williams does show signed
#      diagnostics on their true negative ranges -- frontogenesis is simply
#      not one of them. The smooth decay from a finite first bin is the
#      density of |X| for a signed X; a clipped max(X,0) would put a ~50 %
#      point mass in that bin instead.
#   2. Measured p97/median: variant A 754, variant C 22.7, published 13.6.
#      Only C is in family with the other twenty diagnostics.
#   3. Sharman Table B1's units (m^2 s^-3 K^-2) require the un-normalised
#      form, which A, B and C all satisfy and D does not.
#
# This also DISSOLVES the sign question that opened the whole investigation:
# under a magnitude, A9's leading minus is irrelevant. That is presumably why
# Sharman's printed inconsistency between A9's two sides never mattered
# operationally.
#
# Anything computed before this date used A. The variant is recorded in each
# output's attributes, so a zarr can always be traced to the reading that
# produced it, and `--f2d-variant A` reproduces the old behaviour exactly.
F2D_DEFAULT_VARIANT = "C"


def frontogenesis_isentropic(ds: xr.Dataset, target_level: int = 200,
                             variant: str = F2D_DEFAULT_VARIANT) -> xr.DataArray:
    r"""Isentropic frontogenesis, Sharman (2006) A9 (Q-F2D-5).

        A9 = (1/2) D/Dt [ (du/dtheta)^2 + (dv/dtheta)^2 ]

    with D/Dt = d/dt + u d/dx + v d/dy the material derivative, and
    Q = (du/dtheta)^2 + (dv/dtheta)^2 throughout.

    WHICH FORM IS THIS? -- FOUR CANDIDATES, ONE FLAG (FORMULA_AUDIT.md §4)
    ---------------------------------------------------------------------
    The expression above is NOT what A9 prints. A9 reads

        F_theta  ~  -D/Dt [ Q ]^(1/2)
                 =  |dv/dtheta|^-1 [ du/dtheta D/Dt(du/dtheta)
                                   + dv/dtheta D/Dt(dv/dtheta) ]     (A9)

    and differs from what this function computed before 2026-08-29 in two
    independent ways:

      1. A NORMALISATION. A9 carries |dv/dtheta|^-1; the expression above does
         not. Since D/Dt(Q^1/2) = D/Dt(Q) / (2 Q^1/2), ours equals A9's times
         |dv/dtheta| -- a spatially VARYING field, not a constant. It changes
         RANKS, which is STATUS.md §5's category-1 error, the one a percentile
         calibration cannot absorb.

      2. A SIGN. A9's left-hand side carries a leading minus that its own
         right-hand side does not; the two sides as printed are not equal. The
         surrounding text ("...invoking continuity gives F_theta ~ -D/Dt...")
         puts the minus on the physics, and physically it belongs: on an
         isentropic surface a strengthening |grad theta| corresponds to a
         WEAKENING |dv/dtheta|, so frontogenesis is a decrease of the bracket.

    Sharman's Table B1 lists A9 in m^2 s^-3 K^-2, which are the units of the
    UN-normalised form and not of literal A9 (m s^-2 K^-1). That settles the
    algebra in favour of variants A/B/C -- and says nothing whatever about the
    sign, which is why the sign stayed open. Note also that the un-normalised
    form is structurally A10 (the constant-pressure midlevel form) evaluated in
    theta rather than p.

    Two further pieces of evidence, neither of which can be read off the paper:

      - Williams (2017) Table 2 p97 / Williams & Joshi (2013) median = 770/56.6
        = 13.6 for this diagnostic, from the same model, box, season and
        daily-mean sampling. A SIGNED material derivative cannot produce a
        finite ratio like that: D/Dt of a positive quantity in a statistically
        stationary atmosphere is centred on zero, so its median goes to zero
        and the ratio diverges. 13.6 is what a one-sided quantity looks like,
        and it fits the algebraic-degree pattern of the other 20 exactly.
      - Prosser (2023) Figure S5's frontogenesis panel is near-zero to slightly
        negative over the North Atlantic, where variant A gives the second
        largest positive trend of the 21.

    So rather than guess, all four readings are available here and
    `ada/check_f2d_variants.py` measured which one reproduces the published
    distribution.

    **RESOLVED 2026-08-30: the answer is C, and the default is now C.**
    Williams (2017) Figure 1 plots this diagnostic's histogram on 0..300,
    anchored at zero and decaying from an 18 % first bin, in the same figure
    where Negative Richardson runs -300..0 and Colson-Panofsky -45..-25. It is
    a magnitude, not a signed tendency. Measured p97/median confirms it:
    A 754, C 22.7, published 13.6. See the note on F2D_DEFAULT_VARIANT above
    and FORMULA_AUDIT.md 10.4.

    Args:
        variant: one of F2D_VARIANTS.
            "A"  +0.5 D/Dt[Q]        (default, unchanged behaviour)
            "B"  -0.5 D/Dt[Q]        A9's leading minus
            "C"  |0.5 D/Dt[Q]|       magnitude
            "D"  -D/Dt[sqrt(Q)]      literal A9, including the normalisation.
                 NOTE: variant D has DIFFERENT UNITS (m s^-2 K^-1, not
                 m^2 s^-3 K^-2), so its magnitude is not comparable with the
                 published tables even though its ranks are meaningful.

    Q-F2D-2/3 (Opus, resolved 2026-07-03) confirmed the Miller physical-
    space form above (frontogenesis_2d) is the WRONG diagnostic -- Prosser
    (2023) inherits the isentropic A9 via Williams (2017) -> Sharman
    (2006), confirmed independently by the citation chain and by units
    (m^2 s^-3 K^-2 is only derivable from A9; K^2 m^-2 s^-1 cannot produce
    it). frontogenesis_2d() is kept above for history/documentation but is
    no longer the target implementation for #20.

    du/dtheta, dv/dtheta: chain rule (du/dp)/(dtheta/dp) via
    theta_derivative_on_pressure_level -- NOT a full isentropic remap,
    computed directly on the existing 175/200/225 hPa levels, same bulk
    pattern as every other vertical derivative in this module.

    D/Dt: the FIRST diagnostic in this set needing an actual time
    derivative (every other diagnostic is purely spatial). Time
    derivative via xarray .differentiate("time") -- centered in the
    interior of the time dimension, one-sided at the first/last
    timestep. IMPORTANT CAVEAT verified against the 2-timestep
    (00Z/03Z) trial subset: with only 2 timesteps, .differentiate
    returns the SAME one-sided slope at both points -- there is no
    "interior" to center on. A true centered d/dt requires >= 3
    consecutive timesteps, i.e. a timestep strictly before AND after
    every evaluation time. See Q-DATA-1 / Q-F2D-5 data-requirements note.

    Units: m^2 s^-3 K^-2 for variants A/B/C; m s^-2 K^-1 for variant D.
    """
    if variant not in F2D_VARIANTS:
        raise ValueError(
            f"unknown f2d variant {variant!r}; expected one of "
            f"{sorted(F2D_VARIANTS)}. " +
            " | ".join(f"{k}: {d}" for k, d in F2D_VARIANTS.items())
        )

    u = ds["eastward_wind"]
    v = ds["northward_wind"]
    t = ds["temperature"]
    theta = _potential_temperature(t, t["pressure_level"])

    du_dtheta = theta_derivative_on_pressure_level(u, theta)
    dv_dtheta = theta_derivative_on_pressure_level(v, theta)
    q_field = du_dtheta ** 2 + dv_dtheta ** 2   # (du/dtheta)^2 + (dv/dtheta)^2, all levels/times

    q_lvl = q_field.sel(pressure_level=target_level)
    u_lvl = _sel_level(ds, "eastward_wind", target_level)
    v_lvl = _sel_level(ds, "northward_wind", target_level)

    # Variant D differentiates |dv/dtheta| = sqrt(Q) DIRECTLY rather than
    # forming D/Dt(Q)/(2 sqrt(Q)). Algebraically the same; numerically better,
    # and it avoids dividing by a sqrt(Q) that can be arbitrarily small in
    # low-shear cells, which would manufacture a spurious tail exactly where
    # the diagnostic is least meaningful.
    scalar = np.sqrt(q_lvl) if variant == "D" else q_lvl

    d_dt = scalar.differentiate("time", datetime_unit="s")
    d_dx, d_dy = _grad(scalar)
    material_derivative = d_dt + u_lvl * d_dx + v_lvl * d_dy

    if variant == "A":
        a9 = 0.5 * material_derivative
    elif variant == "B":
        a9 = -0.5 * material_derivative
    elif variant == "C":
        a9 = np.abs(0.5 * material_derivative)
    else:  # "D"
        a9 = -material_derivative

    a9 = a9.rename("f2d")
    a9.attrs.update({
        "long_name": "Isentropic frontogenesis (A9)",
        "units": "m s-2 K-1" if variant == "D" else "m2 s-3 K-2",
        "sharman_eq": "A9",
        "f2d_variant": variant,
        "f2d_variant_formula": F2D_VARIANTS[variant],
        "note": "Q-F2D-5: replaces Miller-form f2d as #20 target. "
                "du/dtheta,dv/dtheta via chain rule (du/dp)/(dtheta/dp); "
                "material derivative d/dt + V.grad; d/dt needs >=3 "
                "consecutive timesteps for a true centered difference "
                "(2-timestep data gives a one-sided estimate only). "
                "Variant recorded above -- see FORMULA_AUDIT.md 4 and "
                "ada/check_f2d_variants.py; A is the historical default and "
                "does NOT match A9 as printed.",
    })
    return a9


# ---------------------------------------------------------------------------
# #21 — NCSU1  (Sharman A36)
# ---------------------------------------------------------------------------
def ncsu1(ds: xr.Dataset, target_level: int = 200, ri_floor: float = 1e-5) -> xr.DataArray:
    r"""NCSU1 index, Sharman (2006) A36 / Kaplan et al. (2005).

        NCSU1 = [ 1 / max(Ri, 1e-5) ]
                · max( u ∂u/∂x + v ∂v/∂y, 0 )
                · |∇ζ|                                             (A36)

    Two things must match rojak exactly for this to be correct:
      1. du/dx and dv/dy are derivatives of VECTOR COMPONENTS, not scalars,
         so they need rojak's vector_derivatives() (which applies extra
         map-projection correction terms for a sphere), not the plain
         scalar-field _grad().
      2. Ri is floored from BELOW: max(Ri, floor), matching rojak's own
         `ri.clip(min=RI_THRESHOLD)` -- not an abs-value floor.

    Units: s⁻³.
    """
    u    = _sel_level(ds, "eastward_wind", target_level)
    v    = _sel_level(ds, "northward_wind", target_level)
    zeta = _sel_level(ds, "vorticity", target_level)

    # Projection-corrected vector-component derivatives (see docstring).
    vd = vector_derivatives(u, v, "deg",
                            components=[VelocityDerivative.DU_DX, VelocityDerivative.DV_DY])
    du_dx = vd[VelocityDerivative.DU_DX]
    dv_dy = vd[VelocityDerivative.DV_DY]

    advection = (u * du_dx + v * dv_dy).clip(min=0)

    dzeta_dx, dzeta_dy = _grad(zeta)
    grad_zeta = np.sqrt(dzeta_dx ** 2 + dzeta_dy ** 2)

    ri = richardson(ds, negative=False).sel(pressure_level=target_level)
    ri_safe = ri.clip(min=ri_floor)     # max(Ri, floor) -- matches rojak exactly

    out = (advection * grad_zeta) / ri_safe
    out = out.rename("ncsu1")
    out.attrs.update({
        "long_name": "NCSU1 index",
        "units": "s-3",
        "sharman_eq": "A36",
        "note": "corrected Ri; projection-corrected du_dx,dv_dy; floor=max(Ri,1e-5)",
    })
    return out


# ---------------------------------------------------------------------------
# #12 — Relative vorticity advection magnitude  (not in rojak)
# ---------------------------------------------------------------------------
def rva(ds: xr.Dataset, target_level: int = 200) -> xr.DataArray:
    r"""Magnitude of relative vorticity advection.

        RVA = | u ∂ζ/∂x + v ∂ζ/∂y |

    Not provided by rojak. Computed with geospatial (metre-based)
    derivatives. Units: s⁻².
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


# ---------------------------------------------------------------------------
# #15 — Brown energy dissipation rate (Brown2)  (Sharman A14)
# ---------------------------------------------------------------------------
def brown2(ds: xr.Dataset, brown1: xr.DataArray) -> xr.DataArray:
    r"""Brown Index 2 / Brown EDR, Sharman (2006) A14.

        ε_brown = (1/24) Φ Sv²                                     (A14)

    rojak transcribes A14 faithfully; the only wrinkle is dimensional --
    A14 yields s⁻³, whereas W&J Table 1 lists J kg⁻¹ s⁻¹ = m² s⁻³, implying
    an unspecified length² factor. We keep native s⁻³ and treat this as a
    rank-only diagnostic (Prosser's percentile-threshold method is
    invariant to a constant multiplicative L²); no length scale invented.
    """
    u = ds["eastward_wind"]
    v = ds["northward_wind"]

    du_dz = altitude_derivative_on_pressure_level(u, ds["geopotential"])
    dv_dz = altitude_derivative_on_pressure_level(v, ds["geopotential"])
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
def compute_all_21(
    catdata: CATData, target_level: int = 200,
    f2d_variant: str = F2D_DEFAULT_VARIANT,
) -> tuple[dict[str, xr.DataArray], list[DiagnosticFailure]]:
    """Compute all 21 diagnostics: 14 from rojak, 7 hand-written.

    TWO THINGS THIS FUNCTION DOES THAT compute_rojak_diagnostics DOES NOT
    ---------------------------------------------------------------------
    `compute_rojak_diagnostics` is a faithful passthrough of rojak and is
    deliberately left that way, because `tests/test_analytic.py` and
    `tests/report_errors.py` verify rojak against manufactured solutions
    through it. This function is the PRODUCTION entry point -- everything that
    reaches disk goes through here -- so the two corrections below live here
    and nowhere else:

      1. `deformation` is un-squared. See the block below.
      2. `f2d`'s variant is selectable. See F2D_VARIANTS.

    Returns ({short_key: xr.DataArray}, [DiagnosticFailure, ...]) --
    Q-INTEG-4: the failures list is empty in the normal case, and carries
    one entry per rojak diagnostic that failed and was replaced with an
    all-NaN placeholder for this chunk (see compute_rojak_diagnostics).
    Callers looping over chunks (the 1979-2020 driver) should accumulate
    this list across chunks into a run-level failure summary -- that
    accumulation happens outside this function, which only sees one chunk.

    Diagnostics that still carry a pressure_level dimension are reduced to
    `target_level` at the end.
    """
    ds = catdata._dataset

    out, failures = compute_rojak_diagnostics(catdata)

    # -----------------------------------------------------------------------
    # #8 DEFORMATION -- persist DEF, not DEF^2   (FORMULA_AUDIT.md §5)
    # -----------------------------------------------------------------------
    # rojak's DEF diagnostic returns the SQUARE of total deformation
    # (`DeformationSquared`), while Sharman A17 and every published table
    # define DEF = (D_SH^2 + D_ST^2)^(1/2). Until 2026-08-29 the square root
    # was applied only in 3_pipeline.py's PRETRANSFORM, i.e. only on the path
    # that builds the W&J comparison table -- so every zarr written by
    # ada/diagnostics_global.py held DEF^2 under the name `deformation`.
    #
    # WHY THIS WAS INVISIBLE, AND WHY IT STILL MATTERS. Squaring is strictly
    # increasing on a non-negative field, so by STATUS.md §5.5 the exceedance
    # field is IDENTICAL and nothing in the §11 calibration check or the §12
    # trend check is affected -- which is exactly why it survived this long.
    # It stops being harmless the moment magnitudes are used instead of ranks:
    # STATUS.md §6 phase 5 plans peaks-over-threshold and GPD fitting on
    # retained excesses, and a GPD fitted to DEF^2 is not a GPD fitted to DEF
    # (the tail index doubles under squaring). Fixing it after the production
    # run means re-deriving ~287 GB; fixing it here costs nothing.
    #
    # Note this does NOT touch TI1/TI2/NGM1/NGM2, which never had the problem:
    # rojak hands those `CATData.total_deformation()`, which is
    # `magnitude_of_vector(..., is_squared=False)` -- already un-squared.
    # Verified by reading rojak at the pinned rev 25b8685, and pinned by
    # tests/test_analytic.py::TestDeformationFamily.
    if "deformation" in out:
        _def_squared = out["deformation"]
        _def = np.sqrt(np.abs(_def_squared)).rename("deformation")
        _def.attrs.update(dict(_def_squared.attrs))     # np.sqrt drops attrs
        _def.attrs.update({
            "long_name": "Total deformation",
            "units": "s-1",
            "sharman_eq": "A17",
            "note": "sqrt of rojak's DeformationSquared, so this is DEF not "
                    "DEF^2 (FORMULA_AUDIT.md 5). abs() guards float noise "
                    "only -- rojak's value is non-negative by construction. "
                    "Zarr written before 2026-08-29 holds DEF^2 under this "
                    "name; exceedance fields are unaffected, magnitudes are "
                    "not.",
        })
        out["deformation"] = _def

    brown1 = DiagnosticFactory(catdata).create(TurbulenceDiagnostics.BROWN1).computed_value
    builders = {
        "negative_richardson": lambda: richardson(ds, negative=True),
        "colson_panofsky":     lambda: colson_panofsky(ds),
        "ubf":                 lambda: ubf(ds, target_level=target_level),
        "f2d":                 lambda: frontogenesis_isentropic(ds, target_level=target_level, variant=f2d_variant),  # Q-UNITS-1 caught this: Q-F2D-5 added the isentropic A9 impl and wired it into 4_verify.py's hand-dict, but never updated THIS dispatch (used by 3_pipeline.py/cat_pipeline.py) -- was still silently calling the old Miller-form frontogenesis_2d()
        "ncsu1":               lambda: ncsu1(ds, target_level=target_level),
        "rva_magnitude":       lambda: rva(ds, target_level=target_level),
        "brown2":              lambda: brown2(ds, brown1),
    }
    for key, builder in builders.items():
        out[key] = builder().rename(key)

    for k, da in list(out.items()):
        if "pressure_level" in da.dims:
            try:
                out[k] = da.sel(pressure_level=target_level)
            except KeyError:
                # Q-INTEG-3 (Bug 2): do NOT silently substitute a different
                # real level's value here. `da.isel(pressure_level=len//2)`
                # used to fall back to whatever level happened to sit at the
                # midpoint of THIS diagnostic's own (possibly level-reduced,
                # e.g. colson_panofsky's [200, 175]) coordinate -- a real,
                # valid value at the WRONG level, returned with no warning.
                # That is strictly worse than a crash: every other
                # config-driven level request (ubf, f2d, ncsu1, rva_magnitude,
                # ...) shares this dispatch loop, so a silent substitution
                # here is a landmine for all of them, not just CP.
                #
                # Deliberately scoped to fail loudly instead of guessing.
                # This does NOT add arbitrary-target-level support to any
                # individual diagnostic (e.g. colson_panofsky() itself is
                # unchanged) -- whether CP specifically should support
                # arbitrary levels is an open scope question (reopened
                # pressure-level question, Q-RI-1 / Q-MODEL-3-CHECK-B) and is
                # deliberately NOT decided here.
                available = list(da.pressure_level.values)
                raise KeyError(
                    f"compute_all_21: diagnostic '{k}' has no value at "
                    f"target_level={target_level!r} -- its pressure_level "
                    f"coordinate only contains {available}. Refusing to "
                    f"silently substitute a different level's value. If "
                    f"{target_level!r} is a level this diagnostic is "
                    f"expected to support, that is a scope decision "
                    f"(see Q-RI-1 / Q-MODEL-3-CHECK-B in STATUS), not a bug "
                    f"to paper over here."
                ) from None
    return out, failures