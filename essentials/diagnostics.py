"""
diagnostics.py -- Stage 2: one ERA5 month in, 21 CAT diagnostics at 200 hPa out.

Reading guide
-------------
    1. open_era5 / prepare          read a GRIB file and put it in one standard layout
    2. grid_geometry                grid spacing (m) and map-scale factors
    3. operators                    d/dx, d/dy, wind derivatives, d/dp, d/dz, theta, f
    4. flow_fields                  vorticity, divergence, PV: archived or computed
    5. the 21 diagnostics           one short function each, formula in the docstring
    6. compute_21 / compute_month   all 21 for one file, in 1-day chunks

Only numpy, xarray and pyproj are used. The finite differences and the
spherical corrections are written out here, and they reproduce rojak @ 25b8685
(which copies MetPy) exactly. tests/compare_with_original.py checks this against the
original 2_diagnostics.py on real ERA5 data (bit-for-bit, both conventions).

Notation used throughout
------------------------
    u, v       wind components [m s-1]         T   temperature [K]
    z          geopotential Phi = g*Z [m2 s-2]  p   pressure [hPa]
    zeta       relative vorticity [s-1]         f   Coriolis parameter 2*Omega*sin(phi)
    delta      horizontal divergence [s-1]      theta  potential temperature [K]
    phi, lam   latitude, longitude
    x, y       east and north distance on the Earth
    k_x, k_y   map-scale factors (see grid_geometry)

Every array is (latitude, longitude, time, pressure_level) until the very end,
where each diagnostic is evaluated at 200 hPa. The three levels 175/200/225 hPa
exist only to give the vertical derivatives a centred stencil at 200 hPa.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import xarray as xr
from pyproj import CRS, Geod, Proj

import config as C


# ===========================================================================
# 1. Reading and preparing one ERA5 file
# ===========================================================================
RENAME = {
    "isobaricInhPa": "pressure_level",
    "u": "eastward_wind", "v": "northward_wind", "t": "temperature",
    "z": "geopotential", "d": "divergence_of_wind", "vo": "vorticity",
    "pv": "potential_vorticity",
}


def open_era5(path) -> xr.Dataset:
    """Open one ERA5 GRIB (or NetCDF) file without writing a cfgrib .idx file."""
    path = Path(path)
    if path.suffix in (".grib", ".grb", ".grb2"):
        return xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
    return xr.open_dataset(path)


def icao_altitude(p_hpa):
    """ICAO standard-atmosphere height [m] of a pressure level.

        z(p) = T0/L * (1 - (p/p0)^(R L / g))

    Only used for the Colson-Panofsky length scale, where it enters as a
    constant (it depends on p alone).
    """
    p = np.asarray(p_hpa, dtype=float) * 100.0
    p0, T0, L, g, R = 101325.0, 288.15, 0.0065, 9.80665, 287.0531
    return T0 / L * (1.0 - (p / p0) ** ((R * L) / g))


def prepare(ds: xr.Dataset) -> xr.Dataset:
    """Give a raw ERA5 dataset the one layout everything below assumes.

    - standard variable names (u -> eastward_wind, ...)
    - longitudes in [-180, 180), ascending (global files come as 0..360)
    - an ICAO `altitude` coordinate on the pressure levels
    - dimension order (latitude, longitude, time, pressure_level), time sorted
    """
    if "time" not in ds.dims and "valid_time" in ds.dims:     # newer CDS NetCDF layout
        ds = ds.rename({"valid_time": "time"})
    ds = ds.rename({k: v for k, v in RENAME.items() if k in ds.variables and k != v})

    if float(ds.longitude.max()) > 180.0:
        ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180)).sortby("longitude")
    ds = ds.assign_coords(altitude=("pressure_level", icao_altitude(ds.pressure_level.values)))
    ds = ds.transpose("latitude", "longitude", "time", "pressure_level").sortby("time")

    # Two silent assumptions made loud (audit F8):
    levels = np.asarray(ds["pressure_level"].values, dtype=float)
    if levels.size >= 3 and not np.allclose(np.diff(levels), np.diff(levels)[0], rtol=1e-9, atol=1e-9):
        raise ValueError(f"pressure levels {levels.tolist()} are not evenly spaced; "
                         "the wrapped wind-direction difference (Endlich) needs them to be")
    if float(ds["latitude"].max()) <= 4.0:
        raise ValueError("domain lies entirely south of 4N; the degree/radian test in the "
                         "Coriolis parameter would misread it (audit 7.10)")
    return ds


# ===========================================================================
# 2. Geometry: grid spacing and map-scale factors
# ===========================================================================
def grid_geometry(latitude, longitude, dy_method: str) -> dict:
    """Grid spacings and map-scale factors for a regular lat-lon grid.

    The grid is treated as a map (plate carree) with coordinates
        X = a * lam,        Y = a * phi         (a = WGS84 equatorial radius)
    so a finite difference along the grid is a derivative with respect to X or
    Y. A physical derivative on the ellipsoid then follows by multiplying by a
    map-scale factor:

        d/dx = k_x d/dX,   k_x = a / (N(phi) cos phi)      (parallel scale)
        d/dy = k_y d/dY,   k_y = a / M(phi)                (meridional scale)

    with N and M the prime-vertical and meridional radii of curvature,
        N = a / sqrt(1 - e^2 sin^2 phi),   M = a (1 - e^2) / (1 - e^2 sin^2 phi)^(3/2).

    k_x and k_y are taken from PROJ (`+proj=latlon`), exactly as rojak/MetPy
    do. PROJ evaluates them numerically; they agree with the closed forms above
    to ~1e-11 (tests/test_operators.py), and at the poles PROJ returns a large
    finite k_x (~1e5) instead of infinity.

    dx = a * dlam between neighbouring columns (identical to the geodesic
    distance along the equator that rojak uses).

    dy depends on the convention -- this is audit fix F1:
      "geodesic"  dy = the true meridian arc between rows (~ M dphi), which is
                  what rojak/MetPy use. Multiplying that by k_y = a/M applies
                  the ellipsoid correction twice: a latitude-dependent error of
                  +0.42 % at 30N to -0.27 % at 75N in every d/dy.
      "map"       dy = a * dphi, the map distance, so k_y is applied once.
                  Error against exact ellipsoidal gradients: 2e-7.
    dy is negative because ERA5 latitudes run north to south.

    Also returned: dk_x/dY and dk_y/dX, the rate at which the scale factors
    change across the grid. They are the curvature terms in wind_derivatives.
    """
    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)

    dx = C.A_WGS84 * np.deg2rad(np.abs(np.diff(lon)))
    if dy_method == "map":
        dy = C.A_WGS84 * np.deg2rad(np.diff(lat))
    elif dy_method == "geodesic":
        zeros = np.zeros_like(lat)
        azimuth, _, dy = Geod(ellps="WGS84").inv(zeros[:-1], lat[:-1], zeros[1:], lat[1:])
        dy[(azimuth < -90.0) | (azimuth > 90.0)] *= -1      # southward step -> negative
    else:
        raise ValueError(f"dy_method must be 'map' or 'geodesic', got {dy_method!r}")

    lon_grid, lat_grid = np.meshgrid(lon, lat)
    factors = Proj(CRS("+proj=latlon")).get_factors(lon_grid, lat_grid)
    coords = {"latitude": latitude, "longitude": longitude}
    k_x = xr.DataArray(factors.parallel_scale, dims=("latitude", "longitude"), coords=coords)
    k_y = xr.DataArray(factors.meridional_scale, dims=("latitude", "longitude"), coords=coords)

    geo = {"dx": dx, "dy": dy, "k_x": k_x, "k_y": k_y}
    geo["dkx_dY"] = grid_difference(k_x, dy, "latitude")
    geo["dky_dX"] = grid_difference(k_y, dx, "longitude")
    return geo


# ===========================================================================
# 3. Operators
# ===========================================================================
def grid_difference(f: xr.DataArray, spacing_m: np.ndarray, dim: str) -> xr.DataArray:
    """d f / d(grid distance) along `dim`: np.gradient on the cumulative distance.

    Second-order centred in the interior,
        f'_i = [h_-^2 f_{i+1} + (h_+^2 - h_-^2) f_i - h_+^2 f_{i-1}] / [h_- h_+ (h_- + h_+)],
    first-order one-sided at the two edge rows/columns. (For equal spacings this
    is (f_{i+1} - f_{i-1}) / 2h.) Keeps the input dtype, like np.gradient.
    """
    distance = np.cumsum(np.insert(spacing_m, 0, [0]))
    return f.copy(data=np.gradient(f.values, distance, axis=f.get_axis_num(dim)))


def d_dx(f: xr.DataArray, geo: dict) -> xr.DataArray:
    """Eastward derivative of a SCALAR field: df/dx = k_x * df/dX."""
    return grid_difference(f, geo["dx"], "longitude") * geo["k_x"]


def d_dy(f: xr.DataArray, geo: dict) -> xr.DataArray:
    """Northward derivative of a SCALAR field: df/dy = k_y * df/dY."""
    return grid_difference(f, geo["dy"], "latitude") * geo["k_y"]


def wind_derivatives(u: xr.DataArray, v: xr.DataArray, geo: dict) -> dict:
    """The four horizontal derivatives of a VECTOR field (u, v) on the sphere.

    Differentiating a vector component is not the same as differentiating a
    scalar: the local east/north unit vectors rotate as you move, which adds
    curvature terms (MetPy's `vector_derivative`):

        du/dx = k_x du/dX - v * (k_y/k_x) dk_x/dY
        du/dy = k_y du/dY + v * (k_x/k_y) dk_y/dX
        dv/dx = k_x dv/dX + u * (k_y/k_x) dk_x/dY
        dv/dy = k_y dv/dY - u * (k_x/k_y) dk_y/dX

    On a sphere (k_y/k_x) dk_x/dY = tan(phi)/R, so e.g. a purely zonal flow
    with no shear still has vorticity u tan(phi)/R. Without these terms the
    vorticity is 100 % wrong (Stokes' theorem test, STATUS 4g).
    dk_y/dX is zero on this grid, so the last terms vanish in practice.

    Standing rule: any derivative OF u or v, and any divergence of a vector,
    goes through this function -- never through d_dx / d_dy.
    """
    k_x, k_y = geo["k_x"], geo["k_y"]
    x_curv = (k_y / k_x) * geo["dkx_dY"]
    y_curv = (k_x / k_y) * geo["dky_dX"]
    return {
        "du_dx": k_x * grid_difference(u, geo["dx"], "longitude") - v * x_curv,
        "du_dy": k_y * grid_difference(u, geo["dy"], "latitude") + v * y_curv,
        "dv_dx": k_x * grid_difference(v, geo["dx"], "longitude") + u * x_curv,
        "dv_dy": k_y * grid_difference(v, geo["dy"], "latitude") - u * y_curv,
    }


def d_dp(f: xr.DataArray) -> xr.DataArray:
    """df/dp on the pressure levels [per hPa]: centred at 200 hPa,
    (f(175) - f(225)) / (175 - 225), one-sided at 175 and 225."""
    return f.differentiate("pressure_level")


def d_dz(f: xr.DataArray, z: xr.DataArray) -> xr.DataArray:
    """Vertical derivative with respect to height, on pressure levels.

        df/dz = (df/dp) / (dZ/dp) = g * (df/dp) / (dPhi/dp)       (Phi = g Z)
    """
    return C.G * (d_dp(f) / d_dp(z))


def potential_temperature(T: xr.DataArray) -> xr.DataArray:
    """theta = T / (p / p0)^kappa,  p0 = 1000 hPa, kappa = R_d/c_p = 2/7."""
    return T / ((T["pressure_level"] / C.P0_HPA) ** C.KAPPA)


def coriolis(latitude: xr.DataArray) -> xr.DataArray:
    """f = 2 Omega sin(phi)."""
    return 2 * C.OMEGA * np.sin(np.deg2rad(latitude))


def wind_direction(u: xr.DataArray, v: xr.DataArray) -> xr.DataArray:
    """Meteorological wind direction psi in [0, 2pi): the direction the wind
    blows FROM, clockwise from north. psi = pi/2 - atan2(-v, -u)."""
    psi = np.pi / 2 - np.arctan2(-v, -u)
    psi = xr.where(psi <= 0.0, psi + 2 * np.pi, psi)
    return xr.where((u == 0) & (v == 0), 0, psi)


def _shorter_arc(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """|a - b| for two angles, taken the short way round: min(|a-b|, 2pi-|a-b|)."""
    diff = np.abs(a - b)
    other_way = (2 * np.pi) - diff
    return np.where(diff < other_way, diff, other_way)


def d_angle_dp(psi: xr.DataArray) -> xr.DataArray:
    """|d psi / dp| for an angle, with the difference wrapped to the shorter arc.

    Same stencil as d_dp (centred in the middle, one-sided at the ends). Note
    the result is the ABSOLUTE turning divided by the (negative) level spacing,
    so only its magnitude is meaningful. Requires evenly spaced levels.
    """
    p = psi["pressure_level"].values
    step = np.diff(p)[0]
    a = np.moveaxis(psi.values, psi.get_axis_num("pressure_level"), 0)
    out = np.empty_like(a)
    out[1:-1] = _shorter_arc(a[2:], a[:-2]) / (2.0 * step)
    out[0] = _shorter_arc(a[1], a[0]) / step
    out[-1] = _shorter_arc(a[-1], a[-2]) / step
    return psi.copy(data=np.moveaxis(out, 0, psi.get_axis_num("pressure_level")))


# ===========================================================================
# 4. Vorticity, divergence and PV: archived or computed (audit fix F12)
# ===========================================================================
def flow_fields(ds: xr.Dataset, wd: dict, source: str):
    """Return (zeta, delta, PV) on all levels.

    "archived": ERA5's own vo, d and pv. These are the IFS model's spectral
                variables (u, v are derived from them), i.e. closest to the model.
                pv is the full Ertel PV.
    "computed": rebuilt from u, v, T, z -- the four fields Prosser (2023) used:
                    zeta  = dv/dx - du/dy
                    delta = du/dx + dv/dy
                    PV    = -g (zeta + f) dtheta/dp          (Sharman 2006 A18)
                with the wind derivatives of section 3. A18 keeps only the
                vertical term of Ertel PV and is 100x the SI value because p is
                in hPa (a constant factor, harmless for percentile thresholds).
    """
    if source == "archived":
        return ds["vorticity"], ds["divergence_of_wind"], ds["potential_vorticity"]
    if source != "computed":
        raise ValueError(f"fields must be 'archived' or 'computed', got {source!r}")
    zeta = wd["dv_dx"] - wd["du_dy"]
    delta = wd["du_dx"] + wd["dv_dy"]
    theta = potential_temperature(ds["temperature"])
    pv = -C.G * (zeta + coriolis(zeta["latitude"])) * d_dp(theta)
    return zeta, delta, pv


# ===========================================================================
# 5. The 21 diagnostics
# ===========================================================================
# Every function takes the dict `s` built in compute_21 (the prepared fields and
# the shared building blocks) and returns an xarray DataArray. Numbers in
# brackets are the Williams & Joshi (2013) Table 1 numbers; "A.." are equations
# in Sharman et al. (2006) Appendix A. Functions that differ between the two
# conventions say so.

def building_blocks(ds: xr.Dataset, geo: dict, convention: dict) -> dict:
    """Fields shared by several diagnostics, computed once."""
    u, v, T, z = ds["eastward_wind"], ds["northward_wind"], ds["temperature"], ds["geopotential"]
    wd = wind_derivatives(u, v, geo)
    zeta, delta, pv = flow_fields(ds, wd, convention["fields"])

    D_sh = wd["dv_dx"] + wd["du_dy"]           # shearing deformation
    D_st = wd["du_dx"] - wd["dv_dy"]           # stretching deformation
    du_dz, dv_dz = d_dz(u, z), d_dz(v, z)
    theta = potential_temperature(T)
    N2 = (C.G / theta) * d_dz(theta, z)        # Brunt-Vaisala frequency squared (A2)
    return dict(ds=ds, geo=geo, conv=convention, u=u, v=v, T=T, z=z, wd=wd,
                zeta=zeta, delta=delta, pv=pv, D_sh=D_sh, D_st=D_st,
                DEF=np.hypot(D_sh, D_st),      # total deformation (A17)
                du_dz=du_dz, dv_dz=dv_dz,
                Sv=np.hypot(du_dz, dv_dz),     # vertical wind shear (A3)
                Sv2=du_dz ** 2 + dv_dz ** 2,
                theta=theta, N2=N2, lev=C.TARGET_LEVEL_HPA)


def _at(da: xr.DataArray, s: dict) -> xr.DataArray:
    return da.sel(pressure_level=s["lev"])


def magnitude_pv(s):
    """[1] |PV|. Archived: |full Ertel PV|. Computed: |-g (zeta+f) dtheta/dp| (A18)."""
    return np.abs(s["pv"])


def colson_panofsky(s):
    """[2] Colson-Panofsky index (A4), in its reduced form

        CP = lambda^2 (Sv^2 - N^2 / Ri_crit),       Ri_crit = 0.5

    which equals lambda^2 Sv^2 (1 - Ri/Ri_crit) without forming Ri. lambda is the
    ICAO height difference between 200 and 225 hPa: a constant.
    """
    alt = s["ds"]["altitude"]
    lam = alt.diff("pressure_level", label="upper").sel(pressure_level=s["lev"])
    return np.square(lam) * (_at(s["Sv2"], s) - _at(s["N2"], s) / C.RI_CRIT)


def brown1(s):
    """[3] Brown index (A13): Phi = sqrt(0.3 zeta_a^2 + D_sh^2 + D_st^2), zeta_a = zeta + f."""
    zeta_a = s["zeta"] + coriolis(s["zeta"]["latitude"])
    return np.sqrt(0.3 * np.square(zeta_a) + np.square(s["D_sh"]) + np.square(s["D_st"]))


def temperature_gradient(s):
    """[4] |grad_H T| = sqrt((dT/dx)^2 + (dT/dy)^2)   (A23)."""
    return np.hypot(d_dx(s["T"], s["geo"]), d_dy(s["T"], s["geo"]))


def horizontal_divergence(s):
    """[5] |delta|, delta = du/dx + dv/dy   (A33). Archived ERA5 d, or computed (convention)."""
    return np.abs(s["delta"])


def vertical_wind_shear(s):
    """[6] Sv = |dV/dz| = sqrt((du/dz)^2 + (dv/dz)^2)   (A3)."""
    return s["Sv"]


def endlich(s):
    """[7] Endlich index (A25): |V| * |d psi / dz|, psi = wind direction.

    "angle":      |d psi/dz| = g * |wrapped psi difference / dp| / (dPhi/dp)
    "components": |d psi/dz| = |u dv/dz - v du/dz| / |V|^2  (exact for psi = atan2)
    The two are different estimators on a 50 hPa stencil (audit F10).
    """
    u, v, z = s["u"], s["v"], s["z"]
    if s["conv"]["endlich"] == "angle":
        dpsi_dz = C.G * d_angle_dp(wind_direction(u, v)) / d_dp(z)
        return np.hypot(u, v) * np.abs(dpsi_dz)
    speed = np.hypot(u, v)
    return np.abs(u * s["dv_dz"] - v * s["du_dz"]) / xr.where(speed > 0, speed, np.nan)


def deformation(s):
    """[8] DEF = sqrt(D_sh^2 + D_st^2)  (A17), D_sh = dv/dx + du/dy, D_st = du/dx - dv/dy.

    Written as sqrt(|DEF^2|) because that is how the original pipeline undid
    rojak's squared output; kept so the stored numbers are bit-identical.
    """
    return np.sqrt(np.abs(np.square(s["DEF"])))


def wind_speed(s):
    """[9] |V| = sqrt(u^2 + v^2)   (A24)."""
    return np.hypot(s["u"], s["v"])


def ngm2(s):
    """[10] NGM2 = |dT/dz| * DEF   (A29)."""
    return np.abs(d_dz(s["T"], s["z"])) * s["DEF"]


def negative_richardson(s):
    """[11] -Ri, with the gradient Richardson number Ri = N^2 / Sv^2  (A1)."""
    return -(s["N2"] / s["Sv2"])


def rva_magnitude(s):
    """[12] |relative vorticity advection| = |u dzeta/dx + v dzeta/dy|
    (Williams & Joshi 2013; not in Sharman Appendix A)."""
    zeta = _at(s["zeta"], s)
    return np.abs(_at(s["u"], s) * d_dx(zeta, s["geo"]) + _at(s["v"], s) * d_dy(zeta, s["geo"]))


def ubf(s):
    """[13] Unbalanced flow (A30, Koch & Caracena 2002): the residual of the
    nonlinear balance equation

        UBF = | lap(Phi) - 2 J(u,v) - f zeta + beta u |

    lap(Phi) = divergence of grad(Phi)   (a vector divergence: wind_derivatives)
    J(u,v)   = du/dx dv/dy - du/dy dv/dx
    beta     = df/dy = 2 Omega cos(phi) / R
    """
    geo, lev = s["geo"], s["lev"]
    u, zeta, phi = _at(s["u"], s), _at(s["zeta"], s), _at(s["z"], s)
    lat_rad = np.deg2rad(s["ds"]["latitude"])
    f = (2 * C.OMEGA * np.sin(lat_rad)).broadcast_like(u.isel(time=0))
    beta = (2 * C.OMEGA * np.cos(lat_rad) / C.R_EARTH).broadcast_like(u.isel(time=0))

    grad_phi = wind_derivatives(d_dx(phi, geo), d_dy(phi, geo), geo)
    laplacian = grad_phi["du_dx"] + grad_phi["dv_dy"]
    wd = {k: w.sel(pressure_level=lev) for k, w in s["wd"].items()}
    jacobian = wd["du_dx"] * wd["dv_dy"] - wd["du_dy"] * wd["dv_dx"]

    return np.abs(laplacian - 2 * jacobian - f * zeta + beta * u)


def nva(s):
    """[14] Negative absolute vorticity advection (A37):
        NVA = max( -(u d(zeta+f)/dx + v d(zeta+f)/dy), 0 )."""
    zeta_a = s["zeta"] + coriolis(s["zeta"]["latitude"])
    advection_x = s["u"] * d_dx(zeta_a, s["geo"])
    advection_y = s["v"] * d_dy(zeta_a, s["geo"])
    return (-advection_x - advection_y).clip(min=0)


def brown2(s, brown1_value):
    """[15] Brown energy dissipation rate (A14): eps = Phi * Sv^2 / 24.
    Units s-3 (the published m2 s-3 implies an unstated length^2; constant)."""
    return (1.0 / 24.0) * brown1_value * s["Sv2"]


def vorticity_squared(s):
    """[16] zeta^2   (A21). Archived or computed zeta (convention)."""
    return np.square(s["zeta"])


def ti1(s):
    """[17] Ellrod TI1 = Sv * DEF   (A15)."""
    return s["Sv"] * s["DEF"]


def ngm1(s):
    """[18] NGM1 = |V| * DEF   (A28)."""
    return np.hypot(s["u"], s["v"]) * s["DEF"]


def ti2(s):
    """[19] Ellrod TI2 = Sv * (DEF + CVG), convergence CVG = -delta   (A16)."""
    return s["Sv"] * (s["DEF"] + (-s["delta"]))


F2D_VARIANTS = {
    "A": "+0.5 D/Dt[Q]     signed (Williams & Storer 2022 Eq. 3)",
    "B": "-0.5 D/Dt[Q]     A9's leading minus",
    "C": "|0.5 D/Dt[Q]|    magnitude",
    "D": "-D/Dt[sqrt(Q)]   literal A9 incl. normalisation",
}


def f2d(s):
    """[20] Isentropic frontogenesis (A9):

        F = 1/2 * D/Dt [ Q ],   Q = (du/dtheta)^2 + (dv/dtheta)^2
        D/Dt = d/dt + u d/dx + v d/dy          (material derivative at 200 hPa)
        du/dtheta = (du/dp) / (dtheta/dp)      (chain rule on the pressure levels)

    d/dt is centred in time (3-hourly), one-sided at the first and last step of
    the file. The variant (convention) decides the sign/magnitude reading.
    """
    variant = s["conv"]["f2d_variant"]
    u, v, theta = s["u"], s["v"], s["theta"]
    Q = (d_dp(u) / d_dp(theta)) ** 2 + (d_dp(v) / d_dp(theta)) ** 2
    scalar = _at(Q, s)
    if variant == "D":
        scalar = np.sqrt(scalar)
    d_dt = scalar.differentiate("time", datetime_unit="s")
    material = d_dt + _at(u, s) * d_dx(scalar, s["geo"]) + _at(v, s) * d_dy(scalar, s["geo"])
    if variant == "A":
        return 0.5 * material
    if variant == "B":
        return -0.5 * material
    if variant == "C":
        return np.abs(0.5 * material)
    if variant == "D":
        return -material
    raise ValueError(f"unknown f2d variant {variant!r}; choose from {sorted(F2D_VARIANTS)}")


def ncsu1(s):
    """[21] NCSU1 (A36, Kaplan et al. 2005):

        NCSU1 = max(u du/dx + v dv/dy, 0) * |grad zeta| / max(Ri, 1e-5)
    """
    geo = s["geo"]
    u, v, zeta = _at(s["u"], s), _at(s["v"], s), _at(s["zeta"], s)
    du_dx, dv_dy = _at(s["wd"]["du_dx"], s), _at(s["wd"]["dv_dy"], s)
    advection = (u * du_dx + v * dv_dy).clip(min=0)
    grad_zeta = np.sqrt(d_dx(zeta, geo) ** 2 + d_dy(zeta, geo) ** 2)
    ri = _at(s["N2"] / s["Sv2"], s)
    return (advection * grad_zeta) / ri.clip(min=C.RI_FLOOR)


# ===========================================================================
# 6. All 21 for one file
# ===========================================================================
def compute_21(ds: xr.Dataset, convention: dict) -> dict:
    """All 21 diagnostics at 200 hPa for one prepared dataset."""
    geo = grid_geometry(ds["latitude"], ds["longitude"], convention["dy"])
    s = building_blocks(ds, geo, convention)
    b1 = brown1(s)
    out = {
        "magnitude_pv": magnitude_pv(s),
        "colson_panofsky": colson_panofsky(s),
        "brown1": b1,
        "temperature_gradient": temperature_gradient(s),
        "horizontal_divergence": horizontal_divergence(s),
        "vertical_wind_shear": vertical_wind_shear(s),
        "endlich": endlich(s),
        "deformation": deformation(s),
        "wind_speed": wind_speed(s),
        "ngm2": ngm2(s),
        "negative_richardson": negative_richardson(s),
        "rva_magnitude": rva_magnitude(s),
        "ubf": ubf(s),
        "nva": nva(s),
        "brown2": brown2(s, b1),
        "vorticity_squared": vorticity_squared(s),
        "ti1": ti1(s),
        "ngm1": ngm1(s),
        "ti2": ti2(s),
        "f2d": f2d(s),
        "ncsu1": ncsu1(s),
    }
    for key, da in out.items():
        if "pressure_level" in da.dims:
            da = da.sel(pressure_level=s["lev"])
        out[key] = da.rename(key)
    return out


def compute_month(ds_raw: xr.Dataset, convention: dict, chunk_days: int = 1,
                  verbose: bool = True) -> xr.Dataset:
    """All 21 diagnostics for one file, `chunk_days` days at a time, as float32.

    Chunking only limits memory. Each chunk is computed with one extra time
    step on either side and then trimmed, so the time derivative in f2d sees
    the same neighbours as in a single pass; the result is identical to
    chunk_days=0 (whole file). Files are NOT stitched to each other: the first
    and last step of every file get a one-sided d/dt.
    """
    n = ds_raw.sizes["time"]
    times = ds_raw["time"].values
    step = chunk_days * C.STEPS_PER_DAY if chunk_days > 0 else n
    parts: dict[str, list] = {}
    for a in range(0, n, step):
        b = min(a + step, n)
        chunk = prepare(ds_raw.isel(time=slice(max(0, a - 1), min(n, b + 1))))
        for key, da in compute_21(chunk, convention).items():
            parts.setdefault(key, []).append(da.reindex(time=times[a:b]).astype(np.float32))
        if verbose:
            print(f"    steps {a:4d}-{b:4d} of {n}   peak memory {peak_memory_gb():.1f} GB", flush=True)
    return xr.Dataset({k: xr.concat(p, dim="time") if len(p) > 1 else p[0]
                       for k, p in parts.items()})


# ===========================================================================
# Driver for one month (called by main.py)
# ===========================================================================
def peak_memory_gb() -> float:
    """Peak resident memory of this process (sacct does not report it on ADA)."""
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def convention_attrs(name: str) -> dict:
    """What gets written into every output store, and compared on re-runs."""
    conv = C.CONVENTIONS[name]
    return {"convention": name, **{f"convention_{k}": v for k, v in conv.items()},
            "f2d_variant_formula": F2D_VARIANTS[conv["f2d_variant"]],
            "target_level_hPa": C.TARGET_LEVEL_HPA, "written_by": "essentials/diagnostics.py"}


def store_is_current(path: Path, convention_name: str) -> bool:
    """True if a COMPLETE store built under this convention already exists."""
    zattrs = Path(path) / ".zattrs"
    if not (Path(path) / ".zmetadata").exists() or not zattrs.exists():
        return False
    have = json.loads(zattrs.read_text())
    want = convention_attrs(convention_name)
    return all(have.get(k) == v for k, v in want.items() if k.startswith("convention"))


def run_month(in_path, out_path, convention_name: str, chunk_days: int = 1,
              skip_if_current: bool = True) -> int:
    """Read one GRIB month, write one zarr store of 21 float32 diagnostics."""
    in_path, out_path = Path(in_path), Path(out_path)
    if skip_if_current and store_is_current(out_path, convention_name):
        print(f"SKIP {out_path.name}: already computed under '{convention_name}'")
        return 0
    if not in_path.exists():
        print(f"!! input not found: {in_path}")
        return 1
    t0 = time.time()
    convention = C.CONVENTIONS[convention_name]
    print(f">>> {in_path.name} -> {out_path}   convention={convention_name} {convention}")
    ds_out = compute_month(open_era5(in_path), convention, chunk_days)
    ds_out = ds_out[[k for k in C.DIAGNOSTIC_KEYS]]
    ds_out.attrs.update(convention_attrs(convention_name))
    ds_out.attrs["source_file"] = in_path.name

    empty = [k for k in ds_out.data_vars if not np.isfinite(ds_out[k].values).any()]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds_out.to_zarr(out_path, mode="w", zarr_format=2, consolidated=True)
    print(f"    wrote 21 diagnostics in {(time.time() - t0) / 60:.1f} min, "
          f"peak memory {peak_memory_gb():.1f} GB")
    if empty:
        print(f"!! no finite values in: {empty}")
        return 2
    return 0
