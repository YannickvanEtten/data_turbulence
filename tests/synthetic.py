"""
tests/synthetic.py
==================
Manufactured-solution machinery for the analytic verification suite.

WHY THIS EXISTS
---------------
14 of the 21 diagnostics ("tier B") were verified by reading rojak's source
against Sharman (2006) Appendix A. 2 more ("tier C" -- rva_magnitude, f2d)
have no second implementation anywhere to compare against. Neither group has
ever been checked against a number that is known to be right.

A cross-implementation check cannot help here: for tier B the second
implementation *is* rojak, and for tier C there is no second implementation.
The only remaining tool is a manufactured solution -- build an atmosphere whose
fields are analytic functions, work out by hand what each diagnostic must
equal, and compare.

DESIGN
------
Fields are analytic in (lambda, phi, p, t) with lambda/phi in radians. Three
properties make the expected values exact rather than approximate:

  * u, v and theta are LINEAR in pressure, so the vertical derivatives that
    2_diagnostics.py takes with `.differentiate("pressure_level")` over three
    levels are exact, not second-order approximations.
  * The horizontal structure is separable and p-independent, so horizontal and
    vertical derivatives can be reasoned about one at a time.
  * The p-dependent part is written as (p - P_TARGET), so it VANISHES at the
    target level: u and v reduce to their purely horizontal parts there, while
    du/dp and dv/dp stay non-zero. Advecting velocities and vertical shears are
    therefore independently controllable.

THE METRIC
----------
rojak differentiates on the WGS84 ellipsoid (nominal grid spacing from
pyproj's Geod, times PROJ's parallel/meridional scale factors), not on a
sphere. Expected values here use the exact ellipsoidal metric factors instead
-- the radius of the parallel N(phi)cos(phi), and the meridional radius of
curvature M(phi) -- computed directly from the WGS84 defining constants and
independent of both PROJ and pyproj.

Measured agreement between the two routes (see verification/):
    d/dx : 6 significant figures
    d/dy : ~0.2 %  (rojak's meridional_scale is very slightly redundant with
                    its geodesic dy; harmless at this magnitude)

That 0.2 % sets the noise floor for these tests. Tolerances are set at 1-2 %,
which is far above the floor and far below any real formula error -- a wrong
term, a missing factor or a flipped sign moves a diagnostic by tens of percent
at least, usually by a factor of 2 or more.
"""
from __future__ import annotations

import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# WGS84 defining constants (independent of pyproj / PROJ)
# ---------------------------------------------------------------------------
WGS84_A = 6378137.0                # semi-major axis, m
WGS84_F = 1.0 / 298.257223563      # flattening
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)   # first eccentricity squared

# rojak's own constants, for quantities that must match it exactly
G = 9.80665
KAPPA = 0.28571428571428564        # R/cp, exactly as rojak's potential_temperature
P_REF = 1000.0                     # hPa

# ---------------------------------------------------------------------------
# Grid and field parameters. Values are physically plausible for 200 hPa but
# their exact size does not matter -- what matters is that every derivative
# below is analytic.
# ---------------------------------------------------------------------------
P_TARGET = 200
LEVELS = [175, 200, 225]

U0, V0 = 25.0, -6.0        # m/s, constant part of the wind
UH, VH = 12.0, 8.0         # m/s, horizontal structure amplitude
UP, VP = 0.12, -0.09       # m/s/hPa, vertical shear amplitude
K_LON, N_LAT, M_LON = 5.0, 4.0, 3.0     # wavenumbers (per radian)

# Time frequency of the manufactured signal. DELIBERATELY SLOW (10-day period,
# so q -- which varies at 2*OMEGA -- turns over in 5 days).
#
# The first version of this suite used a diurnal OMEGA and F2D "failed" with a
# 28 % error. It was not a code error. A centred difference of a signal at
# angular frequency W sampled at spacing h returns the true derivative times
# sin(Wh)/(Wh); for a diurnal q on 3-hourly data that factor is 0.637, i.e.
# 36 % low. Slowing the signal to 10 days puts the factor at 0.996 and isolates
# the FORMULA, which is what this file is for.
#
# The physical caveat that episode uncovered is real and is recorded in
# STATUS.md: on 3-hourly ERA5, F2D's material derivative systematically
# understates anything evolving on sub-diurnal timescales. That is a property
# of the data, not of the code, and no test here can fix it.
OMEGA = 2.0 * np.pi / (10 * 24 * 3600.0)

THETA0, GAMMA_THETA, TH_H = 460.0, -0.55, 6.0   # K, K/hPa, K
Z200, ZP = 11800.0, -55.0   # m, m/hPa  (altitude falls as pressure rises)

Z0_VORT, D0_DIV, PV0 = 8.0e-5, 3.0e-5, 4.0e-6   # s^-1, s^-1, PV units


def meridional_radius(lat_deg):
    """M(phi): meridional radius of curvature on WGS84, metres."""
    s = np.sin(np.deg2rad(lat_deg))
    return WGS84_A * (1.0 - WGS84_E2) / (1.0 - WGS84_E2 * s * s) ** 1.5


def parallel_radius(lat_deg):
    """N(phi)cos(phi): radius of the parallel of latitude on WGS84, metres."""
    lat = np.deg2rad(lat_deg)
    s = np.sin(lat)
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * s * s)
    return n * np.cos(lat)


def make_grid(dlat=0.25, dlon=0.25, lat0=35.0, lat1=55.0, lon0=-40.0, lon1=-10.0,
              n_times=5, dt_hours=3):
    lat = np.arange(lat0, lat1 + 1e-9, dlat)
    lon = np.arange(lon0, lon1 + 1e-9, dlon)
    # minutes, so the convergence tests can halve dt below one hour
    step_min = int(round(dt_hours * 60))
    time = np.array([np.datetime64("2016-07-01T00:00:00") + np.timedelta64(i * step_min, "m")
                     for i in range(n_times)])
    return lat, lon, time


def _coord_arrays(lat, lon, time, levels):
    """Broadcast (lat, lon, time, level) -> radians / seconds / hPa arrays,
    in the (latitude, longitude, time, pressure_level) order rojak wants."""
    phi = np.deg2rad(lat)[:, None, None, None]
    lam = np.deg2rad(lon)[None, :, None, None]
    t_s = ((time - time[0]) / np.timedelta64(1, "s")).astype(float)[None, None, :, None]
    p = np.asarray(levels, dtype=float)[None, None, None, :]
    return phi, lam, t_s, p


def analytic_fields(lat, lon, time, levels=LEVELS):
    """The manufactured atmosphere, as plain numpy arrays.

    u  = U0 + UH sin(K lam) cos(N phi) + UP (p - 200) sin(K lam) cos(w t)
    v  = V0 + VH cos(M lam) sin(N phi) + VP (p - 200) cos(N phi) sin(w t)
    th = TH0 + GAMMA (p - 200) + TH_H sin(K lam) cos(N phi)
    Phi= g [ Z200 + ZP (p - 200) ]
    vo, d, pv: independent analytic fields (ERA5 supplies these directly, so
    the diagnostics that consume them never re-derive them from u and v).
    """
    phi, lam, t_s, p = _coord_arrays(lat, lon, time, levels)
    dp = p - P_TARGET

    u = U0 + UH * np.sin(K_LON * lam) * np.cos(N_LAT * phi) \
        + UP * dp * np.sin(K_LON * lam) * np.cos(OMEGA * t_s)
    v = V0 + VH * np.cos(M_LON * lam) * np.sin(N_LAT * phi) \
        + VP * dp * np.cos(N_LAT * phi) * np.sin(OMEGA * t_s)

    theta = THETA0 + GAMMA_THETA * dp + TH_H * np.sin(K_LON * lam) * np.cos(N_LAT * phi)
    temperature = theta * (p / P_REF) ** KAPPA      # inverse of rojak's Poisson eq.

    geopotential = G * (Z200 + ZP * dp) * np.ones_like(u)

    vo = Z0_VORT * np.sin(K_LON * lam) * np.cos(N_LAT * phi) * np.ones_like(t_s)
    d = D0_DIV * np.cos(K_LON * lam) * np.sin(N_LAT * phi) * np.ones_like(t_s)
    pv = PV0 * np.sin(K_LON * lam + N_LAT * phi) * np.ones_like(t_s)

    shape = np.broadcast(u, v, theta, vo).shape
    return {k: np.broadcast_to(a, shape).copy() for k, a in {
        "u": u, "v": v, "t": temperature, "z": geopotential,
        "vo": vo, "d": d, "pv": pv, "theta": theta,
    }.items()}


def make_dataset(dlat=0.25, dlon=0.25, n_times=5, levels=LEVELS, **grid_kw):
    """An ERA5-shaped xr.Dataset of the manufactured atmosphere, using the
    GRIB short names that `prepare_for_rojak` expects."""
    lat, lon, time = make_grid(dlat=dlat, dlon=dlon, n_times=n_times, **grid_kw)
    f = analytic_fields(lat, lon, time, levels)
    dims = ("latitude", "longitude", "time", "pressure_level")
    coords = {"latitude": lat, "longitude": lon, "time": time,
              "pressure_level": np.asarray(levels)}
    ds = xr.Dataset(
        {name: (dims, f[name]) for name in ("u", "v", "t", "z", "vo", "d", "pv")},
        coords=coords,
    )
    ds["pressure_level"].attrs["units"] = "hPa"
    return ds


# ---------------------------------------------------------------------------
# Analytic derivatives at the TARGET LEVEL, as xr.DataArrays on
# (latitude, longitude, time). Every one is differentiated by hand -- nothing
# here calls a numerical derivative, and nothing here calls rojak.
# ---------------------------------------------------------------------------
class Expect:
    """Closed-form fields and derivatives at P_TARGET."""

    def __init__(self, lat, lon, time, levels=LEVELS):
        self.lat, self.lon, self.time = lat, lon, time
        phi = np.deg2rad(lat)[:, None, None]
        lam = np.deg2rad(lon)[None, :, None]
        t_s = ((time - time[0]) / np.timedelta64(1, "s")).astype(float)[None, None, :]
        self.phi, self.lam, self.t = np.broadcast_arrays(phi, lam, t_s)
        self.mx = parallel_radius(lat)[:, None, None]     # metres per radian of longitude
        self.my = meridional_radius(lat)[:, None, None]   # metres per radian of latitude
        self.dims = ("latitude", "longitude", "time")
        self.coords = {"latitude": lat, "longitude": lon, "time": time}

    def wrap(self, arr):
        return xr.DataArray(np.broadcast_to(arr, self.phi.shape).copy(),
                            dims=self.dims, coords=self.coords)

    # -- wind at the target level (the (p - 200) terms vanish there) --------
    @property
    def u(self):
        return U0 + UH * np.sin(K_LON * self.lam) * np.cos(N_LAT * self.phi)

    @property
    def v(self):
        return V0 + VH * np.cos(M_LON * self.lam) * np.sin(N_LAT * self.phi)

    # -- horizontal derivatives of the wind at the target level ------------
    @property
    def du_dx(self):
        return UH * K_LON * np.cos(K_LON * self.lam) * np.cos(N_LAT * self.phi) / self.mx

    @property
    def du_dy(self):
        return -UH * N_LAT * np.sin(K_LON * self.lam) * np.sin(N_LAT * self.phi) / self.my

    @property
    def dv_dx(self):
        return -VH * M_LON * np.sin(M_LON * self.lam) * np.sin(N_LAT * self.phi) / self.mx

    @property
    def dv_dy(self):
        return VH * N_LAT * np.cos(M_LON * self.lam) * np.cos(N_LAT * self.phi) / self.my

    # -- vertical (altitude) shear ----------------------------------------
    # dZ/dp = ZP, so du/dZ = (du/dp)/ZP.
    @property
    def du_dz(self):
        return UP * np.sin(K_LON * self.lam) * np.cos(OMEGA * self.t) / ZP

    @property
    def dv_dz(self):
        return VP * np.cos(N_LAT * self.phi) * np.sin(OMEGA * self.t) / ZP

    @property
    def vws(self):
        return np.sqrt(self.du_dz ** 2 + self.dv_dz ** 2)

    # -- vector-component derivatives (curvature-corrected) ----------------
    # Differentiating a VECTOR COMPONENT on a curved earth is not the same as
    # differentiating a scalar: the unit vectors rotate with position, adding
    # metric terms. rojak follows metpy, correcting by the map-factor
    # gradients; on a sphere those reduce to a single tan(phi)/M term:
    #
    #     du/dx = (grad u)_x - v tan(phi)/M        dv/dx = (grad v)_x + u tan(phi)/M
    #     du/dy = (grad u)_y                       dv/dy = (grad v)_y
    #
    # Sanity check on the pair: these reproduce the standard spherical forms
    #     div  = du/dx + dv/dy = (grad u)_x + (grad v)_y - v tan(phi)/M
    #     zeta = dv/dx - du/dy = (grad v)_x - (grad u)_y + u tan(phi)/M
    # which is how we know the correction is the right one and not a fudge.
    #
    # This is not a small term. At 45N with v ~ 10 m/s it is ~1.6e-6 s^-1,
    # around a tenth of the raw gradient -- and because deformation is built
    # from DIFFERENCES of these derivatives, omitting it shifts DEF by tens of
    # percent, which is exactly what an early version of this suite saw.
    @property
    def _tan_corr(self):
        return np.tan(np.deg2rad(self.lat))[:, None, None] / meridional_radius(self.lat)[:, None, None]

    @property
    def du_dx_vec(self):
        return self.du_dx - self.v * self._tan_corr

    @property
    def dv_dx_vec(self):
        return self.dv_dx + self.u * self._tan_corr

    @property
    def du_dy_vec(self):
        return self.du_dy

    @property
    def dv_dy_vec(self):
        return self.dv_dy

    # -- kinematic combinations -------------------------------------------
    @property
    def wind_speed(self):
        return np.sqrt(self.u ** 2 + self.v ** 2)

    @property
    def stretching_deformation(self):
        return self.du_dx_vec - self.dv_dy_vec

    @property
    def shearing_deformation(self):
        return self.dv_dx_vec + self.du_dy_vec

    @property
    def wind_divergence(self):
        """Divergence computed from the wind, as TI2 uses it -- NOT the
        independent ERA5 `d` field (`self.divergence`). rojak's TI2 takes
        du/dx and dv/dy from vector_derivatives, so it never sees ERA5's own
        divergence product."""
        return self.du_dx_vec + self.dv_dy_vec

    @property
    def total_deformation(self):
        return np.sqrt(self.stretching_deformation ** 2 + self.shearing_deformation ** 2)

    # -- ERA5-supplied fields ---------------------------------------------
    @property
    def vorticity(self):
        return Z0_VORT * np.sin(K_LON * self.lam) * np.cos(N_LAT * self.phi)

    @property
    def dvort_dx(self):
        return Z0_VORT * K_LON * np.cos(K_LON * self.lam) * np.cos(N_LAT * self.phi) / self.mx

    @property
    def dvort_dy(self):
        return -Z0_VORT * N_LAT * np.sin(K_LON * self.lam) * np.sin(N_LAT * self.phi) / self.my

    @property
    def divergence(self):
        return D0_DIV * np.cos(K_LON * self.lam) * np.sin(N_LAT * self.phi)

    @property
    def potential_vorticity(self):
        return PV0 * np.sin(K_LON * self.lam + N_LAT * self.phi)

    # -- temperature / potential temperature at the target level -----------
    @property
    def theta(self):
        return THETA0 + TH_H * np.sin(K_LON * self.lam) * np.cos(N_LAT * self.phi)

    @property
    def temperature(self):
        return self.theta * (P_TARGET / P_REF) ** KAPPA

    @property
    def dT_dx(self):
        scale = (P_TARGET / P_REF) ** KAPPA
        return scale * TH_H * K_LON * np.cos(K_LON * self.lam) * np.cos(N_LAT * self.phi) / self.mx

    @property
    def dT_dy(self):
        scale = (P_TARGET / P_REF) ** KAPPA
        return -scale * TH_H * N_LAT * np.sin(K_LON * self.lam) * np.sin(N_LAT * self.phi) / self.my

    @property
    def grad_T(self):
        return np.sqrt(self.dT_dx ** 2 + self.dT_dy ** 2)

    # -- vertical temperature / theta gradients ----------------------------
    # theta(p) = THETA0 + GAMMA (p-200) + TH_H S,  T(p) = theta (p/1000)^kappa
    #   dT/dp = (p/1000)^kappa [ GAMMA + theta kappa / p ]
    # and dZ/dp = ZP, so dT/dZ = (dT/dp) / ZP. Unlike u, v and theta, T is NOT
    # linear in p, so rojak's 3-level `.differentiate` carries a genuine
    # O(dp^2) truncation error here -- expect ~1 %, not ~0.02 %.
    @property
    def dT_dz(self):
        scale = (P_TARGET / P_REF) ** KAPPA
        dT_dp = scale * (GAMMA_THETA + self.theta * KAPPA / P_TARGET)
        return dT_dp / ZP

    @property
    def dtheta_dz(self):
        return GAMMA_THETA / ZP

    # -- directional shear (Endlich) ---------------------------------------
    # alpha = atan2(v, u);  d(alpha)/dz = (u dv/dz - v du/dz) / (u^2 + v^2)
    # Endlich = |V| * |d(alpha)/dz| = |u dv/dz - v du/dz| / |V|
    @property
    def endlich(self):
        return np.abs(self.u * self.dv_dz - self.v * self.du_dz) / self.wind_speed

    # -- tier C: relative vorticity advection ------------------------------
    @property
    def rva(self):
        return np.abs(self.u * self.dvort_dx + self.v * self.dvort_dy)

    # -- tier C: isentropic frontogenesis, Sharman A9 ----------------------
    # q = (du/dtheta)^2 + (dv/dtheta)^2 with dtheta/dp = GAMMA_THETA, so
    #   du/dtheta = UP sin(K lam) cos(w t) / GAMMA
    #   dv/dtheta = VP cos(N phi) sin(w t) / GAMMA
    # q is independent of p, which is what makes D/Dt tractable in closed form.
    @property
    def q(self):
        a = UP * np.sin(K_LON * self.lam) * np.cos(OMEGA * self.t) / GAMMA_THETA
        b = VP * np.cos(N_LAT * self.phi) * np.sin(OMEGA * self.t) / GAMMA_THETA
        return a ** 2 + b ** 2

    @property
    def dq_dt(self):
        ca = (UP / GAMMA_THETA) ** 2 * np.sin(K_LON * self.lam) ** 2
        cb = (VP / GAMMA_THETA) ** 2 * np.cos(N_LAT * self.phi) ** 2
        # d/dt [ca cos^2(wt) + cb sin^2(wt)] = w sin(2wt) (cb - ca)
        return OMEGA * np.sin(2 * OMEGA * self.t) * (cb - ca)

    @property
    def dq_dx(self):
        ca = (UP / GAMMA_THETA) ** 2 * np.cos(OMEGA * self.t) ** 2
        return ca * K_LON * np.sin(2 * K_LON * self.lam) / self.mx

    @property
    def dq_dy(self):
        cb = (VP / GAMMA_THETA) ** 2 * np.sin(OMEGA * self.t) ** 2
        return -cb * N_LAT * np.sin(2 * N_LAT * self.phi) / self.my

    @property
    def f2d_isentropic(self):
        return 0.5 * (self.dq_dt + self.u * self.dq_dx + self.v * self.dq_dy)


def relative_error(got, expected, interior=3):
    """Median and 95th-percentile |got-expected| / rms(expected), computed on
    the grid interior only.

    The edges are excluded because np.gradient falls back to a one-sided
    difference there, which is first-order rather than second-order accurate.
    That is a known and accepted property of the scheme, not a defect, and
    including those points would swamp the signal this suite is looking for.

    Normalising by the RMS of the expected field rather than pointwise avoids
    manufacturing enormous "relative errors" at the zero crossings of a
    sinusoid, where the true value is legitimately ~0.
    """
    g = np.asarray(got, dtype=float)
    e = np.asarray(expected, dtype=float)
    g, e = np.broadcast_arrays(g, e)
    sl = (slice(interior, -interior), slice(interior, -interior), Ellipsis)
    g, e = g[sl], e[sl]
    m = np.isfinite(g) & np.isfinite(e)
    scale = np.sqrt(np.mean(e[m] ** 2))
    if scale == 0:
        scale = 1.0
    err = np.abs(g[m] - e[m]) / scale
    return float(np.median(err)), float(np.percentile(err, 95))
