"""
Self-contained checks of the building blocks (no original repo, no rojak, no ERA5).

    pixi run python -m pytest essentials/tests/test_operators.py -v

Each test states the property it checks. Together with compare_with_original.py
(bit-for-bit against the original pipeline) these are the evidence that the
essentials code is right.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import config as C              # noqa: E402
import diagnostics as D         # noqa: E402
import thresholds as T          # noqa: E402

F = 1 / 298.257223563
E2 = F * (2 - F)
LAT = np.arange(60, 29.99, -0.25)
LON = np.arange(-75, 0.001, 0.25)


def radii(lat_deg):
    s = np.sin(np.deg2rad(lat_deg))
    N = C.A_WGS84 / np.sqrt(1 - E2 * s ** 2)
    M = C.A_WGS84 * (1 - E2) / (1 - E2 * s ** 2) ** 1.5
    return N, M


def as_field(values):
    return xr.DataArray(values, dims=("latitude", "longitude"),
                        coords={"latitude": LAT, "longitude": LON})


def test_scale_factors_match_closed_form():
    """PROJ's k_x, k_y equal a/(N cos phi) and a/M(phi) to 1e-10."""
    geo = D.grid_geometry(LAT, LON, "map")
    N, M = radii(LAT)
    k_x = C.A_WGS84 / (N * np.cos(np.deg2rad(LAT)))
    k_y = C.A_WGS84 / M
    assert np.allclose(geo["k_x"].values, k_x[:, None], rtol=1e-10, atol=0)
    assert np.allclose(geo["k_y"].values, k_y[:, None], rtol=1e-10, atol=0)


def test_dx_is_equatorial_arc():
    """dx = a * dlambda (0.25 deg = 27 829.87 m)."""
    geo = D.grid_geometry(LAT, LON, "map")
    assert np.allclose(geo["dx"], 27829.87269831839, rtol=1e-15)


@pytest.mark.parametrize("dy_method", ["map", "geodesic"])
def test_meridional_derivative_against_exact_ellipsoid(dy_method):
    """For f = phi (radians) the exact northward derivative is 1/M(phi), and a
    centred difference of a linear function has no truncation error.
    'map' (audit, F1) reproduces 1/M to ~1e-10; 'geodesic' (baseline) applies
    the ellipsoid correction twice and is off by up to ~0.4 %."""
    phi = np.deg2rad(LAT)
    f = as_field(np.repeat(phi[:, None], LON.size, axis=1))
    _, M = radii(LAT)
    geo = D.grid_geometry(LAT, LON, dy_method)
    rel = np.abs(D.d_dy(f, geo).values[1:-1, 0] * M[1:-1] - 1)
    if dy_method == "map":
        assert rel.max() < 1e-9
    else:
        assert 3e-3 < rel.max() < 5e-3


def test_zonal_flow_has_curvature_vorticity():
    """A uniform eastward wind has no shear but still has vorticity ~ u tan(phi)/R.
    wind_derivatives sees it; a scalar derivative of u would give exactly zero."""
    U = 20.0
    u = as_field(np.full((LAT.size, LON.size), U))
    v = as_field(np.zeros((LAT.size, LON.size)))
    wd = D.wind_derivatives(u, v, D.grid_geometry(LAT, LON, "map"))
    zeta = (wd["dv_dx"] - wd["du_dy"]).values[1:-1, 5]
    _, M = radii(LAT)
    expected = (U * np.tan(np.deg2rad(LAT)) / M)[1:-1]
    assert np.allclose(zeta, expected, rtol=1e-2)
    assert np.allclose(D.d_dy(u, D.grid_geometry(LAT, LON, "map")).values, 0, atol=1e-18)


def test_shorter_arc_wraps():
    """Turning from 350 deg to 10 deg is 20 deg, not 340 deg."""
    a, b = np.deg2rad(np.float32(10)), np.deg2rad(np.float32(350))
    assert np.isclose(D._shorter_arc(a, b), np.deg2rad(20), atol=1e-6)


def test_weighted_percentile_equal_weights_is_hazen():
    """With equal weights the weighted percentile is numpy's 'hazen' method."""
    x = np.random.default_rng(1).lognormal(size=10_001)
    q = [50, 97, 99.9]
    assert np.allclose(T.weighted_percentile(x, np.ones_like(x), q),
                       np.percentile(x, q, method="hazen"), rtol=1e-12)


def test_tail_method_equals_full_sort():
    """Keeping only the top ~12 % plus W_below reproduces the full weighted percentile."""
    rng = np.random.default_rng(2)
    lat = np.linspace(80, -80, 41)
    w_lat = T.latitude_weights(lat)
    x = rng.lognormal(size=(41, 60, 500)).astype(np.float32)
    x[:, :, :50] = 0.0                                       # a tie mass
    weights = np.broadcast_to(w_lat[:, None, None], x.shape)
    q = np.array(list(C.SEVERITIES.values()))
    full = T.weighted_percentile(x, weights, q)

    cut = float(T.weighted_percentile(x[:, :, ::7], weights[:, :, ::7], [88])[0])
    keep = x >= cut
    rows = np.broadcast_to(np.arange(41, dtype=np.int16)[:, None, None], x.shape)[keep]
    w_below = float(np.dot(w_lat, (~keep).sum(axis=(1, 2))))
    w_total = float(np.dot(w_lat, np.full(41, 60 * 500)))
    tail, *_ = T.thresholds_from_tail(x[keep], rows, w_lat, w_below, w_total, q, "test", 0.005)
    assert np.allclose(tail, full, rtol=1e-7)   # float32 input; interp amplifies last bits


def test_guard_refuses_a_cut_that_is_too_high():
    """If the kept tail starts above p97 the calibration stops instead of guessing."""
    x = np.random.default_rng(3).lognormal(size=100_000).astype(np.float32)
    w_lat = np.ones(1)
    cut = np.quantile(x, 0.98)
    keep = x >= cut
    with pytest.raises(SystemExit):
        T.thresholds_from_tail(x[keep], np.zeros(keep.sum(), np.int16), w_lat,
                               float((~keep).sum()), float(x.size), [97.0], "test", 0.005)


def test_ols_recovers_a_known_line():
    import trends
    years = np.arange(1979, 2021, dtype=float)
    fit = trends.ols(years, 0.01 + 2e-4 * (years - 1979))
    assert np.isclose(fit["slope"], 2e-4) and fit["r2"] > 0.999999
    row = trends.trend_row(years, 0.01 + 2e-4 * (years - 1979))
    assert np.isclose(row["change"], 41 * 2e-4 / 0.01)


def test_season_lengths():
    import trends
    assert trends.season_days(2000, "djf") == 91 and trends.season_days(2001, "djf") == 90
    assert trends.season_days(2000, "annual") == 366 and trends.season_days(1979, "annual") == 365
