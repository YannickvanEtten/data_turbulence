"""
tests/test_convergence.py
=========================
Grid-refinement (convergence) tests.

WHY THESE MATTER MORE THAN THE TOLERANCE TESTS
-----------------------------------------------
A tolerance test says "the answer is close". A convergence test says "the answer
is close FOR THE RIGHT REASON". If a diagnostic uses a correct second-order
scheme, halving the grid spacing must cut the error by roughly four. A formula
that happens to land near the truth on one grid -- a plausible-looking wrong
coefficient, a term that is small at this resolution -- will NOT show that
scaling. This is the test that distinguishes the two, and it is the reason a
manufactured solution is worth building at all.

Each test refines until the error reaches the ~0.2 % metric floor, below which
the comparison is measuring the difference between two ellipsoid
implementations rather than anything about the diagnostic. The assertions
therefore check the FIRST refinement step, where truncation error still
dominates, and require a clear improvement rather than a precise factor of 4.
"""
from __future__ import annotations

import numpy as np
import pytest

import synthetic as S

MIN_IMPROVEMENT = 2.5   # second order would give 4; allow room before the floor


def _prepared(diag, **kw):
    return diag.prepare_for_rojak(S.make_dataset(**kw))


def test_rva_converges_in_space(diag):
    """Refine dx; the error in |u.grad zeta| must fall roughly as dx^2."""
    errors = []
    for dx in (1.0, 0.5, 0.25):
        ds = _prepared(diag, dlat=dx, dlon=dx)
        got = diag.rva(ds._dataset).transpose("latitude", "longitude", "time").values
        lat, lon, time = S.make_grid(dlat=dx, dlon=dx)
        median, _ = S.relative_error(got, S.Expect(lat, lon, time).rva)
        errors.append(median)

    improvement = errors[0] / errors[1]
    assert improvement > MIN_IMPROVEMENT, (
        f"RVA error fell only {improvement:.2f}x when dx halved "
        f"({errors[0]:.4%} -> {errors[1]:.4%}); a correct second-order scheme "
        f"should give ~4x. Errors at 1.0/0.5/0.25 deg: "
        f"{[f'{e:.4%}' for e in errors]}"
    )
    assert errors[-1] < 0.005, f"finest-grid error {errors[-1]:.4%} is above the metric floor"


def test_f2d_converges_in_time(diag):
    """Refine dt; the error in the material derivative must fall as dt^2.

    This is the test that would have caught the aliasing episode described in
    synthetic.py: a signal too fast for the sampling produces a large error
    that does NOT improve on refinement in space, only in time."""
    errors = []
    for dt in (6, 3, 1.5):
        ds = _prepared(diag, n_times=9, dt_hours=dt)
        got = diag.frontogenesis_isentropic(ds._dataset)
        got = got.transpose("latitude", "longitude", "time").isel(time=slice(1, -1)).values
        lat, lon, time = S.make_grid(n_times=9, dt_hours=dt)
        median, _ = S.relative_error(got, S.Expect(lat, lon, time).f2d_isentropic[:, :, 1:-1])
        errors.append(median)

    improvement = errors[0] / errors[1]
    assert improvement > MIN_IMPROVEMENT, (
        f"F2D error fell only {improvement:.2f}x when dt halved "
        f"({errors[0]:.4%} -> {errors[1]:.4%}). Errors at 6/3/1.5 h: "
        f"{[f'{e:.4%}' for e in errors]}"
    )


def test_ngm2_vertical_error_is_truncation_not_formula(diag):
    """NGM2 sits at ~2.8 % on the project's real 175/200/225 hPa levels.

    This test demonstrates that the 2.8 % is the O(dp^2) truncation error of a
    three-level derivative of a temperature field that is not linear in
    pressure -- bringing the levels closer together drives it to zero. If it
    were a formula error it would stay put, and this test would fail.

    The practical consequence is recorded in STATUS.md: NGM2 as computed on the
    real level spacing carries a few percent of vertical-discretisation error,
    which the percentile calibration largely absorbs (a near-constant factor
    moves the threshold by the same factor) but which is worth knowing about.
    """
    errors = []
    for dp in (25, 12.5, 6.25):
        ds = _prepared(diag, levels=[200 - dp, 200, 200 + dp])
        out, failures = diag.compute_rojak_diagnostics(ds)
        assert not failures
        got = out["ngm2"].sel(pressure_level=200, method="nearest")
        got = got.transpose("latitude", "longitude", "time").values
        lat, lon, time = S.make_grid()
        e = S.Expect(lat, lon, time)
        median, _ = S.relative_error(got, np.abs(e.dT_dz) * e.total_deformation)
        errors.append(median)

    improvement = errors[0] / errors[1]
    assert improvement > MIN_IMPROVEMENT, (
        f"NGM2 error fell only {improvement:.2f}x when the level spacing halved "
        f"({errors[0]:.4%} -> {errors[1]:.4%}). That would mean the residual is "
        f"NOT vertical truncation error, i.e. it is a formula problem. "
        f"Errors at dp=25/12.5/6.25 hPa: {[f'{e:.4%}' for e in errors]}"
    )
    assert errors[0] < 0.05, f"NGM2 at the real 25 hPa spacing is {errors[0]:.4%}, worse than expected"


def test_gradient_metric_matches_exact_ellipsoid(diag):
    """The foundation under every gradient-based diagnostic.

    rojak reaches d/dx by a nominal equatorial grid spacing times PROJ's
    parallel_scale; this compares that against the exact WGS84 metric computed
    from the defining constants, via a field whose derivative is known in
    closed form (longitude in radians has d/dx = 1 / (N cos phi) exactly).

    If this test fails, nothing else in the suite means anything.
    """
    import xarray as xr
    from rojak.core.derivatives import GradientMode, spatial_gradient

    lat = np.arange(35.0, 55.01, 0.25)
    lon = np.arange(-40.0, -9.99, 0.25)
    lat_grid, lon_grid = np.meshgrid(lat, lon, indexing="ij")

    f = xr.DataArray(np.deg2rad(lon_grid), dims=("latitude", "longitude"),
                     coords={"latitude": lat, "longitude": lon})
    got_x = spatial_gradient(f, "deg", GradientMode.GEOSPATIAL)["dfdx"].values
    exact_x = 1.0 / S.parallel_radius(lat_grid)
    rel_x = np.abs(got_x - exact_x) / exact_x
    assert rel_x.max() < 1e-4, f"zonal metric off by {rel_x.max():.2e} (expected < 1e-4)"

    g = xr.DataArray(np.deg2rad(lat_grid), dims=("latitude", "longitude"),
                     coords={"latitude": lat, "longitude": lon})
    grad_y = spatial_gradient(g, "deg", GradientMode.GEOSPATIAL)
    got_y = grad_y["dfdy"].values
    exact_y = 1.0 / S.meridional_radius(lat_grid)
    rel_y = np.abs(got_y - exact_y) / exact_y
    assert rel_y.max() < 5e-3, f"meridional metric off by {rel_y.max():.2e} (expected < 5e-3)"

    # Cross terms must vanish: longitude has no meridional gradient and vice versa.
    assert np.abs(grad_y["dfdx"].values).max() < 1e-15
