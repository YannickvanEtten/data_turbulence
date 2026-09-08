"""
tests/test_metric_terms.py
==========================
The pass/fail criterion for audit fix F1, written before the fix was trusted.

WHY THIS FILE EXISTS
--------------------
`test_geometry_theorems.py` settles WHICH operator is physical (scalar gradient
versus vector derivative) by theorem. It does not check whether the operator's
METRIC TERMS are right -- whether d/dx really carries 1/(a cos phi) and d/dy
really carries 1/M. That question cannot be settled by comparing two
implementations, because rojak is the only implementation; it has to be settled
against a field whose gradient is known in closed form on the WGS84 ellipsoid.

Two such fields are used here, and neither involves a derivative operator:

  * f = the eastward arc length along a parallel, N(phi) cos(phi) * dlambda.
    By construction df/dx == 1 everywhere.
  * f = the northward meridional arc length, integral of M(phi) dphi.
    By construction df/dy == 1 everywhere.

WHAT THEY FOUND (2026-09-08, AUDIT_diagnostics_vs_literature.md 3.1)
--------------------------------------------------------------------
  * d/dx is EXACT. The 1/cos(phi) metric term is present and correct on the
    ellipsoid. `nominal_grid_spacing` takes dx at the equator, so it is the
    plate-carree MAP distance a*dlambda, and parallel_scale = a/(N cos phi)
    converts it in one step.
  * d/dy is DOUBLE-CORRECTED, by exactly meridional_scale = a/M(phi).
    `nominal_grid_spacing` takes dy as the true geodesic meridional arc
    M*dphi -- already d/dy -- and `spatial_gradient` then multiplies by a/M
    anyway. Measured: +0.4213 % at 30N falling to -0.2678 % at 75N. Small,
    monotone in latitude, and therefore invisible in a global median and
    fatal-in-principle in a regional exceedance frequency.

The two errors are different because the two axes are not in the same
coordinate: dx is a map distance and dy is a ground distance. That asymmetry is
the bug, and `meridional_metric_patch` removes it by making dy a map distance
too.

These tests are the reason F1 can be merged on evidence rather than on argument.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from rojak.core.derivatives import (
    GradientMode,
    LatLonUnits,
    VelocityDerivative,
    spatial_gradient,
    vector_derivatives,
)

# WGS84, the ellipsoid pyproj's scale factors are defined against.
A = 6378137.0
F = 1.0 / 298.257223563
E2 = F * (2.0 - F)

LATS_CHECKED = (30.0, 45.0, 60.0, 75.0)


def _prime_vertical_radius(phi):        # N(phi)
    return A / np.sqrt(1.0 - E2 * np.sin(phi) ** 2)


def _meridional_radius(phi):            # M(phi)
    return A * (1.0 - E2) / (1.0 - E2 * np.sin(phi) ** 2) ** 1.5


def _grid(descending: bool):
    lats = np.arange(80.0, 19.999, -0.25)
    if not descending:
        lats = lats[::-1].copy()
    lons = np.arange(-60.0, 0.001, 0.25)
    return lats, lons


def _wrap(field, lats, lons):
    return xr.DataArray(field, dims=("latitude", "longitude"),
                        coords={"latitude": lats, "longitude": lons})


def _zonal_arc_field(lats, lons):
    """f = eastward arc along the parallel. Exactly linear in longitude, so the
    finite difference is exact and any error is the metric term's."""
    phi, lam = np.deg2rad(lats), np.deg2rad(lons)
    scale = _prime_vertical_radius(phi) * np.cos(phi)
    return _wrap(scale[:, None] * (lam - lam[0])[None, :], lats, lons)


def _meridional_arc_field(lats, lons):
    """f = northward meridional arc measured from the southernmost row.

    The arc comes from `pyproj.Geod.inv` along a meridian, which is the EXACT
    ellipsoidal geodesic -- not a quadrature of M(phi), which would leave an
    O(dphi^2) residual of order 5e-8 and blunt the test. Geod is independent of
    every derivative operator under test here, so this stays an external
    reference rather than a circular one."""
    from pyproj import Geod

    geod = Geod(ellps="WGS84")
    base = float(np.min(lats))
    zeros = np.zeros_like(lats)
    _, _, arc = geod.inv(zeros, np.full_like(lats, base), zeros, lats)
    arc = np.asarray(arc, dtype=float)
    return _wrap(np.tile(arc[:, None], (1, len(lons))), lats, lons)


def _at(values, lats, lat_deg, col=120):
    return float(values[int(np.argmin(np.abs(lats - lat_deg))), col])


# ---------------------------------------------------------------------------
# The zonal derivative is already exact, with and without the fix
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("descending", [True, False])
def test_zonal_derivative_is_exact(descending):
    lats, lons = _grid(descending)
    dfdx = spatial_gradient(_zonal_arc_field(lats, lons),
                            LatLonUnits.DEG, GradientMode.GEOSPATIAL)["dfdx"].values
    for lat in LATS_CHECKED:
        assert _at(dfdx, lats, lat) == pytest.approx(1.0, abs=1e-9), (
            f"d/dx at {lat}N should be exactly 1; the 1/(a cos phi) metric term "
            f"is correct in the shipped code and must stay correct")


@pytest.mark.parametrize("descending", [True, False])
def test_the_fix_does_not_disturb_the_zonal_derivative(diag, descending):
    """F1 touches dy only. If it moves dx, it is the wrong patch."""
    lats, lons = _grid(descending)
    field = _zonal_arc_field(lats, lons)
    with diag.meridional_metric_patch(enabled=True):
        dfdx = spatial_gradient(field, LatLonUnits.DEG,
                                GradientMode.GEOSPATIAL)["dfdx"].values
    for lat in LATS_CHECKED:
        assert _at(dfdx, lats, lat) == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# The defect, stated as a measurement rather than a claim
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("descending", [True, False])
def test_meridional_derivative_is_double_corrected_without_the_fix(descending):
    """Baseline d/dy equals a/M(phi) -- exactly the meridional map-scale factor,
    applied a second time to a delta that was already a ground distance.

    This test PASSES on the broken behaviour. It exists so that the defect is
    pinned as a number, and so that a silent upstream fix (a new rojak pin)
    announces itself here instead of quietly changing 15 diagnostics.
    """
    lats, lons = _grid(descending)
    dfdy = spatial_gradient(_meridional_arc_field(lats, lons),
                            LatLonUnits.DEG, GradientMode.GEOSPATIAL)["dfdy"].values
    for lat in LATS_CHECKED:
        expected = A / _meridional_radius(np.deg2rad(lat))
        assert _at(dfdy, lats, lat) == pytest.approx(expected, rel=2e-5), (
            f"baseline d/dy at {lat}N should be inflated by a/M = {expected:.6f}")


def test_the_defect_is_latitude_structured_not_a_constant():
    """The reason this matters at all: it is monotone in latitude, so it cannot
    be absorbed by a percentile calibration the way a constant can."""
    lats, lons = _grid(descending=True)
    dfdy = spatial_gradient(_meridional_arc_field(lats, lons),
                            LatLonUnits.DEG, GradientMode.GEOSPATIAL)["dfdy"].values
    errs = [_at(dfdy, lats, lat) - 1.0 for lat in LATS_CHECKED]
    assert all(a > b for a, b in zip(errs, errs[1:])), (
        f"error should decrease monotonically with latitude, got {errs}")
    assert errs[0] > 0 > errs[-1], f"should change sign across the band, got {errs}"
    assert (errs[0] - errs[-1]) == pytest.approx(0.0069, abs=5e-4), (
        f"peak-to-peak spread over 30-75N should be ~0.69 %, got {errs}")


# ---------------------------------------------------------------------------
# The fix
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("descending", [True, False])
def test_meridional_metric_patch_makes_the_derivative_exact(diag, descending):
    lats, lons = _grid(descending)
    field = _meridional_arc_field(lats, lons)
    with diag.meridional_metric_patch(enabled=True) as active:
        assert active is True
        dfdy = spatial_gradient(field, LatLonUnits.DEG,
                                GradientMode.GEOSPATIAL)["dfdy"].values
    for lat in LATS_CHECKED:
        # 2e-7, not 1e-9. The residual is the SECOND-ORDER TRUNCATION ERROR of
        # a centred difference applied to an arc that is nonlinear in phi
        # (~dphi^2/6 * d3(arc)/dphi3, measured at 3e-8 on a 0.25 deg grid), not
        # a metric error -- the zonal test above is exact to 1e-9 precisely
        # because its field is exactly linear in lambda. 2e-7 is still four
        # orders of magnitude below the 4.2e-3 defect being fixed, so the test
        # discriminates by a wide margin.
        assert _at(dfdy, lats, lat) == pytest.approx(1.0, abs=2e-7), (
            f"patched d/dy at {lat}N should be 1 to within finite-difference "
            f"truncation")


def test_patch_is_a_no_op_when_disabled(diag):
    lats, lons = _grid(descending=True)
    field = _meridional_arc_field(lats, lons)
    plain = spatial_gradient(field, LatLonUnits.DEG,
                             GradientMode.GEOSPATIAL)["dfdy"].values
    with diag.meridional_metric_patch(enabled=False) as active:
        assert active is False
        gated = spatial_gradient(field, LatLonUnits.DEG,
                                 GradientMode.GEOSPATIAL)["dfdy"].values
    assert np.array_equal(plain, gated)


def test_patch_is_reverted_on_exit(diag):
    """In-memory and reversible. A leaked patch would silently change every
    later diagnostic in the same process -- including 4_verify.py's oracle."""
    import rojak.core.derivatives as der
    before = der.nominal_grid_spacing
    with diag.meridional_metric_patch(enabled=True):
        assert der.nominal_grid_spacing is not before
    assert der.nominal_grid_spacing is before

    with pytest.raises(RuntimeError):
        with diag.meridional_metric_patch(enabled=True):
            raise RuntimeError("boom")
    assert der.nominal_grid_spacing is before, "patch leaked through an exception"


# ---------------------------------------------------------------------------
# The patch reaches vector_derivatives too, and does not break it
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("enabled", [False, True])
def test_solid_body_zonal_flow_stays_non_divergent(diag, enabled):
    """u = U cos(phi), v = 0 has zero divergence on a sphere. `spatial_gradient`
    is called from inside `vector_derivatives`, so patching
    `nominal_grid_spacing` reaches deformation, the Jacobian and the spherical
    Laplacian as well -- this checks the patch does not damage that path."""
    lats, lons = _grid(descending=True)
    phi = np.deg2rad(lats)
    u = _wrap(np.tile((30.0 * np.cos(phi))[:, None], (1, len(lons))), lats, lons)
    v = _wrap(np.zeros((len(lats), len(lons))), lats, lons)
    with diag.meridional_metric_patch(enabled=enabled):
        vd = vector_derivatives(u, v, LatLonUnits.DEG)
    div = (vd[VelocityDerivative.DU_DX] + vd[VelocityDerivative.DV_DY]).values
    assert np.nanmax(np.abs(div)) < 1e-15, (
        f"solid-body zonal flow must stay non-divergent, got "
        f"{np.nanmax(np.abs(div)):.3e}")
