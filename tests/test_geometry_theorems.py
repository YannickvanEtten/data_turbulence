"""
tests/test_geometry_theorems.py
===============================
Which derivative operator is the physical one, settled by theorem.

WHY THIS FILE EXISTS
--------------------
`_grad()`'s docstring in 2_diagnostics.py has always said that scalar gradients
are wrong for differentiating a wind component. A docstring is not evidence, and
"the comment says so" is not a reason to accept a term or to reject one. These
tests replace that assertion with two theorems that hold independently of
anyone's convention:

  * STOKES: the circulation of a wind field around a closed loop equals the
    area integral of its vorticity. Compute the circulation by a pure line
    integral -- no derivative operator anywhere -- and it decides which
    candidate vorticity is real.

  * DIVERGENCE THEOREM: the flux of grad(f) out through a closed boundary
    equals the area integral of the Laplacian of f. Same idea: the boundary
    integral involves no second derivative, so it arbitrates.

WHAT THEY FOUND (2026-08-26)
----------------------------
Both are far more than refinements at mid-latitudes:

  * vorticity from scalar gradients is ~100% wrong -- WRONG SIGN. A zonal flow
    on a sphere carries curvature vorticity u*tan(phi)/M even with no shear at
    all, and the scalar operator cannot see it.
  * the flat Laplacian Phi_xx + Phi_yy is ~99% wrong. The missing
    -tan(phi)/M * dPhi/dy term is the bulk of the answer.

Both errors were present in ubf() (diagnostic #13) and in rojak, which is why
4_verify.py's cross-check could never have caught them: the two implementations
shared the mistake. That is the general lesson -- agreement between two
implementations is only as good as their common assumptions, and a theorem has
no common assumptions with either.

The residual ~0.3% in the "correct" rows is the sphere-vs-WGS84 difference in
the closed-form tan(phi)/M used for the expectation; the shipped code obtains
the same correction through rojak's ellipsoidal scale factors instead.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import synthetic as S


def _metric(lat_deg):
    return S.parallel_radius(lat_deg), S.meridional_radius(lat_deg)


def test_stokes_says_vorticity_needs_the_curvature_term():
    """Circulation / area is the physical vorticity. Compare both operators."""
    u_of = lambda lam, phi: S.U0 + S.UH * np.sin(S.K_LON * lam) * np.cos(S.N_LAT * phi)
    v_of = lambda lam, phi: S.V0 + S.VH * np.cos(S.M_LON * lam) * np.sin(S.N_LAT * phi)

    p1, p2 = np.deg2rad([42.0, 48.0])
    l1, l2 = np.deg2rad([-30.0, -22.0])
    lam = np.linspace(l1, l2, 20001)
    phi = np.linspace(p1, p2, 20001)
    mx = lambda p: S.parallel_radius(np.rad2deg(p))
    my = lambda p: S.meridional_radius(np.rad2deg(p))

    circulation = (
        np.trapezoid(u_of(lam, p1) * mx(p1), lam)
        - np.trapezoid(u_of(lam, p2) * mx(p2), lam)
        + np.trapezoid(v_of(l2, phi) * my(phi), phi)
        - np.trapezoid(v_of(l1, phi) * my(phi), phi)
    )
    area = np.trapezoid(mx(phi) * my(phi), phi) * (l2 - l1)
    zeta_physical = circulation / area

    LAM, PHI = np.meshgrid(np.linspace(l1, l2, 1501), np.linspace(p1, p2, 1501))
    MX, MY = mx(PHI), my(PHI)
    dv_dlam = -S.VH * S.M_LON * np.sin(S.M_LON * LAM) * np.sin(S.N_LAT * PHI)
    du_dphi = -S.UH * S.N_LAT * np.sin(S.K_LON * LAM) * np.sin(S.N_LAT * PHI)

    zeta_scalar = dv_dlam / MX - du_dphi / MY
    zeta_vector = zeta_scalar + u_of(LAM, PHI) * np.tan(PHI) / MY

    def box_mean(f):
        return np.trapezoid(np.trapezoid(f * MX * MY, LAM[0], axis=1), PHI[:, 0]) / area

    err_scalar = abs(box_mean(zeta_scalar) - zeta_physical) / abs(zeta_physical)
    err_vector = abs(box_mean(zeta_vector) - zeta_physical) / abs(zeta_physical)

    assert err_vector < 0.01, (
        f"curvature-corrected vorticity is {err_vector:.2%} from the Stokes value; "
        f"expected under 1%")
    assert err_scalar > 0.5, (
        f"scalar-gradient vorticity is only {err_scalar:.2%} from the Stokes value. "
        f"It was ~100% wrong when this was established; if that is no longer true, "
        f"the geometry assumption behind vector_derivatives() has changed and every "
        f"diagnostic that differentiates u or v needs re-examining.")


def test_divergence_theorem_says_laplacian_needs_the_curvature_term():
    """Boundary flux / area is the physical Laplacian.

    Uses a latitude-only field, so the east/west edges carry zero flux and the
    correction term cannot average itself away -- an earlier version of this
    check used a field whose sin(N*phi) happened to vanish at the test latitude
    and reported, wrongly, that the two candidates agreed.
    """
    PH, N = 900.0, 2.0
    my = lambda p: S.meridional_radius(np.rad2deg(p))
    mx = lambda p: S.parallel_radius(np.rad2deg(p))
    Phi_y = lambda p: -PH * N * np.sin(N * p) / my(p)

    p0 = np.deg2rad(45.0)
    assert abs(Phi_y(p0)) > 1e-6, "test field is degenerate at the sample latitude"

    h = np.deg2rad(0.03125)
    p1, p2 = p0 - h, p0 + h
    phis = np.linspace(p1, p2, 200001)
    physical = (Phi_y(p2) * mx(p2) - Phi_y(p1) * mx(p1)) / np.trapezoid(mx(phis) * my(phis), phis)

    eps = 1e-7
    phi_yy = (Phi_y(p0 + eps) - Phi_y(p0 - eps)) / (2 * eps * my(p0))
    flat = phi_yy
    spherical = phi_yy - np.tan(p0) / my(p0) * Phi_y(p0)

    err_flat = abs(flat - physical) / abs(physical)
    err_sph = abs(spherical - physical) / abs(physical)

    assert err_sph < 0.01, f"spherical Laplacian is {err_sph:.2%} off; expected under 1%"
    assert err_flat > 0.5, (
        f"flat Laplacian is only {err_flat:.2%} off. It was ~99% wrong when this was "
        f"established; if that has changed, revisit ubf()'s Laplacian.")


def test_ubf_uses_the_spherical_laplacian_and_vector_jacobian(diag, prepared):
    """Regression guard on ubf() itself.

    Recomputes UBF with the OLD flat geometry and asserts the shipped function
    does not reproduce it. Cheap insurance: this is a correction that is easy to
    lose in a refactor, produces no error when lost, and shifts 41% of cells by
    more than 10% on real ERA5.
    """
    ds = prepared._dataset
    lvl = 200
    u = diag._sel_level(ds, "eastward_wind", lvl)
    v = diag._sel_level(ds, "northward_wind", lvl)
    phi_g = diag._sel_level(ds, "geopotential", lvl)
    zeta = diag._sel_level(ds, "vorticity", lvl)
    lat_rad = np.deg2rad(ds["latitude"])
    f = (2 * diag.OMEGA * np.sin(lat_rad)).broadcast_like(u.isel(time=0))
    beta = (2 * diag.OMEGA * np.cos(lat_rad) / diag.R_EARTH).broadcast_like(u.isel(time=0))

    dpx, dpy = diag._grad(phi_g)
    flat_lap = diag._grad(dpx)[0] + diag._grad(dpy)[1]
    du_dx, du_dy = diag._grad(u)
    dv_dx, dv_dy = diag._grad(v)
    scalar_jac = du_dx * dv_dy - du_dy * dv_dx
    old = np.abs(flat_lap - 2 * scalar_jac - f * zeta + beta * u)

    new = diag.ubf(ds)
    diff = np.abs(np.asarray(new) - np.asarray(old))[3:-3, 3:-3]
    scale = np.nanmedian(np.abs(np.asarray(new)[3:-3, 3:-3]))
    assert np.nanmedian(diff) / scale > 1e-4, (
        "ubf() now matches the old flat-geometry formula. The spherical "
        "Laplacian correction and/or the vector-derivative Jacobian has been "
        "reverted -- see STATUS.md section 4g."
    )
