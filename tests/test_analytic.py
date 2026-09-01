"""
tests/test_analytic.py
======================
Analytic (manufactured-solution) verification of all 21 CAT diagnostics.

WHAT THIS SUITE ESTABLISHES, AND WHAT IT DOES NOT
-------------------------------------------------
It establishes that each diagnostic CORRECTLY COMPUTES ITS STATED FORMULA. Every
expected value below is derived by hand from the formula in rojak's own
docstring (each of which cites its source paper: Ellrod & Knapp 1992 for TI1 and
TI2, Reap 1996 / Sharman 2016 for NGM1 and NGM2, Brown 1973 for Brown1,
Sharman 2006 for NVA and NCSU1, Endlich 1964 for Endlich) or, for the seven
hand-written diagnostics, from the formula in 2_diagnostics.py's own docstring.

It does NOT establish that those formulas are the ones Prosser (2023) intended.
That is a question about the literature, not about the code, and it needs the
PDFs in the Articles folder rather than a test runner. The two questions are
genuinely separate and conflating them is how a project convinces itself it has
verified more than it has.

WHY MANUFACTURED SOLUTIONS
--------------------------
The 14 rojak diagnostics could not be checked by cross-implementation, because
the second implementation would be rojak itself. rva_magnitude and f2d could not
be checked that way either, because no second implementation exists anywhere.
Building an atmosphere with known analytic derivatives is the only remaining
option, and it is the stronger one: it compares against truth, not against
another guess.

TOLERANCES
----------
Two independent error sources set the floor:

  * ~0.2 % from the metric. rojak differentiates on the WGS84 ellipsoid via
    PROJ scale factors; the expectations here use the exact ellipsoidal metric
    computed from the WGS84 defining constants. The two agree to 6 significant
    figures zonally and ~0.2 % meridionally.
  * O(h^2) truncation from the centred differences, which is what
    `test_convergence.py` measures directly.

Real formula errors are nothing like this small. A wrong term, a missing factor,
a flipped sign or an omitted curvature correction moves a diagnostic by tens of
percent -- as this suite demonstrated twice while it was being written (see the
notes on DEF and on TI2 below). A 2 % tolerance separates the two cleanly.
"""
from __future__ import annotations

import numpy as np
import pytest

import synthetic as S

OMEGA_EARTH = 7.292115e-05

# Tolerances on the 95th-percentile relative error over the grid interior.
TOL = 0.02          # 2 % -- the normal case
TOL_LOOSE = 0.05    # 5 % -- diagnostics with a genuine O(dp^2) vertical term


def at_target(da, level=200):
    if "pressure_level" in da.dims:
        da = da.sel(pressure_level=level, method="nearest")
    return da.transpose("latitude", "longitude", "time").values


def check(got, expected, tol=TOL, label=""):
    got = np.asarray(got, dtype=float)
    expected = np.broadcast_to(np.asarray(expected, dtype=float), got.shape)
    median, p95 = S.relative_error(got, expected)
    assert p95 < tol, (
        f"{label}: p95 relative error {p95:.3%} exceeds {tol:.1%} "
        f"(median {median:.3%}). Expected and computed fields disagree by far "
        f"more than the ~0.2% metric floor, which points at the FORMULA, not "
        f"at discretisation."
    )
    return median, p95


# ===========================================================================
# Tier C -- the two diagnostics that had no verification of any kind
# ===========================================================================
class TestTierC:
    """rva_magnitude and f2d. No second implementation exists for either, so
    before this suite the only evidence for them was that their units looked
    plausible -- and for f2d even that was computed over an empty array."""

    def test_rva_magnitude(self, diag, prepared, expect):
        """RVA = |u dzeta/dx + v dzeta/dy|, advection of ERA5's own vorticity
        field by the wind."""
        got = at_target(diag.rva(prepared._dataset))
        check(got, expect.rva, label="rva_magnitude")

    def test_f2d_isentropic(self, diag, prepared, expect):
        """Isentropic frontogenesis, Sharman A9:
            0.5 D/Dt [ (du/dtheta)^2 + (dv/dtheta)^2 ]
        Interior timesteps only -- d/dt is one-sided at the first and last
        step, which is first-order and would dominate the comparison.

        VARIANT A EXPLICITLY, not the default. synthetic.py's manufactured
        solution is the SIGNED closed form, and since 2026-08-30 the default is
        variant C = |A| (FORMULA_AUDIT.md 10.4). Calling the default here would
        compare a magnitude against a signed expectation and report ~200 %
        error on every cell where the true value is negative -- which is
        exactly what it did, in the first run after the default changed.

        Pinning variant A keeps this test doing the job it was built for:
        verifying the material-derivative arithmetic against a closed form.
        The production quantity is covered by the next test.
        """
        got = diag.frontogenesis_isentropic(prepared._dataset, variant="A")
        got = got.transpose("latitude", "longitude", "time").isel(time=slice(1, -1)).values
        check(got, expect.f2d_isentropic[:, :, 1:-1], label="f2d (variant A)")

    def test_f2d_default_is_the_magnitude(self, diag, prepared, expect):
        """#20 as it now reaches disk, against the same closed form under an
        absolute value.

        Splitting the two matters for diagnosis, not just for bookkeeping: if
        this passes and the one above fails, the material derivative is wrong;
        if the one above passes and this fails, the variant wiring is wrong.
        A single combined assertion could not tell those apart."""
        got = diag.frontogenesis_isentropic(prepared._dataset)
        got = got.transpose("latitude", "longitude", "time").isel(time=slice(1, -1)).values
        check(got, np.abs(expect.f2d_isentropic[:, :, 1:-1]), label="f2d (default)")

    def test_f2d_material_derivative_parts(self, diag, prepared, expect):
        """Decompose A9 so a failure localises to d/dt, d/dx or d/dy rather
        than just reporting that the total is wrong."""
        ds = prepared._dataset
        theta = diag._potential_temperature(ds["temperature"], ds["temperature"]["pressure_level"])
        du = diag.theta_derivative_on_pressure_level(ds["eastward_wind"], theta)
        dv = diag.theta_derivative_on_pressure_level(ds["northward_wind"], theta)
        q = (du ** 2 + dv ** 2).sel(pressure_level=200).transpose("latitude", "longitude", "time")

        check(q.values, expect.q, label="q = (du/dtheta)^2 + (dv/dtheta)^2")
        check(q.differentiate("time", datetime_unit="s").isel(time=slice(1, -1)).values,
              expect.dq_dt[:, :, 1:-1], label="dq/dt")
        gx, gy = diag._grad(q)
        check(gx.transpose("latitude", "longitude", "time").values, expect.dq_dx, label="dq/dx")
        check(gy.transpose("latitude", "longitude", "time").values, expect.dq_dy, label="dq/dy")


# ===========================================================================
# Tier B -- the 14 taken from rojak, previously verified by reading only
# ===========================================================================
class TestPassthrough:
    """Diagnostics that consume an ERA5 product directly rather than deriving
    it. Cheap to check and worth checking: a silent |.| or a missing square
    here would be invisible in any plot."""

    def test_wind_speed(self, rojak_out, expect):
        check(at_target(rojak_out["wind_speed"]), expect.wind_speed, label="wind_speed")

    def test_vorticity_squared(self, rojak_out, expect):
        check(at_target(rojak_out["vorticity_squared"]), expect.vorticity ** 2,
              label="vorticity_squared")

    def test_horizontal_divergence_is_absolute(self, rojak_out, expect):
        """Note the absolute value: the diagnostic is |delta|, so it cannot
        distinguish convergence from divergence. That is intended (W&J list it
        as a magnitude) but it is worth having pinned down by a test."""
        check(at_target(rojak_out["horizontal_divergence"]), np.abs(expect.divergence),
              label="horizontal_divergence")

    def test_magnitude_pv_is_absolute(self, rojak_out, expect):
        check(at_target(rojak_out["magnitude_pv"]), np.abs(expect.potential_vorticity),
              label="magnitude_pv")


class TestGradientBased:
    def test_temperature_gradient(self, rojak_out, expect):
        check(at_target(rojak_out["temperature_gradient"]), expect.grad_T,
              label="temperature_gradient")

    def test_vertical_wind_shear(self, rojak_out, expect):
        """Sv = sqrt((du/dZ)^2 + (dv/dZ)^2), altitude from geopotential."""
        check(at_target(rojak_out["vertical_wind_shear"]), expect.vws,
              label="vertical_wind_shear")

    def test_endlich(self, rojak_out, expect):
        """|v| |dpsi/dz| with psi the wind direction, which reduces to
        |u dv/dz - v du/dz| / |v|."""
        check(at_target(rojak_out["endlich"]), expect.endlich, label="endlich")


class TestDeformationFamily:
    """DEF and the four diagnostics built on it.

    These five all failed on the first run of this suite, by 20-40 %. The cause
    was in the EXPECTATION, not in rojak: differentiating a vector component on
    a curved earth adds metric terms that differentiating a scalar does not, and
    the first version of synthetic.py omitted them. rojak applies them (metpy's
    map-factor correction, which on a sphere reduces to tan(phi)/M); once the
    expectation did too, all five dropped to well under 0.2 %.

    Worth recording because it cuts both ways: it is real evidence that rojak
    handles spherical geometry correctly, and it is a warning that any future
    hand-written diagnostic differentiating u or v MUST use
    vector_derivatives() rather than the scalar _grad()."""

    def test_deformation_is_squared(self, rojak_out, expect):
        """rojak's DEF returns total deformation SQUARED.

        This asserts rojak's RAW output and stays that way deliberately --
        `rojak_out` is compute_rojak_diagnostics, which is a faithful
        passthrough. As of 2026-08-29 the square root is taken one level up,
        in compute_all_21, so that everything reaching disk holds DEF rather
        than DEF^2 (FORMULA_AUDIT.md §5). Before that it was taken in
        3_pipeline.py, on the comparison-table path only.

        The pairing this test used to guard now lives in
        tests/test_audit_fixes.py::TestDeformationUnsquared, which checks both
        halves at once."""
        check(at_target(rojak_out["deformation"]), expect.total_deformation ** 2,
              label="deformation (squared)")

    def test_ti1(self, rojak_out, expect):
        """Ellrod & Knapp (1992): TI1 = Sv * DEF, with DEF unsquared."""
        check(at_target(rojak_out["ti1"]), expect.vws * expect.total_deformation, label="ti1")

    def test_ti2_uses_era5_divergence(self, rojak_out, expect):
        """Ellrod & Knapp (1992): TI2 = Sv * (DEF - delta).

        `delta` is ERA5's OWN divergence product, not the divergence derived
        from u and v -- rojak's TI2 accepts du_dx and dv_dy but never uses
        them in _compute. On real ERA5 the two are nearly identical (the
        divergence product comes from the same spectral wind), so this is not
        a defect; but the distinction is invisible in the source and a test is
        the right place to record it."""
        check(at_target(rojak_out["ti2"]),
              expect.vws * (expect.total_deformation - expect.divergence), label="ti2")

    def test_ngm1(self, rojak_out, expect):
        check(at_target(rojak_out["ngm1"]), expect.wind_speed * expect.total_deformation,
              label="ngm1")

    def test_ngm2(self, rojak_out, expect):
        """NGM2 = |dT/dz| * DEF.

        TOL_LOOSE, deliberately. T is not linear in pressure, so rojak's
        three-level `.differentiate` carries a real O(dp^2) truncation error --
        about 2.8 % at the project's actual 175/200/225 hPa spacing. That is
        discretisation, not a formula error, and test_convergence.py proves it
        by showing the error fall as the levels are brought closer together."""
        check(at_target(rojak_out["ngm2"]), np.abs(expect.dT_dz) * expect.total_deformation,
              tol=TOL_LOOSE, label="ngm2")


class TestVorticityFamily:
    def test_brown1(self, rojak_out, expect, grid):
        """Brown (1973): sqrt(0.3 zeta_a^2 + D_sh^2 + D_st^2), zeta_a = zeta + f."""
        lat = grid[0]
        f = 2 * OMEGA_EARTH * np.sin(np.deg2rad(lat))[:, None, None]
        zeta_a = expect.vorticity + f
        check(at_target(rojak_out["brown1"]),
              np.sqrt(0.3 * zeta_a ** 2
                      + expect.shearing_deformation ** 2
                      + expect.stretching_deformation ** 2),
              label="brown1")

    def test_nva(self, rojak_out, expect, grid):
        """Sharman (2006): NVA = max(-u d(zeta+f)/dx - v d(zeta+f)/dy, 0).

        f varies with latitude only, so d(f)/dy = beta = 2 Omega cos(phi)/M
        and d(f)/dx = 0. Getting beta wrong here is exactly the bug that was
        found in rojak's UBF, so it is worth exercising in the expectation."""
        lat = grid[0]
        beta = (2 * OMEGA_EARTH * np.cos(np.deg2rad(lat))
                / S.meridional_radius(lat))[:, None, None]
        advection = -(expect.u * expect.dvort_dx + expect.v * (expect.dvort_dy + beta))
        check(at_target(rojak_out["nva"]), np.maximum(advection, 0.0), label="nva")


# ===========================================================================
# Coverage guard
# ===========================================================================
def test_every_diagnostic_is_covered(diag):
    """Fail if a diagnostic is added to REFERENCE_TABLE without a test.

    Without this, the suite silently stops being complete the moment the
    project grows a 22nd diagnostic -- the same shape of failure as a
    diagnostic being dropped from a dict and averaged over as if it were
    never there."""
    covered = {
        # tier C
        "rva_magnitude", "f2d",
        # tier B
        "wind_speed", "vorticity_squared", "horizontal_divergence", "magnitude_pv",
        "temperature_gradient", "vertical_wind_shear", "endlich", "deformation",
        "ti1", "ti2", "ngm1", "ngm2", "brown1", "nva",
    }
    # The remaining five are covered by 4_verify.py's cross-check against a
    # corrected rojak, which is independent evidence of the same strength.
    cross_checked = {"negative_richardson", "colson_panofsky", "ubf", "ncsu1", "brown2"}

    all_keys = set(diag.REFERENCE_TABLE)
    missing = all_keys - covered - cross_checked
    assert not missing, (
        f"{len(missing)} diagnostic(s) have no verification at all: {sorted(missing)}. "
        f"Add an analytic test here, or a cross-check in 4_verify.py."
    )
    assert len(all_keys) == 21, f"expected 21 diagnostics, REFERENCE_TABLE has {len(all_keys)}"
