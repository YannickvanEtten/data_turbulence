"""
tests/test_audit_fixes.py
=========================
Pins the two code changes made from FORMULA_AUDIT.md on 2026-08-29.

Neither change is visible in the replication results -- that is precisely why
each needs a test. The deformation fix is invisible because squaring is a
monotone transform and the exceedance field is identical (STATUS.md §5.5); the
f2d variant flag is invisible because its default reproduces the previous
behaviour exactly. Both would therefore drift silently.

No ERA5 required: everything here runs on the manufactured atmosphere in
tests/synthetic.py, like the rest of the suite.
"""
from __future__ import annotations

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# #8 DEFORMATION -- compute_all_21 must un-square what rojak returns
# ---------------------------------------------------------------------------
class TestDeformationUnsquared:
    """FORMULA_AUDIT.md §5.

    rojak's DEF diagnostic returns DEF^2. Sharman A17 and every published table
    define DEF = (D_SH^2 + D_ST^2)^(1/2). Until 2026-08-29 the square root was
    applied only in 3_pipeline.py's comparison-table path, so every zarr
    written by ada/diagnostics_global.py held DEF^2 under the name
    `deformation`.

    The two assertions below have to move together, and that is the point of
    putting them in one class: un-squaring in compute_all_21 while LEAVING the
    square root in 3_pipeline would take the fourth root and the comparison
    table would drift, which is the same failure the old
    test_deformation_is_squared docstring warned about from the other side.
    """

    def test_compute_all_21_returns_unsquared_deformation(self, diag, prepared):
        out, failures = diag.compute_all_21(prepared)
        assert not failures, f"diagnostics failed: {[f.key for f in failures]}"
        assert "deformation" in out

        rojak_raw, _ = diag.compute_rojak_diagnostics(prepared)
        raw = rojak_raw["deformation"]
        if "pressure_level" in raw.dims:
            raw = raw.sel(pressure_level=200)

        got = np.asarray(out["deformation"].values, dtype=np.float64)
        want = np.sqrt(np.abs(np.asarray(raw.values, dtype=np.float64)))
        ok = np.isfinite(got) & np.isfinite(want) & (want > 0)
        assert ok.any(), "no finite deformation values to compare"
        rel = np.abs(got[ok] - want[ok]) / want[ok]
        assert float(np.median(rel)) < 1e-6, (
            f"compute_all_21's deformation is not sqrt(rojak's): "
            f"median relative error {float(np.median(rel)):.3e}"
        )

    def test_compute_all_21_deformation_records_the_convention(self, diag, prepared):
        """The attribute is how a zarr written later can be told apart from one
        written before the fix. A file with no such attribute is old."""
        out, _ = diag.compute_all_21(prepared)
        attrs = out["deformation"].attrs
        assert attrs.get("units") == "s-1"
        assert attrs.get("sharman_eq") == "A17"

    def test_pipeline_does_not_also_take_a_square_root(self):
        """3_pipeline.PRETRANSFORM must NOT carry a deformation entry any more.

        If it does, the comparison table takes the fourth root of DEF^2 and
        every deformation ratio in cat_outputs/comparison_table.csv is wrong by
        a square root -- while still looking entirely plausible.
        """
        from conftest import _load
        pipeline = _load("pipeline", "3_pipeline.py")
        assert "deformation" not in pipeline.PRETRANSFORM, (
            "3_pipeline.PRETRANSFORM still square-roots deformation, but "
            "compute_all_21 already did. See FORMULA_AUDIT.md §5."
        )


# ---------------------------------------------------------------------------
# #20 F2D -- the four readings of Sharman A9
# ---------------------------------------------------------------------------
class TestF2dVariants:
    """FORMULA_AUDIT.md §4.

    A9 as printed is internally inconsistent: its left-hand side carries a
    leading minus its right-hand side does not, and it carries a
    |dv/dtheta|^-1 normalisation the implementation never had. Rather than
    guess, all four readings are selectable and ada/check_f2d_variants.py
    measures which reproduces the published distribution shape.

    These tests fix the ALGEBRA between the variants. WHICH ONE IS CORRECT WAS
    an empirical question when this class was written, and is not any more:
    Williams & Storer (2022) Eq. (3), p. 1427 writes the diagnostic as the
    signed D/Dt|du/dtheta|^2, i.e. variant A. The default assertion below
    tracks that citation; the algebra assertions are independent of it.
    """

    def test_default_is_A(self, diag, prepared):
        """The default is A as of 2026-09-08, superseding the 2026-08-30 C.

        This assertion was `== "C"`, justified by Williams (2017) Fig. 1's
        0..300 axis and a measured p97/median of 22.7 against a published 13.6.
        Both are inferences from a figure, and they were the best available
        while no equation for this diagnostic could be found anywhere in the
        lineage.

        Williams & Storer (2022), Q. J. R. Meteorol. Soc. 148, 1424-1438,
        Eq. (3), p. 1427 states it directly:

            F_theta = D/Dt |du/dtheta|^2

        signed, un-normalised, no absolute value, no clip -- variant A up to a
        constant 1/2. Same author, one year before Prosser (2023), who lists
        Williams as a co-author. A written equation from the replicated lineage
        outranks an inference from an axis, so the pin moved.

        The counter-evidence is NOT resolved by this and must not be dropped:
        Fig. 1's axis really is anchored at zero, and a p97/median of 13.6 is
        not something a zero-centred material derivative produces. That is a
        conflict between two Williams papers. It is recorded in
        AUDIT_diagnostics_vs_literature.md 6 F6 and belongs in any write-up.

        Pinned as a test because the default is what every production run
        silently inherits, and because changing it changes what 504 months of
        output mean. `--f2d-variant C` regenerates the pre-2026-09-08 archive
        exactly."""
        assert diag.F2D_DEFAULT_VARIANT == "A"
        ds = prepared._dataset
        default = diag.frontogenesis_isentropic(ds)
        explicit = diag.frontogenesis_isentropic(ds, variant="A")
        np.testing.assert_allclose(default.values, explicit.values,
                                   rtol=0, atol=0)

    def test_variant_C_still_reproduces_the_archived_run(self, diag, prepared):
        """Switching the default must not remove the ability to regenerate the
        archive. The direction of travel reversed on 2026-09-08: C was the
        default and A the escape hatch, and now A is the default and C is the
        escape hatch for every zarr written between 2026-08-30 and 2026-09-08.
        Either way the relationship is the same and exact -- C is |A| -- so one
        test covers both directions."""
        ds = prepared._dataset
        a = diag.frontogenesis_isentropic(ds, variant="A").values
        c = diag.frontogenesis_isentropic(ds, variant="C").values
        np.testing.assert_allclose(c, np.abs(a), rtol=0, atol=0)

    def test_B_is_exactly_minus_A(self, diag, prepared):
        ds = prepared._dataset
        a = diag.frontogenesis_isentropic(ds, variant="A").values
        b = diag.frontogenesis_isentropic(ds, variant="B").values
        np.testing.assert_allclose(b, -a, rtol=0, atol=0)

    def test_C_is_exactly_abs_A(self, diag, prepared):
        ds = prepared._dataset
        a = diag.frontogenesis_isentropic(ds, variant="A").values
        c = diag.frontogenesis_isentropic(ds, variant="C").values
        np.testing.assert_allclose(c, np.abs(a), rtol=0, atol=0)

    def test_A_and_B_flag_disjoint_tails(self, diag, prepared):
        """The consequence that makes the sign question matter.

        A and B have mirror-image distributions, so every symmetric statistic
        agrees and no magnitude comparison can tell them apart. What differs is
        WHICH CELLS a one-tailed percentile threshold selects -- and those sets
        are essentially disjoint. That is why the choice cannot be left to a
        default, and why FORMULA_AUDIT.md §4 treats it as a §5 category-2
        error rather than a cosmetic one.
        """
        ds = prepared._dataset
        a = np.asarray(diag.frontogenesis_isentropic(ds, variant="A").values,
                       dtype=np.float64).ravel()
        b = -a
        finite = np.isfinite(a)
        a, b = a[finite], b[finite]
        if a.size < 100:
            pytest.skip("synthetic atmosphere too small for a tail comparison")
        ta, tb = np.quantile(a, 0.99), np.quantile(b, 0.99)
        ma, mb = a >= ta, b >= tb
        overlap = np.count_nonzero(ma & mb) / max(np.count_nonzero(ma | mb), 1)
        assert overlap < 0.10, (
            f"A and B flag overlapping cells (Jaccard {overlap:.3f}); expected "
            f"near-disjoint tails"
        )

    def test_D_has_different_units(self, diag, prepared):
        """Variant D is literal A9 including the normalisation, so it is a
        different physical quantity: m s^-2 K^-1, not m^2 s^-3 K^-2. That is
        exactly why Sharman's own Table B1 units argue against it, and the
        attribute has to say so or a later magnitude comparison will silently
        use the wrong published row."""
        ds = prepared._dataset
        d = diag.frontogenesis_isentropic(ds, variant="D")
        a = diag.frontogenesis_isentropic(ds, variant="A")
        assert d.attrs["units"] == "m s-2 K-1"
        assert a.attrs["units"] == "m2 s-3 K-2"

    def test_variant_recorded_in_attrs(self, diag, prepared):
        ds = prepared._dataset
        for v in "ABCD":
            out = diag.frontogenesis_isentropic(ds, variant=v)
            assert out.attrs["f2d_variant"] == v

    def test_unknown_variant_raises(self, diag, prepared):
        with pytest.raises(ValueError, match="unknown f2d variant"):
            diag.frontogenesis_isentropic(prepared._dataset, variant="Z")

    def test_compute_all_21_threads_the_variant(self, diag, prepared):
        out_a, _ = diag.compute_all_21(prepared, f2d_variant="A")
        out_b, _ = diag.compute_all_21(prepared, f2d_variant="B")
        np.testing.assert_allclose(out_b["f2d"].values, -out_a["f2d"].values,
                                   rtol=0, atol=0)
        assert out_b["f2d"].attrs["f2d_variant"] == "B"


# ---------------------------------------------------------------------------
# F12 -- the four-variable provenance substitution
# ---------------------------------------------------------------------------
class TestFourVariableProvenance:
    """AUDIT §5.5 generalised; STATUS.md §18.

    Prosser (2023) p. 2 downloads u, v, T and z and nothing else, so every
    vorticity, divergence and potential vorticity in his 21 is a finite
    difference of that wind. This project downloads seven fields and takes
    ERA5's archived `vo`, `d` and `pv` because CDS offers them. F12 removes
    that difference in one place.

    These tests do NOT assert that the substitution is more correct. It is
    not: ERA5's archived zeta and delta are the IFS's own prognostic spectral
    state (IFS CY49R1 Part III §2.2.7) and u, v are derived from THEM. What
    they assert is that the flag does exactly what it claims -- one operator,
    every consumer, nothing touched when it is off.
    """

    def test_flag_defaults_off_and_baseline_is_unchanged(self, diag):
        assert diag.FixSet().four_variable_provenance is False
        assert diag.BASELINE_FIXES.four_variable_provenance is False
        assert diag.FixSet().label() == "baseline"

    def test_label_and_attrs_carry_the_flag(self, diag):
        f = diag.FixSet(four_variable_provenance=True)
        assert f.label() == "four_variable_provenance"
        assert f.as_attrs()["audit_fix_four_variable_provenance"] == 1
        assert diag.FixSet().as_attrs()["audit_fix_four_variable_provenance"] == 0

    def test_substitution_marks_all_three_fields(self, diag, prepared):
        sub = diag.substitute_computed_fields(prepared)
        ds = sub._dataset
        for name in ("vorticity", "divergence_of_wind", "potential_vorticity"):
            assert ds[name].attrs["field_source"] == "computed_from_u_v_t_z", name
        assert ds.attrs["four_variable_provenance"] == 1

    def test_substituted_fields_keep_shape_and_are_finite(self, diag, prepared):
        ds0 = prepared._dataset
        ds1 = diag.substitute_computed_fields(prepared)._dataset
        for name in ("vorticity", "divergence_of_wind", "potential_vorticity"):
            assert ds1[name].dims == ds0[name].dims, name
            assert ds1[name].shape == ds0[name].shape, name
            assert np.isfinite(np.asarray(ds1[name].values)).any(), name

    def test_zeta_and_delta_come_from_one_call(self, diag, prepared):
        """The point of substituting on the dataset rather than per
        diagnostic: zeta and delta must be mutually consistent, i.e. built
        from the same four velocity derivatives. Checked by rebuilding them
        independently from the same operator and requiring bit equality."""
        from rojak.core.derivatives import vector_derivatives, VelocityDerivative
        ds = prepared._dataset
        lev = float(np.asarray(ds["pressure_level"].values)[0])
        vd = vector_derivatives(
            ds["eastward_wind"].sel(pressure_level=lev),
            ds["northward_wind"].sel(pressure_level=lev), "deg",
            components=[VelocityDerivative.DU_DX, VelocityDerivative.DU_DY,
                        VelocityDerivative.DV_DX, VelocityDerivative.DV_DY])
        expect_zeta = (vd[VelocityDerivative.DV_DX] - vd[VelocityDerivative.DU_DY])
        expect_div = (vd[VelocityDerivative.DU_DX] + vd[VelocityDerivative.DV_DY])
        got = diag.substitute_computed_fields(prepared)._dataset
        np.testing.assert_array_equal(
            np.asarray(got["vorticity"].sel(pressure_level=lev).values),
            np.asarray(expect_zeta.transpose(*got["vorticity"].sel(pressure_level=lev).dims).values))
        np.testing.assert_array_equal(
            np.asarray(got["divergence_of_wind"].sel(pressure_level=lev).values),
            np.asarray(expect_div.transpose(*got["divergence_of_wind"].sel(pressure_level=lev).dims).values))

    def test_off_leaves_the_archived_fields_untouched(self, diag, prepared):
        ds = prepared._dataset
        before = {n: np.array(ds[n].values, copy=True)
                  for n in ("vorticity", "divergence_of_wind")}
        diag.compute_all_21(prepared, fixes=diag.FixSet())
        for n, arr in before.items():
            np.testing.assert_array_equal(np.asarray(ds[n].values), arr)

    def test_flag_changes_the_archived_field_consumers(self, diag, prepared):
        """#5 horizontal_divergence is PURE archived `d` in the baseline
        (STATUS.md §15.6), so it is the sharpest single test that the
        substitution reaches rojak's own diagnostics and not only the
        hand-written seven."""
        base, _ = diag.compute_all_21(prepared, fixes=diag.FixSet())
        var, _ = diag.compute_all_21(
            prepared, fixes=diag.FixSet(four_variable_provenance=True))
        assert set(base) == set(var)
        a = np.asarray(base["horizontal_divergence"].values)
        b = np.asarray(var["horizontal_divergence"].values)
        m = np.isfinite(a) & np.isfinite(b)
        assert m.any()
        assert not np.allclose(a[m], b[m]), (
            "horizontal_divergence is unchanged under F12, which would mean "
            "the substitution never reached rojak's diagnostics")

    def test_attrs_record_the_run(self, diag, prepared):
        var, _ = diag.compute_all_21(
            prepared, fixes=diag.FixSet(four_variable_provenance=True))
        for da in var.values():
            assert da.attrs["audit_fix_four_variable_provenance"] == 1
            assert "four_variable_provenance" in da.attrs["audit_fixes"]
