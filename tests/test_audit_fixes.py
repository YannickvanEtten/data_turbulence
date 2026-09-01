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

    These tests fix the ALGEBRA between the variants. Which one is correct is
    an empirical question and is deliberately not asserted here.
    """

    def test_default_is_C(self, diag, prepared):
        """The default is C as of 2026-08-30 (FORMULA_AUDIT.md §10.4).

        Williams (2017) Fig. 1 plots this diagnostic on 0..300, anchored at
        zero, in the same figure where Negative Richardson runs -300..0 — so
        the published quantity is a magnitude, not a signed tendency. The
        measured p97/median agrees: A 754, C 22.7, published 13.6.

        Pinned as a test because the default is the thing every production run
        silently inherits, and because reverting it would quietly change what
        504 months of output mean."""
        assert diag.F2D_DEFAULT_VARIANT == "C"
        ds = prepared._dataset
        default = diag.frontogenesis_isentropic(ds)
        explicit = diag.frontogenesis_isentropic(ds, variant="C")
        np.testing.assert_allclose(default.values, explicit.values,
                                   rtol=0, atol=0)

    def test_variant_A_still_reproduces_the_old_behaviour(self, diag, prepared):
        """Switching the default must not remove the ability to regenerate
        anything computed before 2026-08-30. `--f2d-variant A` is the escape
        hatch that makes the change reversible, so it is worth a test of its
        own: C is the absolute value of A, exactly."""
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
