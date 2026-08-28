"""
test_calibration_roundtrip.py
=============================
The calibration -> aggregation identity, tested without any ERA5.

WHY THIS TEST EXISTS
--------------------
Prosser's severity thresholds are percentiles OF THE CALIBRATION DISTRIBUTION
ITSELF. Applying them back to the calibration data is therefore not an
approximation to be judged by eye -- it is an identity:

    cos(phi)-weighted exceedance rate of the p-th percentile threshold
        ==  (100 - p) / 100

    LOG  = p97.0  -> 3.00 %      MSOG = p99.8 -> 0.20 %
    LMOG = p99.1  -> 0.90 %      SOG  = p99.9 -> 0.10 %
    MOG  = p99.6  -> 0.40 %

Percentile source: Williams (2017) Table 1, adopted verbatim by Prosser et al.
(2023): "the diagnostic values corresponding to the 97th, 99.1st, 99.6th,
99.8th, and 99.9th percentiles were then derived globally for the reference
year 2000".

This identity is the cheapest possible check on the whole Layer-2/3/4 chain.
It catches, before a single byte of ERA5 is downloaded:

  * a comparison in the wrong direction (>= vs <=)
  * latitude weighting applied in calibration but not in aggregation, or twice
  * a percentile convention mismatch ([0,1] vs [0,100])
  * an axis collapsed in the wrong order
  * NaN cells silently counted as non-exceedances, biasing every rate low

It does NOT check that the diagnostics are physically right -- that is what
tests/ (analytic) and 4_verify.py (cross-implementation) are for. This checks
that the calibration plumbing carries them correctly.

The synthetic fields are deliberately LATITUDE-DEPENDENT. A field that is iid
in latitude would let a broken cos(phi) weighting pass unnoticed;
test_weighting_actually_bites asserts the weighting changes the answer.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import aggregate  # noqa: E402
import calibration  # noqa: E402
from calib_weighted_percentile import weighted_percentile  # noqa: E402

# Williams (2017) Table 1 / Prosser (2023) section 2.
SEVERITIES: dict[str, float] = {
    "light":              97.0,
    "light_to_moderate":  99.1,
    "moderate":           99.6,
    "moderate_to_severe": 99.8,
    "severe":             99.9,
}
ORDER = list(SEVERITIES)

# Hazen plotting positions put the empirical quantile at the midpoint of a
# value's weight slot, so the recovered rate is off by ~0.5/(n*(1-p)) -- about
# 0.07 % for the severe threshold at this grid size. 1 % is a comfortable
# ceiling that is still far tighter than any real bug would produce.
TOLERANCE = 0.01


@pytest.fixture(scope="module")
def synthetic():
    """A synthetic 'global year' on a 2-degree grid, 60 timesteps."""
    rng = np.random.default_rng(20260828)
    lat = np.arange(-90.0, 90.01, 2.0)
    lon = np.arange(-180.0, 180.0, 2.0)
    time = np.arange(60)
    shape = (len(time), len(lat), len(lon))
    coords = {"time": time, "latitude": lat, "longitude": lon}
    dims = ("time", "latitude", "longitude")

    lat_weights = xr.DataArray(
        np.cos(np.deg2rad(lat)),
        coords={"latitude": lat}, dims=("latitude",), name="lat_weights",
    )
    lat2d = np.broadcast_to(lat[None, :, None], shape)

    fields = {
        # lognormal with a tail that thickens toward the poles -- the case
        # where cos(phi) weighting matters most
        "diag_lognormal": xr.DataArray(
            np.exp(rng.normal(0.0, 0.5 + 0.4 * np.abs(lat2d) / 90.0, shape)),
            coords=coords, dims=dims, name="diag_lognormal"),
        # positive, latitude-varying mean
        "diag_gamma": xr.DataArray(
            rng.gamma(2.0, 1.0 + np.cos(np.deg2rad(lat2d)), shape),
            coords=coords, dims=dims, name="diag_gamma"),
        # NEGATIVE-valued with sign="+": the colson_panofsky /
        # negative_richardson shape, where the turbulent tail is the upper
        # (least negative) one and the whole threshold ladder is negative
        "diag_negative": xr.DataArray(
            -np.exp(rng.normal(3.0, 0.6, shape)),
            coords=coords, dims=dims, name="diag_negative"),
        # ~2 % NaN, to exercise the skipna / populated_count path
        "diag_with_nan": xr.DataArray(
            np.where(rng.random(shape) < 0.02, np.nan,
                     rng.normal(10.0, 3.0, shape)),
            coords=coords, dims=dims, name="diag_with_nan"),
    }
    reference_table = {
        k: {"sign": "+", "num": i + 1, "units": "-", "name": k, "wj_median": 1.0}
        for i, k in enumerate(fields)
    }
    thresholds, signs, sample_sizes = calibration.compute_thresholds(
        calibration_fields=fields,
        lat_weights=lat_weights,
        severities=SEVERITIES,
        reference_table=reference_table,
        weighted_percentile=weighted_percentile,
    )
    weights_3d = lat_weights.broadcast_like(fields["diag_lognormal"])
    return dict(fields=fields, thresholds=thresholds, signs=signs,
                sample_sizes=sample_sizes, lat_weights=lat_weights,
                weights_3d=weights_3d, reference_table=reference_table)


def _weighted_rate(exceed: xr.DataArray, weights_3d: xr.DataArray) -> float:
    """cos(phi)-weighted mean of a 0/1/NaN exceedance field.

    NaN cells are excluded from BOTH numerator and denominator, matching
    aggregate.exceedance_field's deliberate NaN propagation (Q-INTEG-3).
    """
    populated = exceed.notnull()
    return float(
        (exceed.fillna(0.0) * weights_3d).sum() / weights_3d.where(populated).sum()
    )


@pytest.mark.parametrize("severity", ORDER)
@pytest.mark.parametrize(
    "name", ["diag_lognormal", "diag_gamma", "diag_negative", "diag_with_nan"]
)
def test_threshold_reproduces_its_own_percentile(synthetic, name, severity):
    """The identity, per diagnostic and per severity."""
    expected = (100.0 - SEVERITIES[severity]) / 100.0
    exceed = aggregate.exceedance_field(
        synthetic["fields"][name],
        synthetic["thresholds"][name][severity],
        synthetic["signs"][name],
    )
    observed = _weighted_rate(exceed, synthetic["weights_3d"])
    assert observed == pytest.approx(expected, rel=TOLERANCE), (
        f"{name} @ {severity}: p{SEVERITIES[severity]} threshold should be "
        f"exceeded by {expected:.4%} of the calibration data, got {observed:.4%}"
    )


@pytest.mark.parametrize("severity", ORDER)
def test_ensemble_mean_reproduces_percentile(synthetic, severity):
    """Same identity through the real entry point.

    The diagnostic-mean exceedance field must also average to (100-p)/100,
    because every constituent field does.
    """
    out = aggregate.exceedance_mean_all_severities(
        synthetic["fields"], synthetic["thresholds"], ORDER, synthetic["signs"]
    )
    observed = _weighted_rate(out[severity]["exceedance_mean"],
                              synthetic["weights_3d"])
    expected = (100.0 - SEVERITIES[severity]) / 100.0
    assert observed == pytest.approx(expected, rel=TOLERANCE)


def test_populated_count_sees_the_gaps(synthetic):
    """A diagnostic with NaNs must reduce populated_count, not vote 'no'."""
    out = aggregate.exceedance_mean_all_severities(
        synthetic["fields"], synthetic["thresholds"], ORDER, synthetic["signs"]
    )
    count = out["moderate"]["populated_count"]
    assert int(count.max()) == len(synthetic["fields"])
    assert int(count.min()) == len(synthetic["fields"]) - 1, (
        "the 2%-NaN diagnostic should drop the count by exactly one where "
        "it is gapped -- if the count never drops, NaN is being silently "
        "counted as a non-exceedance"
    )


def test_weighting_actually_bites(synthetic):
    """cos(phi) weighting must change the threshold on a latitude-dependent
    field, and in the right direction.

    diag_lognormal's tail thickens toward the poles. Down-weighting the
    over-sampled polar cells must therefore LOWER the threshold. If weighted
    and unweighted agree, the weights are not reaching the percentile.
    """
    lat = synthetic["lat_weights"]["latitude"]
    flat = xr.DataArray(np.ones(lat.size), coords={"latitude": lat},
                        dims=("latitude",))
    unweighted, _, _ = calibration.compute_thresholds(
        calibration_fields=synthetic["fields"],
        lat_weights=flat,
        severities=SEVERITIES,
        reference_table=synthetic["reference_table"],
        weighted_percentile=weighted_percentile,
    )
    w = synthetic["thresholds"]["diag_lognormal"]["moderate"]
    u = unweighted["diag_lognormal"]["moderate"]
    assert abs(w - u) / abs(u) > 0.05, (
        f"weighted ({w:.4f}) and unweighted ({u:.4f}) thresholds are nearly "
        f"identical on a latitude-dependent field -- the cos(phi) weights are "
        f"not reaching weighted_percentile"
    )
    assert w < u, (
        "diag_lognormal is heavier-tailed toward the poles, so down-weighting "
        "polar cells must lower the threshold"
    )


def test_thresholds_are_monotonic_in_severity(synthetic):
    """light <= light_to_moderate <= ... <= severe, for every diagnostic.

    True even for the negative-valued ladder: Williams (2017) Table 2's
    Colson-Panofsky row runs -29.3 -> -22.2, i.e. monotonically INCREASING
    toward severe. A ladder that is not monotonic means the percentiles were
    applied to the wrong tail.
    """
    for name, per_sev in synthetic["thresholds"].items():
        ladder = [per_sev[s] for s in ORDER]
        assert ladder == sorted(ladder), f"{name} ladder is not monotonic: {ladder}"


def test_empty_calibration_field_is_refused(synthetic):
    """A diagnostic with no finite values must raise, not produce a threshold.

    This is the f2d-on-a-2-timestep-file failure mode (STATUS.md 4a) reaching
    the calibration layer.
    """
    fields = dict(synthetic["fields"])
    fields["diag_empty"] = xr.full_like(fields["diag_gamma"], np.nan)
    ref = dict(synthetic["reference_table"])
    ref["diag_empty"] = {"sign": "+", "num": 99, "units": "-",
                         "name": "diag_empty", "wj_median": 1.0}
    with pytest.raises(ValueError, match="no finite values"):
        calibration.compute_thresholds(
            calibration_fields=fields,
            lat_weights=synthetic["lat_weights"],
            severities=SEVERITIES,
            reference_table=ref,
            weighted_percentile=weighted_percentile,
        )


def test_hardcoding_guard_fires_on_literal_table2():
    """Q-CALIB-6: copy-pasted Williams (2017) Table 2 must be refused."""
    literal = {"negative_richardson":
               dict(zip(ORDER, [-15.4, -9.8, -7.9, -6.7, -5.9]))}
    with pytest.raises(RuntimeError, match="Q-CALIB-6"):
        aggregate.assert_thresholds_not_hardcoded_table2(literal, ORDER)


def test_hardcoding_guard_passes_a_real_calibration():
    """...but a genuine calibration landing near the table must be allowed."""
    close = {"negative_richardson":
             dict(zip(ORDER, [-15.7, -9.9, -8.0, -6.8, -6.0]))}
    aggregate.assert_thresholds_not_hardcoded_table2(close, ORDER)
