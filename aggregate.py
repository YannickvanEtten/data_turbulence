"""
aggregate.py
=============
Q-AGG-1: exceedance-then-average aggregation layer.

Prosser's aggregation, in order:
    1. For each of the 21 diagnostics and each severity level, compute a
       binary 0/1 exceedance field per cell/timestep: 1 if the diagnostic
       value at that cell/timestep crosses that severity's threshold, 0
       otherwise.
    2. Average the 21 binary exceedance fields together -> one
       diagnostic-mean exceedance field per severity level.

This is NOT the same as averaging the 21 raw diagnostic values (or
z-scores) together first and thresholding the averaged field once --
the two orderings are non-commutative in general, because averaging
raw values first can wash out a spatially/temporally localized signal
that only a subset of diagnostics actually flagged, while exceedance
counting preserves "how many diagnostics agree this cell is turbulent"
regardless of how large any one diagnostic's raw magnitude is.

This module implements ONLY the exceedance-first (correct) ordering. The
averaging-first ordering is implemented too, but only inside the test
file, as the deliberately-wrong comparison case.

--- Q-AGG-3: sign="either" (colson_panofsky #2, negative_richardson #11) ---

REFERENCE_TABLE marks these two diagnostics sign="either", implying a
genuinely two-tailed turbulence criterion. Checked directly against
Williams (2017) Table 2 (the actual severity-threshold table Prosser's
methodology descends from) rather than trusting the "either" tag at face
value:

    Negative Richardson number: -15.4  -9.8  -7.9  -6.7  -5.9   (light->severe)
    Colson-Panofsky index:      -29.3 -27.0 -25.2 -23.7 -22.2   (light->severe)

Both are single ORDERED LADDERS of thresholds, monotonically increasing
(becoming less negative) from light to severe -- not two side-by-side
threshold ladders (one per tail), which is what a genuinely two-tailed
diagnostic's calibration table would need to look like. This is exactly
the single-tail signature: light turbulence is the least extreme (most
negative) threshold, severe is the closest to zero/most positive. It
matches the physics directly -- Sharman (2006) section on Ri: "regions of
small Ri should be favored regions of turbulence" is a ONE-SIDED
statement (turbulence favored as Ri gets SMALL, not as it gets small OR
large), and Colson-Panofsky = Sv^2*dz^2*(1-Ri/Ri_crit) increases (toward
positive) as Ri decreases below Ri_crit for the same reason.

**Verdict: neither diagnostic is actually two-tailed.** Both are sign="+"
(value >= threshold), just with negative-valued thresholds -- a negative
MEDIAN does not imply a two-tailed CRITERION. REFERENCE_TABLE's "either"
tag for #2 and #11 appears to be a mislabeling and should be corrected to
"+" with the existing scalar-threshold code path (no threshold-tuple
needed for either of these two).

sign="either" is still implemented below as a real, correctly-designed
capability (per the ask), for any diagnostic that turns out to genuinely
need it in the future -- it is NOT collapsed to abs(value)>=abs(threshold)
(forcing symmetry), since a real two-tailed diagnostic's two tails are not
guaranteed to have the same severity ladder in each direction; threshold
is a (lower, upper) tuple instead.

--- Q-CALIB-6: thresholds MUST be calibration-derived, never Table 2 literal ---

The whole units/sign audit this project has done (Q-AGG-3, Q-UNITS-1, and
the broader REFERENCE_TABLE corrections) rests on an assumption that is
easy to silently violate: that severity thresholds fed into this module
come from Prosser's percentile-calibration pipeline (Q-CALIB-1 through
Q-CALIB-5), NOT copy-pasted directly from Williams (2017) Table 2. Table 2's
literal numbers were derived from a DIFFERENT (older, coarser) climatology
than this project's ERA5-based calibration -- using them directly would
silently reintroduce every one of the scale/units issues this session's
audit was meant to close, without any error or warning, because the
numbers would still "look plausible" (same order of magnitude, same sign).

_KNOWN_TABLE2_LITERAL_VALUES below is a small, explicitly-sourced snapshot
of the few Table 2 rows this project has directly quoted during that audit
(Q-AGG-3). assert_thresholds_not_hardcoded_table2() is called automatically
by exceedance_mean_all_severities() on every invocation -- it is not a
comment or a docstring warning, it is a runtime check that raises if any
diagnostic's supplied thresholds match the literal published table to
within a tight tolerance. A real calibration run should NEVER produce an
exact (or near-exact) match to Table 2's numbers; if it ever does, either
the calibration pipeline has a bug, or someone has hardcoded the table as
a shortcut. Either way this should fail loudly, immediately, not silently
propagate into the trend-regression output.
"""
from __future__ import annotations

from typing import Literal, Union

import numpy as np
import xarray as xr

Sign = Literal["+", "-", "either"]
# For sign="+"/"-": a single scalar threshold.
# For sign="either": a (lower, upper) tuple -- exceedance fires if the value
# is <= lower OR >= upper. NOT forced to be symmetric (i.e. NOT abs(value) >=
# abs(threshold)) -- see the Q-AGG-3 docstring below for why.
Threshold = Union[float, tuple[float, float]]


# ---------------------------------------------------------------------------
# Q-CALIB-6: known Williams (2017) Table 2 LITERAL values, for detecting
# accidental hardcoding. Source: directly quoted during the Q-AGG-3 sign
# audit this session. light -> severe order. Extend this dict if more
# Table 2 rows get directly quoted/verified elsewhere in the project --
# it only needs enough entries to catch a copy-paste, not the full table.
# ---------------------------------------------------------------------------
_KNOWN_TABLE2_LITERAL_VALUES: dict[str, list[float]] = {
    "negative_richardson": [-15.4, -9.8, -7.9, -6.7, -5.9],
    "colson_panofsky":     [-29.3, -27.0, -25.2, -23.7, -22.2],
}


def assert_thresholds_not_hardcoded_table2(
    thresholds: dict[str, dict[str, Threshold]],
    severities: list[str],
    rtol: float = 1e-4,
) -> None:
    """Q-CALIB-6: fail loudly if any diagnostic's thresholds match the
    literal published Williams (2017) Table 2 values, instead of coming
    from this project's own percentile calibration.

    Checks each diagnostic in thresholds that also appears in
    _KNOWN_TABLE2_LITERAL_VALUES. If ALL of that diagnostic's supplied
    severity thresholds match the literal table (within rtol) in the same
    order, raises RuntimeError -- a real calibration run on a different
    climatology should not reproduce the published numbers this closely.

    Only checks sign="+"/"-"' scalar thresholds against the literal list
    (both known entries are scalar, sign="+", per the Q-AGG-3 fix) --
    a sign="either" (lower, upper) tuple threshold can't literal-match a
    single-value table row and is skipped.
    """
    for name, literal_values in _KNOWN_TABLE2_LITERAL_VALUES.items():
        if name not in thresholds:
            continue
        try:
            supplied = [thresholds[name][sev] for sev in severities]
        except KeyError:
            continue  # doesn't cover all the same severities -- can't compare fairly
        if len(supplied) != len(literal_values):
            continue
        if any(isinstance(v, tuple) for v in supplied):
            continue  # sign="either" tuple thresholds can't match a scalar literal row
        if np.allclose(supplied, literal_values, rtol=rtol):
            raise RuntimeError(
                f"Q-CALIB-6 VIOLATION: thresholds for '{name}' exactly match the "
                f"literal Williams (2017) Table 2 values {literal_values} instead "
                f"of coming from this project's own percentile calibration "
                f"(Q-CALIB-1..5). This is almost certainly a hardcoded shortcut, "
                f"not a real calibration result -- a different climatology should "
                f"not reproduce the published table to {rtol:.0e} relative "
                f"tolerance. Refusing to proceed."
            )


def exceedance_field(diagnostic: xr.DataArray, threshold: Threshold, sign: Sign = "+") -> xr.DataArray:
    """Binary 0/1 exceedance field for a single diagnostic at a single
    severity threshold.

    sign="+"      -> exceedance means diagnostic value >= threshold
                     (most W&J diagnostics: larger magnitude = more turbulent)
    sign="-"      -> exceedance means diagnostic value <= threshold
                     (diagnostics whose turbulent regime is the negative tail)
    sign="either" -> Q-AGG-3: genuinely two-tailed. threshold is a
                     (lower, upper) tuple; exceedance means value <= lower
                     OR value >= upper. Deliberately NOT collapsed to
                     abs(value) >= abs(threshold) -- see module docstring
                     for Q-AGG-3: the two tails are not required to be
                     symmetric, and forcing symmetry would silently produce
                     the wrong threshold for a diagnostic whose two tails
                     have different physical severity ladders.

    Q-INTEG-3 (Bug 1, decided option (b), skipna): NaN input is
    DELIBERATELY propagated to NaN output, never silently coerced to a
    boolean. IEEE 754 makes `nan >= threshold`, `nan <= threshold`, etc.
    all evaluate False -- left alone, a masked/gap cell would silently
    read as "does not exceed" and bias every downstream mean low with no
    trace. This function is the single choke point all exceedance
    counting flows through, so the propagate-vs-False decision belongs
    here, not repeated (and possibly gotten wrong) at every call site.
    Downstream aggregation (`exceedance_mean_single_severity`) is what
    turns this NaN into "excluded from this cell's mean, and the
    populated-count drops by one" -- this function's only job is to make
    sure the NaN survives long enough to be excluded on purpose, not by
    accident.
    """
    if sign == "+":
        exceeds = diagnostic >= threshold
    elif sign == "-":
        exceeds = diagnostic <= threshold
    elif sign == "either":
        if not (isinstance(threshold, tuple) and len(threshold) == 2):
            raise ValueError(
                f"sign='either' requires threshold=(lower, upper), got {threshold!r}"
            )
        lower, upper = threshold
        if lower > upper:
            raise ValueError(
                f"sign='either' threshold=(lower, upper) must have lower <= upper, "
                f"got lower={lower}, upper={upper}"
            )
        exceeds = (diagnostic <= lower) | (diagnostic >= upper)
    else:
        raise ValueError(f"sign must be '+', '-', or 'either', got {sign!r}")
    # Q-INTEG-3 (Bug 1): re-mask with the input's own NaN pattern. Every
    # branch above independently produces False at NaN cells (IEEE 754);
    # this line is what turns that False back into NaN before it can be
    # silently counted as a non-exceedance.
    exceeds = exceeds.astype(np.float64).where(diagnostic.notnull())
    return exceeds.rename(f"{diagnostic.name}_exceeds")


def exceedance_mean_single_severity(
    diagnostics: dict[str, xr.DataArray],
    thresholds: dict[str, Threshold],
    signs: dict[str, Sign] | None = None,
) -> xr.Dataset:
    """Step 1 + Step 2 for ONE severity level.

    diagnostics: {diagnostic_name: DataArray}, all on the same grid
                 (same cell/timestep coordinates).
    thresholds:  {diagnostic_name: threshold value for this severity} --
                 a scalar for sign="+"/"-", a (lower, upper) tuple for
                 sign="either" (Q-AGG-3).
    signs:       {diagnostic_name: "+"/"-"/"either"}; defaults to "+" for
                 any diagnostic not listed.

    Returns an xr.Dataset with two per-cell/timestep variables:
      - "exceedance_mean":    Q-INTEG-3 (Bug 1, option (b)) -- the MEAN of
        the binary exceedance fields across all diagnostics SUPPLIED,
        skipping any diagnostic that is NaN at that cell (`skipna=True`).
        This is a populated-only mean, not a fixed-21 mean: a cell with
        one gapped diagnostic is averaged over the remaining N-1, not
        silently treated as if the gapped one voted "no."
      - "populated_count":    per-cell/timestep count of how many of the
        `len(diagnostics)` supplied diagnostics actually contributed a
        valid (non-NaN) value at that cell. This is a first-class output,
        not a diagnostic side-channel -- a cell built from 15/21
        diagnostics must stay visibly distinguishable downstream from one
        built from 21/21, rather than the denominator changing silently.
        Callers that want a fixed all-21 semantics instead (Q-INTEG-3
        option (a), NaN-propagate) can derive it trivially from this
        Dataset: `mean.where(count == len(diagnostics))`.

    Exceedance is still computed first, per diagnostic, THEN averaged --
    the ordering Prosser uses (Q-AGG-1) is unchanged by this fix.
    """
    signs = signs or {}
    exceed_fields = [
        exceedance_field(da, thresholds[name], signs.get(name, "+"))
        for name, da in diagnostics.items()
    ]
    stacked = xr.concat(exceed_fields, dim="diagnostic")
    mean = stacked.mean(dim="diagnostic", skipna=True).rename("exceedance_mean")
    populated_count = (
        stacked.notnull().sum(dim="diagnostic").rename("populated_count")
    )
    return xr.Dataset({"exceedance_mean": mean, "populated_count": populated_count})


def exceedance_mean_all_severities(
    diagnostics: dict[str, xr.DataArray],
    thresholds: dict[str, dict[str, Threshold]],
    severities: list[str],
    signs: dict[str, Sign] | None = None,
) -> dict[str, xr.Dataset]:
    """Full Q-AGG-1 aggregation: all 21 diagnostics x all severity levels.

    thresholds: {diagnostic_name: {severity_name: threshold_value}}
                e.g. {"ubf": {"light": 1e-9, "moderate": 5e-9, ...}, ...}
                threshold_value is a scalar, EXCEPT for sign="either"
                diagnostics, where it's a (lower, upper) tuple.
    severities: ordered list of severity names to produce, e.g.
                ["light", "light_to_moderate", "moderate",
                 "moderate_to_severe", "severe"]

    Returns {severity_name: Dataset(exceedance_mean, populated_count)} --
    see exceedance_mean_single_severity for the Q-INTEG-3 (Bug 1) fields.
    """
    assert_thresholds_not_hardcoded_table2(thresholds, severities)  # Q-CALIB-6: mandatory, not optional
    out: dict[str, xr.Dataset] = {}
    for severity in severities:
        per_diag_thresholds = {name: thr[severity] for name, thr in thresholds.items()}
        out[severity] = exceedance_mean_single_severity(diagnostics, per_diag_thresholds, signs)
    return out


# ---------------------------------------------------------------------------
# WRONG ordering, implemented here ONLY for the divergence test below --
# average raw diagnostic values first, threshold the averaged field once.
# Not exported for use anywhere else in the pipeline.
# ---------------------------------------------------------------------------
def _average_then_threshold_WRONG(
    diagnostics: dict[str, xr.DataArray],
    thresholds: dict[str, float],
    signs: dict[str, Sign] | None = None,
) -> xr.DataArray:
    signs = signs or {}
    # Averaging raw values across diagnostics with different units/scales
    # is itself questionable, but the point here is ordering, not units --
    # use z-score-free raw averaging to isolate the ordering effect only.
    stacked = xr.concat(list(diagnostics.values()), dim="diagnostic")
    mean_raw = stacked.mean(dim="diagnostic")
    # thresholding the AVERAGED field once, using the mean of the per-
    # diagnostic thresholds as the single cutoff (again: isolating the
    # ordering effect, not conflating it with a units-mismatch bug)
    mean_threshold = float(np.mean(list(thresholds.values())))
    default_sign = list(signs.values())[0] if signs else "+"
    if default_sign == "+":
        exceeds = mean_raw >= mean_threshold
    else:
        exceeds = mean_raw <= mean_threshold
    return exceeds.astype(np.float64).rename("wrong_averaged_exceedance")
