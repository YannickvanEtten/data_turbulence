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
"""
from __future__ import annotations

from typing import Literal

import numpy as np
import xarray as xr

Sign = Literal["+", "-"]


def exceedance_field(diagnostic: xr.DataArray, threshold: float, sign: Sign = "+") -> xr.DataArray:
    """Binary 0/1 exceedance field for a single diagnostic at a single
    severity threshold.

    sign="+"  -> exceedance means diagnostic value >= threshold
                 (most W&J diagnostics: larger magnitude = more turbulent)
    sign="-"  -> exceedance means diagnostic value <= threshold
                 (diagnostics whose turbulent regime is the negative tail,
                 e.g. negative-Richardson-number-derived indices)
    """
    if sign == "+":
        exceeds = diagnostic >= threshold
    elif sign == "-":
        exceeds = diagnostic <= threshold
    else:
        raise ValueError(f"sign must be '+' or '-', got {sign!r}")
    return exceeds.astype(np.float64).rename(f"{diagnostic.name}_exceeds")


def exceedance_mean_single_severity(
    diagnostics: dict[str, xr.DataArray],
    thresholds: dict[str, float],
    signs: dict[str, Sign] | None = None,
) -> xr.DataArray:
    """Step 1 + Step 2 for ONE severity level.

    diagnostics: {diagnostic_name: DataArray}, all on the same grid
                 (same cell/timestep coordinates).
    thresholds:  {diagnostic_name: threshold value for this severity}.
    signs:       {diagnostic_name: "+"/"-"}; defaults to "+" for any
                 diagnostic not listed.

    Returns the per-cell/timestep MEAN of the binary exceedance fields
    across all diagnostics supplied -- exceedance computed first, per
    diagnostic, THEN averaged. This is the ordering Prosser uses.
    """
    signs = signs or {}
    exceed_fields = [
        exceedance_field(da, thresholds[name], signs.get(name, "+"))
        for name, da in diagnostics.items()
    ]
    stacked = xr.concat(exceed_fields, dim="diagnostic")
    return stacked.mean(dim="diagnostic").rename("exceedance_mean")


def exceedance_mean_all_severities(
    diagnostics: dict[str, xr.DataArray],
    thresholds: dict[str, dict[str, float]],
    severities: list[str],
    signs: dict[str, Sign] | None = None,
) -> dict[str, xr.DataArray]:
    """Full Q-AGG-1 aggregation: all 21 diagnostics x all severity levels.

    thresholds: {diagnostic_name: {severity_name: threshold_value}}
                e.g. {"ubf": {"light": 1e-9, "moderate": 5e-9, ...}, ...}
    severities: ordered list of severity names to produce, e.g.
                ["light", "light_to_moderate", "moderate",
                 "moderate_to_severe", "severe"]

    Returns {severity_name: exceedance-mean DataArray}.
    """
    out: dict[str, xr.DataArray] = {}
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
