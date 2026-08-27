"""
calib_weighted_percentile.py
=============================
Q-CALIB-5: cos(phi) weighted-percentile function for global calibration
(Q-CALIB-3's inferred scheme: cos(phi) per-cell weight, globally
normalized, weighted empirical CDF). This corrects for Gaussian/regular
lat-lon grid pole-density bias (Q-CALIB-1: rojak's own calibration suite
does NOT do this -- confirmed unweighted, no hook).

Implementation: weighted Hazen-type quantile. Sort (value, weight) pairs
by value; each weight w_i occupies a "slot" of width w_i/W (W = total
weight) in cumulative-weight space; its PLOTTING POSITION is the midpoint
of that slot, exactly generalizing Hazen's (i-0.5)/n unweighted plotting
position to unequal weights:

    P_i = (C_i - 0.5*w_i) / W,   C_i = cumulative weight through i

Target percentile p is located by linear interpolation between the two
bracketing plotting positions (flat extrapolation below P_1 / above P_n,
same edge behavior as standard Hazen). This is a genuine superset: with
all weights equal, C_i - 0.5*w_i = i/n - 0.5/n = (i-0.5)/n, IDENTICAL to
unweighted Hazen -- not an approximation, an exact algebraic reduction.
"""
from __future__ import annotations

import numpy as np


def weighted_percentile(values: np.ndarray, weights: np.ndarray, percentiles: np.ndarray | float) -> np.ndarray:
    """Weighted Hazen-type percentile(s).

    values, weights: 1D arrays, same length. weights must be >= 0, not
                      all zero. NaNs in `values` are dropped (with their
                      matching weight) before computing.
    percentiles:      scalar or array, in [0, 100] (same convention as
                      np.percentile, NOT [0,1]).

    Returns: scalar or array matching `percentiles`' shape.
    """
    values = np.asarray(values, dtype=np.float64).ravel()
    weights = np.asarray(weights, dtype=np.float64).ravel()
    if values.shape != weights.shape:
        raise ValueError(f"values and weights must have the same shape, got {values.shape} vs {weights.shape}")
    finite = np.isfinite(values)
    values, weights = values[finite], weights[finite]
    if np.any(weights < 0):
        raise ValueError("weights must be non-negative")
    total_weight = weights.sum()
    if total_weight <= 0:
        raise ValueError("sum of weights must be positive")

    order = np.argsort(values, kind="stable")
    v_sorted = values[order]
    w_sorted = weights[order]

    # Consolidate exact ties: np.interp requires strictly-increasing x and
    # does not do the mathematically-correct thing with duplicate x-values
    # (it silently picks one match rather than combining their weight mass
    # into a single correct plotting position). Sum weights for identical
    # values before computing cumulative weight / plotting positions.
    v_unique, inverse = np.unique(v_sorted, return_inverse=True)
    w_unique = np.zeros(v_unique.shape, dtype=np.float64)
    np.add.at(w_unique, inverse, w_sorted)

    cumw = np.cumsum(w_unique)
    plotting_positions = (cumw - 0.5 * w_unique) / total_weight  # in [0, 1]

    p = np.asarray(percentiles, dtype=np.float64) / 100.0
    scalar_input = p.ndim == 0
    p = np.atleast_1d(p)

    result = np.interp(p, plotting_positions, v_unique, left=v_unique[0], right=v_unique[-1])
    return float(result[0]) if scalar_input else result
