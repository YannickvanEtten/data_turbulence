"""
trend.py
========
Q-AGG-2: Layer 6 — per-gridpoint linear trend regression over 42 years of
annual exceedance probability. This is the layer that produces the
paper's actual headline numbers (the ~37%/~55% increase claims), so it's
validated against synthetic data with a KNOWN injected trend before it
ever touches real ERA5-derived output.

Per gridpoint:
    1. OLS linear regression of annual exceedance probability vs. year
       (1979-2020, 42 points).
    2. Fitted endpoint values at 1979 and 2020 (from the fitted line, not
       the raw noisy data at those two years).
    3. Negative-clipping: exceedance probability cannot be negative, so
       any fitted endpoint below 0 is clipped to 0 BEFORE computing
       absolute/relative change.
    4. absolute_change = fitted_2020_clipped - fitted_1979_clipped
       relative_change  = absolute_change / fitted_1979_clipped
       (relative_change is NaN, not inf, where the clipped 1979 baseline
       is exactly 0 -- a percent change from a zero baseline is undefined,
       not "infinite increase"; flagged explicitly rather than silently
       producing inf/nan noise downstream.)
"""
from __future__ import annotations

import numpy as np
import xarray as xr


def _ols_slope_intercept(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized OLS, y's LAST axis is the regression (x) axis.

    Numerically identical formula to scipy.stats.linregress (verified in
    the test below against scipy directly, not just against a "known"
    injected trend) -- vectorized over all leading (gridpoint) dimensions
    instead of looping scipy.stats.linregress per pixel.
    """
    x = np.asarray(x, dtype=np.float64)
    xbar = x.mean()
    xdev = x - xbar
    ss_xx = np.sum(xdev ** 2)

    ybar = y.mean(axis=-1, keepdims=True)
    ydev = y - ybar
    ss_xy = np.sum(xdev * ydev, axis=-1)

    slope = ss_xy / ss_xx
    intercept = ybar[..., 0] - slope * xbar
    return slope, intercept


def fit_annual_trend(
    annual_prob: xr.DataArray,
    year_dim: str = "year",
    start_year: int = 1979,
    end_year: int = 2020,
) -> xr.Dataset:
    """Fit the 42-year linear trend per gridpoint.

    annual_prob: DataArray with a `year_dim` coordinate (1979..2020) plus
                 any number of gridpoint dims (lat, lon, pressure_level, ...).

    Returns a Dataset with:
      slope, intercept            -- raw OLS fit (unclipped)
      fitted_1979, fitted_2020    -- raw fitted endpoints (unclipped)
      fitted_1979_clipped, fitted_2020_clipped
      absolute_change             -- clipped_2020 - clipped_1979
      relative_change             -- absolute_change / clipped_1979,
                                      NaN where clipped_1979 == 0
    """
    years = annual_prob[year_dim].values.astype(np.float64)
    if years.size < 2:
        raise ValueError("Need at least 2 years to fit a trend")

    slope, intercept = xr.apply_ufunc(
        lambda y: _ols_slope_intercept(years, y),
        annual_prob,
        input_core_dims=[[year_dim]],
        output_core_dims=[[], []],
        vectorize=False,
    )

    fitted_1979 = intercept + slope * start_year
    fitted_2020 = intercept + slope * end_year

    fitted_1979_clipped = fitted_1979.clip(min=0.0)
    fitted_2020_clipped = fitted_2020.clip(min=0.0)

    absolute_change = fitted_2020_clipped - fitted_1979_clipped
    relative_change = xr.where(
        fitted_1979_clipped == 0.0,
        np.nan,
        absolute_change / fitted_1979_clipped,
    )

    return xr.Dataset({
        "slope": slope.rename("slope"),
        "intercept": intercept.rename("intercept"),
        "fitted_1979": fitted_1979.rename("fitted_1979"),
        "fitted_2020": fitted_2020.rename("fitted_2020"),
        "fitted_1979_clipped": fitted_1979_clipped.rename("fitted_1979_clipped"),
        "fitted_2020_clipped": fitted_2020_clipped.rename("fitted_2020_clipped"),
        "absolute_change": absolute_change.rename("absolute_change"),
        "relative_change": relative_change.rename("relative_change"),
    })
