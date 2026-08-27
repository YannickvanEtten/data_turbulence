"""
annual_aggregate.py
====================
Q-AGG-5: the year-grouping bridge between aggregate.py (per-timestep
exceedance-mean fields) and trend.py (expects an already-annual
`annual_prob` series). This layer didn't exist anywhere in the codebase.

Exact spec, Prosser (2023) Methods: "For each year, an average exceedance
field was calculated by taking the mean of the 21 exceedance fields...
Exceedances were converted into percentage probabilities of exceedance,
by normalizing by the number of three-hour periods in each year."

aggregate.py's exceedance_mean_single_severity()/exceedance_mean_all_severities()
already does the "mean of the 21 exceedance fields" part, per timestep.
This module does the remaining piece: group those per-timestep fields by
calendar year, sum over the year, and normalize by the number of 3-hour
periods IN THAT SPECIFIC YEAR -- 2920 for a non-leap year (365*8), 2928
for a leap year (366*8) at 3-hourly resolution. Prosser's paper states
2920 as an illustrative example (the non-leap count) -- summing then
dividing by a fixed 2920 for every year would silently under-normalize
leap years (2928 real periods / 2920 assumed => ~0.27% inflation, small
per-year but systematic and directional across the 11 leap years in
1979-2020, landing on exactly the trend-signal quantity).
"""
from __future__ import annotations

import calendar

import numpy as np
import xarray as xr


def periods_in_year(year: int, hours_per_period: int = 3) -> int:
    """Number of `hours_per_period`-hour periods in a calendar year,
    leap-year aware. E.g. periods_in_year(2000, 3) == 2928 (leap),
    periods_in_year(2001, 3) == 2920 (non-leap)."""
    if 24 % hours_per_period != 0:
        raise ValueError(f"hours_per_period must divide 24 evenly, got {hours_per_period}")
    days = 366 if calendar.isleap(year) else 365
    return days * (24 // hours_per_period)


def annual_exceedance_probability(
    per_timestep_exceedance: xr.DataArray,
    time_dim: str = "time",
    year_dim: str = "year",
    hours_per_period: int = 3,
    validate_complete: bool = True,
) -> xr.DataArray:
    """Group a per-timestep exceedance-mean field (aggregate.py's output)
    by calendar year, sum within each year, and normalize by that year's
    ACTUAL (leap-aware) number of 3-hour periods -- not a fixed constant.

    per_timestep_exceedance: DataArray with a datetime `time_dim` coord
                              (values in [0, 1], aggregate.py's output)
                              plus any gridpoint dims (lat, lon, ...).
    validate_complete: if True (default), raises if a year's ACTUAL
                        timestep count doesn't match the expected
                        leap-aware count -- refuses to silently normalize
                        an incomplete year by the wrong divisor. Set False
                        for legitimate partial-year data (e.g. a trial
                        run), in which case the ACTUAL count is used as
                        the divisor instead of the calendar-expected one.

    Returns a DataArray with a `year_dim` integer coordinate (one value
    per calendar year present) and the same gridpoint dims as the input --
    exactly what trend.py's fit_annual_trend() expects.
    """
    years_coord = per_timestep_exceedance[time_dim].dt.year
    unique_years = sorted(set(int(y) for y in years_coord.values))

    annual_slices = []
    for year in unique_years:
        mask = years_coord == year
        year_data = per_timestep_exceedance.isel({time_dim: mask.values})
        actual_count = year_data.sizes[time_dim]
        expected_count = periods_in_year(year, hours_per_period)

        if validate_complete and actual_count != expected_count:
            kind = "leap" if calendar.isleap(year) else "non-leap"
            raise ValueError(
                f"Year {year} ({kind}): found {actual_count} timesteps in the data, "
                f"expected {expected_count} ({hours_per_period}-hourly). Refusing to "
                f"silently normalize an incomplete year by the wrong divisor. Pass "
                f"validate_complete=False to use the actual count instead (e.g. for "
                f"a legitimate partial-year/trial run)."
            )
        divisor = expected_count if validate_complete else actual_count

        annual_pct = (year_data.sum(dim=time_dim) / divisor) * 100.0
        annual_slices.append(annual_pct)

    out = xr.concat(annual_slices, dim=year_dim)
    out = out.assign_coords({year_dim: np.array(unique_years, dtype=np.int64)})
    return out.rename("annual_exceedance_probability")
