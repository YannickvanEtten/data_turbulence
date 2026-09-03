#!/usr/bin/env python
"""
ada/full_trend_check.py
========================
The DEFINITIVE trend check: the complete 1979-2020 series, all four seasons
plus annual, fitted per severity, against the FULL Prosser (2023) Table 1 --
not the nine-DJF-season proxy in ada/trend_check.py.

    pixi run python ada/full_trend_check.py
    pixi run python ada/full_trend_check.py --season djf
    pixi run python ada/full_trend_check.py --season all --per-diagnostic

RUN ada/check_production_complete.py FIRST. This script trusts that the 504
North Atlantic months and the calibration thresholds are complete and built
under matching conventions; it does not re-check that itself.

WHY A NEW SCRIPT, NOT AN EDIT TO ada/trend_check.py
----------------------------------------------------
ada/trend_check.py's nine-season DJF result is already recorded in
STATUS.md §12, with `--years` specifically documented to reproduce it exactly.
Editing that script to do something bigger risks silently perturbing numbers
that are already cited. This script is additive: same helper shapes
(box-subsetting, weighted rate, OLS with an interval), same season convention,
extended to every season Prosser publishes and to the full 42 years now that
they exist. ada/trend_check.py is unchanged and still reproduces its own
result.

THE ONE THING THAT CHANGES BY GOING FROM n=9 TO n=42
------------------------------------------------------
STATUS.md §12.3-12.4 found the trend check informative but underpowered: MOG's
t-statistic was 2.20 against a threshold of 2.37 at n=9, and interannual
variability (a factor of 1.65 across nine winters) was comparable in size to
the entire 41-year fitted change. Both of those are sample-size artifacts in
part. At n=42 (df=40) the 5% critical t is ~2.02, comfortably below values
that were previously borderline -- so this run is what actually answers
STATUS.md §12.3's open question, not a stronger version of the same proxy.

SEASON DEFINITION -- MATCHES ada/trend_check.py EXACTLY
----------------------------------------------------------
DJF(Y) = Jan(Y) + Feb(Y) + Dec(Y), all from the SAME calendar year Y (Prosser's
record starts 1 Jan 1979, so Dec(1978) was never an option; ada/trend_check.py
established this convention and it is preserved here). MAM/JJA/SON only ever
draw from one calendar year, so they need no special handling. Annual is all
twelve months of Y. Because every season for year Y is a subset of Y's own
twelve months, this script opens each year's twelve NA zarr stores ONCE and
derives all five season groups from that one lazy dataset by masking on
`time.dt.month`, rather than re-opening files per season.

MEMORY
------
xr.open_zarr is lazy; nothing is materialized until a reduction runs. Each
year's twelve stores, box-subset to Prosser's 36-60N/55-10W (a fraction of the
full 30-60N/75-0W NA grid), are processed and discarded before the next year
loads -- so peak memory is bounded by one year's box-subset, not by the whole
504-month series. If this still runs out of memory on the smallest defq nodes,
pass --season djf/mam/jja/son one at a time instead of the default `all`,
which is the same total work spread across separate, smaller jobs.
"""
from __future__ import annotations

import argparse
import calendar
import importlib.util
import sys
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
sys.path.insert(0, str(REPO))

import aggregate                                    # noqa: E402
import calibration                                  # noqa: E402

ORDER = ["light", "light_to_moderate", "moderate", "moderate_to_severe", "severe"]
LABEL = {"light": "LOG", "light_to_moderate": "LMOG", "moderate": "MOG",
         "moderate_to_severe": "MSOG", "severe": "SOG"}

PROSSER_BOX = dict(lat=(36.0, 60.0), lon=(-55.0, -10.0))

SEASONS = ["djf", "mam", "jja", "son", "annual"]
SEASON_MONTHS = {
    "djf": (12, 1, 2),
    "mam": (3, 4, 5),
    "jja": (6, 7, 8),
    "son": (9, 10, 11),
    "annual": tuple(range(1, 13)),
}
SEASON_TITLE = {"djf": "DJF", "mam": "MAM", "jja": "JJA", "son": "SON",
                "annual": "Annual"}

# Prosser (2023) Table 1, in full. Hours per season/year at an average point
# in the North Atlantic box, 1979 -> 2020, and the paper's own rounded
# relative-increase figure. CALIBRATION_REFERENCE.md §4.4.
PROSSER_TABLE = {
    "djf": {
        "light":              (128.9, 155.6, 0.21),
        "light_to_moderate":  (45.6,  59.3,  0.30),
        "moderate":           (22.3,  30.6,  0.37),
        "moderate_to_severe": (12.1,  17.2,  0.43),
        "severe":             (6.4,   9.6,   0.49),
    },
    "mam": {
        "light":              (90.4,  113.4, 0.26),
        "light_to_moderate":  (27.2,  38.9,  0.43),
        "moderate":           (11.8,  18.6,  0.57),
        "moderate_to_severe": (5.7,   9.7,   0.71),
        "severe":             (2.7,   5.0,   0.85),
    },
    "jja": {
        "light":              (114.1, 124.5, 0.09),
        "light_to_moderate":  (36.5,  43.8,  0.20),
        "moderate":           (16.1,  21.1,  0.31),
        "moderate_to_severe": (7.7,   10.9,  0.41),
        "severe":             (3.6,   5.5,   0.52),
    },
    "son": {
        "light":              (133.1, 153.2, 0.15),
        "light_to_moderate":  (43.4,  53.4,  0.23),
        "moderate":           (19.8,  25.8,  0.31),
        "moderate_to_severe": (10.0,  13.9,  0.39),
        "severe":             (5.0,   7.4,   0.47),
    },
    "annual": {
        "light":              (466.5, 546.8, 0.17),
        "light_to_moderate":  (152.7, 195.4, 0.28),
        "moderate":           (70.0,  96.1,  0.37),
        "moderate_to_severe": (35.5,  51.8,  0.46),
        "severe":             (17.7,  27.4,  0.55),
    },
}

# Two-sided 95% critical values of Student's t, by degrees of freedom. n=42
# gives df=40; kept as a lookup for the same reason ada/trend_check.py gives --
# this script's only heavy dependency stays xarray, not scipy.
_TCRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 10: 2.228,
          15: 2.131, 20: 2.086, 25: 2.060, 30: 2.042, 40: 2.021, 60: 2.000,
          120: 1.980}


def tcrit(df: int) -> float:
    if df <= 0:
        return float("nan")
    if df in _TCRIT:
        return _TCRIT[df]
    if df > 120:
        return 1.960
    return _TCRIT[min(k for k in _TCRIT if k >= df)]


def season_days(year: int, season: str) -> int:
    leap = calendar.isleap(year)
    if season == "djf":
        return 31 + (29 if leap else 28) + 31
    if season == "mam":
        return 31 + 30 + 31
    if season == "jja":
        return 30 + 31 + 31
    if season == "son":
        return 30 + 31 + 30
    if season == "annual":
        return 366 if leap else 365
    raise ValueError(season)


def subset_box(ds: xr.Dataset, lat: tuple, lon: tuple) -> xr.Dataset:
    """Subset by bounds regardless of coordinate order (ERA5 latitudes descend)."""
    la, lo = ds["latitude"].values, ds["longitude"].values
    lat_sl = slice(lat[1], lat[0]) if la[0] > la[-1] else slice(lat[0], lat[1])
    lon_sl = slice(lon[1], lon[0]) if lo[0] > lo[-1] else slice(lon[0], lon[1])
    out = ds.sel(latitude=lat_sl, longitude=lon_sl)
    if out.sizes["latitude"] == 0 or out.sizes["longitude"] == 0:
        raise ValueError(f"box {lat} x {lon} selected nothing")
    return out


def lat_weights_for(ds: xr.Dataset) -> xr.DataArray:
    return xr.DataArray(np.cos(np.deg2rad(ds["latitude"].values)),
                        coords={"latitude": ds["latitude"]}, dims=("latitude",))


def weighted_rate(exceed: xr.DataArray, weights: xr.DataArray) -> float:
    w = weights.broadcast_like(exceed)
    populated = exceed.notnull()
    return float((exceed.fillna(0.0) * w).sum() / w.where(populated).sum())


def load_year(base: Path, year: int) -> xr.Dataset:
    """Open all twelve North Atlantic months of `year`, lazily, concatenated."""
    paths = [base / "derived/north_atlantic" / f"diagnostics_na_{year}-{m:02d}.zarr"
             for m in range(1, 13)]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"{year}: missing {missing}")
    return xr.concat([xr.open_zarr(p) for p in paths], dim="time")


def ols(x: np.ndarray, y: np.ndarray) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")

    df = n - 2
    sxx = float(((x - x.mean()) ** 2).sum())
    if df > 0 and sxx > 0:
        se = float(np.sqrt((ss_res / df) / sxx))
    else:
        se = float("nan")
    t = slope / se if se and np.isfinite(se) and se > 0 else float("inf")
    half = tcrit(df) * se if np.isfinite(se) else float("nan")
    return dict(slope=float(slope), intercept=float(intercept), r2=r2,
                se=se, t=float(t), half=float(half), n=n, df=df)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--thresholds", default=None,
                    help="thresholds JSON (default: newest in <base>/calibration)")
    ap.add_argument("--season", default="all",
                    choices=SEASONS + ["all"],
                    help="which season(s) to fit. 'all' reproduces the full "
                         "Prosser Table 1 (five rows); a single season is "
                         "cheaper and useful for splitting the job.")
    ap.add_argument("--start-year", type=int, default=1979)
    ap.add_argument("--end-year", type=int, default=2020)
    ap.add_argument("--year-step", type=int, default=1,
                    help="use every Nth year (default 1 = all 42). Only for "
                         "a fast provisional pass -- the whole point of this "
                         "script is to use the full n=42, so leave at 1 for "
                         "anything you intend to report.")
    ap.add_argument("--per-diagnostic", action="store_true",
                    help="also break the first and last year down by "
                         "diagnostic, for the first season requested")
    args = ap.parse_args()

    years = list(range(args.start_year, args.end_year + 1, args.year_step))
    if len(years) < 3:
        print("!! need at least three years to fit a line with an interval")
        return 1
    seasons = SEASONS if args.season == "all" else [args.season]

    base = Path(args.base)
    if args.thresholds:
        tpath = Path(args.thresholds)
    else:
        found = sorted((base / "calibration").glob("thresholds_*.json"))
        if not found:
            print("!! no thresholds file; run jobs/06_calibration_check.sbatch first")
            return 1
        tpath = found[-1]

    thresholds, signs, provenance = calibration.load_thresholds(tpath)
    print(f">>> thresholds: {tpath.name}")
    print(f"    domain {provenance['calibration_domain']}, "
          f"period {provenance['period']}\n")
    print(f">>> years: {len(years)} ({years[0]}-{years[-1]}, step {args.year_step})")
    print(f">>> seasons: {', '.join(SEASON_TITLE[s] for s in seasons)}\n")
    print("!! Run ada/check_production_complete.py first if you have not --")
    print("   this script does not itself verify every store's provenance.\n")

    diag = importlib.util.spec_from_file_location(
        "diagnostics", REPO / "2_diagnostics.py")
    dmod = importlib.util.module_from_spec(diag)
    sys.modules["diagnostics"] = dmod
    diag.loader.exec_module(dmod)
    names = [k for k in dmod.REFERENCE_TABLE if k in thresholds]

    # rates[season][severity] -> list of per-year weighted exceedance rates
    rates: dict[str, dict[str, list[float]]] = {
        s: {sev: [] for sev in ORDER} for s in seasons
    }
    per_diag: dict[str, dict[int, dict[str, float]]] = {s: {} for s in seasons}
    first_season = seasons[0]

    print("=" * 78)
    print("LOADING — one pass per year, all requested seasons sliced from it")
    print("=" * 78)
    for year in years:
        ds_year = load_year(base, year)
        ds_box = subset_box(ds_year, **PROSSER_BOX)
        w = lat_weights_for(ds_box)

        for season in seasons:
            months = SEASON_MONTHS[season]
            mask = ds_box["time"].dt.month.isin(months)
            ds_s = ds_box.isel(time=mask.values)

            expected = season_days(year, season) * 8
            actual = ds_s.sizes["time"]
            if actual != expected:
                print(f"   !! {year} {SEASON_TITLE[season]}: {actual} timesteps, "
                      f"expected {expected} — using actual count, treat this "
                      f"year/season with caution")

            fields = {k: ds_s[k] for k in names}
            per_sev = aggregate.exceedance_mean_all_severities(
                fields, thresholds, ORDER, signs)
            for sev in ORDER:
                rates[season][sev].append(weighted_rate(
                    per_sev[sev]["exceedance_mean"], w))

            if (args.per_diagnostic and season == first_season
                    and year in (years[0], years[-1])):
                per_diag[season][year] = {
                    n: weighted_rate(
                        aggregate.exceedance_field(
                            fields[n], thresholds[n]["moderate"], signs[n]), w)
                    for n in names
                }
            del ds_s, fields

        print(f"   {year}  done")
        del ds_year, ds_box

    # ------------------------------------------------------------ per season
    all_summaries: dict[str, dict[str, dict]] = {}
    for season in seasons:
        print("\n" + "=" * 78)
        print(f"{SEASON_TITLE[season]} — FITTED TREND vs PROSSER (2023) TABLE 1")
        print("=" * 78)
        x = np.array(years, float)
        summary = {}
        print(f"   {'':<6}{'our 1979':>10}{'our 2020':>10}{'ours %':>9}"
              f"{'95% CI':>17}{'Prosser %':>11}{'ratio':>8}{'R^2':>7}{'t':>7}")
        for sev in ORDER:
            y = np.array(rates[season][sev])
            f = ols(x, y)
            fit_lo = f["slope"] * args.start_year + f["intercept"]
            fit_hi = f["slope"] * args.end_year + f["intercept"]
            rel = (fit_hi - fit_lo) / fit_lo if fit_lo > 0 else np.nan
            span = float(args.end_year - args.start_year)
            rel_half = (span * f["half"] / fit_lo) if fit_lo > 0 else np.nan
            lo, hi = rel - rel_half, rel + rel_half
            p_rel = PROSSER_TABLE[season][sev][2]
            summary[sev] = dict(fit_lo=fit_lo, fit_hi=fit_hi, rel=rel, lo=lo,
                                hi=hi, p_rel=p_rel, **f)
            print(f"   {LABEL[sev]:<6}{fit_lo:>9.3%}{fit_hi:>10.3%}{rel:>8.0%}"
                  f"   [{lo:>5.0%},{hi:>5.0%}]"
                  f"{p_rel:>10.0%}{rel / p_rel if p_rel else float('nan'):>8.2f}"
                  f"{f['r2']:>7.2f}{f['t']:>7.2f}")

        print(f"\n   {'':<6}{'ours (h)':>11}{'':>10}{'Prosser (h)':>13}")
        for sev in ORDER:
            d = summary[sev]
            h_lo = d["fit_lo"] * season_days(args.start_year, season) * 24
            h_hi = d["fit_hi"] * season_days(args.end_year, season) * 24
            p_lo, p_hi, _ = PROSSER_TABLE[season][sev]
            print(f"   {LABEL[sev]:<6}{h_lo:>11.1f}"
                  f" ->{h_hi:>7.1f}{p_lo:>10.1f} ->{p_hi:>6.1f}")

        signs_ok = all(summary[sev]["rel"] > 0 for sev in ORDER)
        sig = [sev for sev in ORDER
               if abs(summary[sev]["t"]) >= tcrit(summary[sev]["df"])]
        inside = [sev for sev in ORDER
                  if summary[sev]["lo"] <= summary[sev]["p_rel"] <= summary[sev]["hi"]]
        df0 = summary["moderate"]["df"]
        print("\n" + "-" * 78)
        print(f"   1. ALL TRENDS POSITIVE                        "
              f"{'PASS' if signs_ok else 'FAIL'}")
        print(f"   2. TRENDS DISTINGUISHABLE FROM ZERO AT 5%      "
              f"{len(sig)}/5   {'PASS' if len(sig) >= 3 else 'WEAK'}"
              f"   (n={len(years)}, df={df0}, |t| must exceed {tcrit(df0):.2f})")
        print(f"   3. PROSSER INSIDE OUR 95% INTERVAL             "
              f"{len(inside)}/5   {'PASS' if len(inside) >= 3 else 'CHECK'}")
        print("-" * 78)

        all_summaries[season] = summary

    # -------------------------------------------------------- master summary
    print("\n" + "=" * 78)
    print("MASTER SUMMARY — all seasons x severities")
    print("=" * 78)
    print(f"   {'season':<8}{'severity':<7}{'ours %':>8}{'Prosser %':>11}"
          f"{'ratio':>8}{'t':>7}{'sig?':>6}{'Prosser in CI?':>16}")
    n_sig = n_inside = n_total = 0
    for season in seasons:
        for sev in ORDER:
            d = all_summaries[season][sev]
            n_total += 1
            is_sig = abs(d["t"]) >= tcrit(d["df"])
            is_inside = d["lo"] <= d["p_rel"] <= d["hi"]
            n_sig += is_sig
            n_inside += is_inside
            print(f"   {SEASON_TITLE[season]:<8}{LABEL[sev]:<7}{d['rel']:>7.0%}"
                  f"{d['p_rel']:>10.0%}"
                  f"{d['rel'] / d['p_rel'] if d['p_rel'] else float('nan'):>8.2f}"
                  f"{d['t']:>7.2f}{'yes' if is_sig else 'no':>6}"
                  f"{'yes' if is_inside else 'no':>16}")
    print(f"\n   significant at 5%: {n_sig}/{n_total}    "
          f"Prosser inside 95% CI: {n_inside}/{n_total}")
    print("   This is the number STATUS.md §12.3 was waiting on n=42 to")
    print("   settle. Compare against the n=9 DJF-only figures already on")
    print("   record (2/5 significant, 5/5 inside-CI) before writing anything")
    print("   into STATUS.md or CAT_overview.tex.")

    # ---------------------------------------------------------- per diagnostic
    if args.per_diagnostic and per_diag.get(first_season):
        print("\n" + "=" * 78)
        print(f"PER-DIAGNOSTIC MOG EXCEEDANCE — {SEASON_TITLE[first_season]}, "
              f"{years[0]} vs {years[-1]}")
        print("=" * 78)
        a, b = years[0], years[-1]
        d0, d1 = per_diag[first_season][a], per_diag[first_season][b]
        print(f"   {'diagnostic':<23}{f'{a}':>11}{f'{b}':>11}{'change':>9}")
        rows = sorted(d0, key=lambda n: d0[n], reverse=True)
        for n in rows:
            r0, r1 = d0[n], d1[n]
            ch = (r1 - r0) / r0 if r0 > 0 else np.nan
            print(f"   {n:<23}{r0:>10.3%}{r1:>11.3%}{ch:>9.0%}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
