#!/usr/bin/env python
"""
ada/trend_check.py
==================
The TREND check: nine DJF seasons across the record, fitted, against
Prosser (2023) Table 1's fitted 1979 -> 2020 change.

    pixi run python ada/trend_check.py [--thresholds PATH] [--per-diagnostic]
    pixi run python ada/trend_check.py --years 1979,1990,2000,2010,2020

WHY THIS AND NOT ANOTHER LEVEL CHECK
------------------------------------
The calibration check (ada/calibration_check.py) verified a turbulence
frequency LEVEL in one season of one year, and landed within 24 % of Prosser.
But the paper's claim -- and the reason this project exists -- is the CHANGE:
DJF MOG rising 22.3 -> 30.6 h, +37 %, over 1979-2020. Layers 5 and 6 (annual
aggregation, per-gridpoint regression) have never touched real data.

Two things make this the right test rather than more of the same:

1. A CONSTANT LEVEL OFFSET CANCELS IN A RELATIVE CHANGE. Our levels came in
   at 0.76x Prosser's, and §11.3 of CALIBRATION_REFERENCE attributes that to
   the wider vertical stencil. If that reading is right, the offset divides
   out of a ratio and THE TREND SHOULD REPLICATE BETTER THAN THE LEVEL DID.

2. IT SEPARATES STENCIL FROM NOISE. Prosser's Table 1 numbers are FITTED
   endpoints from a 42-year regression -- his words, "a guide to the
   underlying turbulence statistics in the absence of interannual
   variability". The calibration check compared a single RAW DJF 1979 against
   a FITTED 1979, so some unknown part of that 24 % is simply 1979 being one
   draw. Several seasons show the scatter and let the comparison be
   fitted-against-fitted.

NINE SEASONS, NOT FIVE
----------------------
1979, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020.

The first version of this script used five. That was not enough, and the
failure mode was subtle: with n=5 the MOG fit gave R^2 = 0.59, i.e. t = 2.1
and p ~ 0.13, with a 95 % interval on the +47 % trend running from about
-25 % to +165 %. Prosser's +37 % lay well inside. A test whose interval
contains both the null and the target cannot distinguish them, so its verdict
carried no information whichever way it printed.

At n=9 and the same R^2, t = 3.2 and p ~ 0.016. Twelve extra months of a
North Atlantic box buys the difference between a number and a result.

Hence the second block below now reports a standard error, a t statistic and
a 95 % interval on every fitted trend, and asks the only question that is
actually decidable: DOES PROSSER'S VALUE LIE INSIDE OUR INTERVAL. The bare
ratio is still printed, because it is what a reader compares by eye, but it
is no longer the verdict.

DJF(Y) = Jan(Y) + Feb(Y) + Dec(Y). Prosser's record starts 1 Jan 1979, so he
cannot have used Dec 1978 either; this matches his convention.
"""
from __future__ import annotations

import argparse
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

YEARS = [1979, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020]
ORDER = ["light", "light_to_moderate", "moderate", "moderate_to_severe", "severe"]
LABEL = {"light": "LOG", "light_to_moderate": "LMOG", "moderate": "MOG",
         "moderate_to_severe": "MSOG", "severe": "SOG"}

PROSSER_BOX = dict(lat=(36.0, 60.0), lon=(-55.0, -10.0))

# Prosser (2023) Table 1, DJF row: fitted hours per season at an average point
# in his North Atlantic box.
PROSSER_DJF = {
    #                  1979    2020   rel. increase
    "light":          (128.9, 155.6, 0.21),
    "light_to_moderate": (45.6, 59.3, 0.30),
    "moderate":         (22.3, 30.6, 0.37),
    "moderate_to_severe": (12.1, 17.2, 0.43),
    "severe":            (6.4,  9.6, 0.49),
}

# Two-sided 95 % critical values of Student's t, by degrees of freedom.
# A lookup rather than scipy.stats: this script's only heavy dependency is
# xarray, and one table of twenty numbers is not worth another import.
_TCRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
          7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179,
          13: 2.160, 14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101,
          19: 2.093, 20: 2.086, 25: 2.060, 30: 2.042}


def tcrit(df: int) -> float:
    if df <= 0:
        return float("nan")
    if df in _TCRIT:
        return _TCRIT[df]
    if df > 30:
        return 1.960
    return _TCRIT[min(k for k in _TCRIT if k >= df)]


def djf_days(year: int) -> int:
    """Jan + Feb + Dec of `year`. 91 in a leap year, 90 otherwise."""
    leap = (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
    return 31 + (29 if leap else 28) + 31


def subset_box(ds: xr.Dataset, lat: tuple, lon: tuple) -> xr.Dataset:
    """Subset by bounds regardless of coordinate order.

    ERA5 latitudes descend, so a naive slice(36, 60) returns an EMPTY
    selection with no error -- which would make every number below
    meaningless rather than wrong-looking.
    """
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


def load_djf(base: Path, year: int) -> xr.Dataset:
    paths = [base / "derived/north_atlantic" / f"diagnostics_na_{year}-{m}.zarr"
             for m in ("01", "02", "12")]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"DJF {year}: missing {missing}")
    return xr.concat([xr.open_zarr(p) for p in paths], dim="time")


def ols(x: np.ndarray, y: np.ndarray) -> dict:
    """Straight-line fit with the uncertainty on the slope.

    Returns slope, intercept, R^2, the standard error of the slope, its t
    statistic against zero, and the half-width of a 95 % interval. With n
    under ten the standard error is the whole point -- a slope estimated from
    five noisy seasons has an interval wide enough to contain almost any
    hypothesis, and quoting it without one invites a false comparison.
    """
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--thresholds", default=None,
                    help="thresholds JSON (default: newest in <base>/calibration)")
    ap.add_argument("--per-diagnostic", action="store_true",
                    help="also break the first and last season down by diagnostic")
    ap.add_argument("--years", default=None,
                    help="comma-separated DJF years (default: the nine above). "
                         "Pass 1979,1990,2000,2010,2020 to reproduce the "
                         "five-season run exactly.")
    args = ap.parse_args()

    years = ([int(v) for v in args.years.split(",")] if args.years else list(YEARS))
    if len(years) < 3:
        print("!! need at least three seasons to fit a line with an interval")
        return 1

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
          f"period {provenance['period']}")
    print(f"    severities {provenance['severities_percentiles']}")
    print(f">>> seasons: {len(years)} -- {', '.join(str(y) for y in years)}\n")

    diag = importlib.util.spec_from_file_location(
        "diagnostics", REPO / "2_diagnostics.py")
    dmod = importlib.util.module_from_spec(diag)
    sys.modules["diagnostics"] = dmod
    diag.loader.exec_module(dmod)
    names = [k for k in dmod.REFERENCE_TABLE if k in thresholds]

    # ---------------------------------------------------------------- per year
    print("=" * 78)
    print("PER-SEASON EXCEEDANCE — Prosser box 36-60N / 55-10W")
    print("=" * 78)
    rates: dict[str, list[float]] = {s: [] for s in ORDER}
    per_diag: dict[int, dict[str, float]] = {}

    header = f"   {'DJF':<8}{'days':>6}" + "".join(f"{LABEL[s]:>10}" for s in ORDER)
    print(header)
    for year in years:
        ds = subset_box(load_djf(base, year), **PROSSER_BOX)
        fields = {k: ds[k] for k in names}
        w = lat_weights_for(ds)
        per_sev = aggregate.exceedance_mean_all_severities(
            fields, thresholds, ORDER, signs)
        row = []
        for s in ORDER:
            r = weighted_rate(per_sev[s]["exceedance_mean"], w)
            rates[s].append(r)
            row.append(r)
        print(f"   {year:<8}{djf_days(year):>6}" + "".join(f"{v:>9.3%}" for v in row))

        if args.per_diagnostic and year in (years[0], years[-1]):
            per_diag[year] = {
                n: weighted_rate(
                    aggregate.exceedance_field(
                        fields[n], thresholds[n]["moderate"], signs[n]), w)
                for n in names
            }
        del ds, fields

    # ------------------------------------------------------------------- trend
    print("\n" + "=" * 78)
    print("FITTED TREND vs PROSSER (2023) TABLE 1, DJF")
    print("=" * 78)
    print(f"   Prosser fits 42 years; this fits {len(years)}. The 95 % interval")
    print("   is on the slope only, with the fitted 1979 level held at its")
    print("   point estimate -- so it understates total uncertainty slightly.\n")
    print(f"   {'':<6}{'our 1979':>10}{'our 2020':>10}{'ours %':>9}"
          f"{'95% CI':>17}{'Prosser %':>11}{'ratio':>8}{'R^2':>7}{'t':>7}")

    x = np.array(years, float)
    summary = {}
    for s in ORDER:
        y = np.array(rates[s])
        f = ols(x, y)
        fit79 = f["slope"] * 1979 + f["intercept"]
        fit20 = f["slope"] * 2020 + f["intercept"]
        rel = (fit20 - fit79) / fit79 if fit79 > 0 else np.nan
        # translate the slope interval into the same relative-change units
        rel_half = (41.0 * f["half"] / fit79) if fit79 > 0 else np.nan
        lo, hi = rel - rel_half, rel + rel_half
        p_rel = PROSSER_DJF[s][2]
        summary[s] = dict(fit79=fit79, fit20=fit20, rel=rel, lo=lo, hi=hi,
                          p_rel=p_rel, **f)
        print(f"   {LABEL[s]:<6}{fit79:>9.3%}{fit20:>10.3%}{rel:>8.0%}"
              f"   [{lo:>5.0%},{hi:>5.0%}]"
              f"{p_rel:>10.0%}{rel / p_rel:>8.2f}{f['r2']:>7.2f}{f['t']:>7.2f}")

    print(f"\n   {'':<6}{'ours (h)':>11}{'':>10}{'Prosser (h)':>13}")
    for s in ORDER:
        d = summary[s]
        print(f"   {LABEL[s]:<6}{d['fit79'] * djf_days(1979) * 24:>11.1f}"
              f" ->{d['fit20'] * djf_days(2020) * 24:>7.1f}"
              f"{PROSSER_DJF[s][0]:>10.1f} ->{PROSSER_DJF[s][1]:>6.1f}")

    # ------------------------------------------------------------------ verdict
    mog = summary["moderate"]
    mog_rel, mog_p = mog["rel"], mog["p_rel"]
    signs_ok = all(summary[s]["rel"] > 0 for s in ORDER)
    monotone = all(summary[ORDER[i]]["rel"] <= summary[ORDER[i + 1]]["rel"] + 0.05
                   for i in range(len(ORDER) - 1))
    sig = [s for s in ORDER if abs(summary[s]["t"]) >= tcrit(summary[s]["df"])]
    inside = [s for s in ORDER
              if summary[s]["lo"] <= summary[s]["p_rel"] <= summary[s]["hi"]]

    print("\n" + "-" * 78)
    print(f"   1. ALL TRENDS POSITIVE       {'PASS' if signs_ok else 'FAIL'}")
    print(f"   2. TRENDS DISTINGUISHABLE FROM ZERO AT 5 %   "
          f"{len(sig)}/5   {'PASS' if len(sig) >= 3 else 'WEAK'}")
    print(f"      (with n={len(years)}, df={mog['df']}, "
          f"|t| must exceed {tcrit(mog['df']):.2f})")
    print(f"   3. PROSSER INSIDE OUR 95 % INTERVAL          "
          f"{len(inside)}/5   {'PASS' if len(inside) >= 3 else 'CHECK'}")
    print(f"   4. STRONGER AT HIGHER SEVERITY   {'PASS' if monotone else 'CHECK'}"
          "   (Prosser: 21/30/37/43/49 %)")

    # ------------------------------------------- level deficit vs trend excess
    # The five-season run found these two mirror each other across severity.
    # If that holds it is not two problems but one: a vertical stencil too
    # wide by ~3x puts us lower on the tail of the EDR distribution, and lower
    # on the tail means the SAME physical intensification produces a LARGER
    # relative change, because exceedance probability is convex in the shift.
    # A constant multiplicative offset would instead leave the trend ratio at
    # 1.00 for every severity. Printing both columns makes the difference
    # between those two stories visible rather than asserted.
    print("\n" + "-" * 78)
    print("   LEVEL DEFICIT vs TREND EXCESS — are they mirror images?")
    print(f"   {'':<6}{'level ratio':>13}{'trend ratio':>13}{'product':>10}")
    for s in ORDER:
        d = summary[s]
        lvl = d["fit79"] / (PROSSER_DJF[s][0] / (djf_days(1979) * 24))
        trd = d["rel"] / d["p_rel"]
        print(f"   {LABEL[s]:<6}{lvl:>13.2f}{trd:>13.2f}{lvl * trd:>10.2f}")
    print("   Level falling and trend rising together across severity is the")
    print("   tail-position signature of a too-wide stencil. Both columns flat")
    print("   at 1.00 would mean a clean multiplicative offset instead; a")
    print("   product drifting far from 1.00 means neither story fits.")

    level_ratio = mog["fit79"] / (PROSSER_DJF["moderate"][0]
                                  / (djf_days(1979) * 24))
    print("\n" + "-" * 78)
    print(f"   MOG level ratio at fitted 1979: {level_ratio:.2f}"
          f"   trend ratio: {mog_rel / mog_p:.2f}")
    print(f"   trend {abs(mog_rel / mog_p - 1):.0%} off vs level "
          f"{abs(level_ratio - 1):.0%} off")
    print("   Read this against criterion 3, not on its own: with a wide")
    print("   interval a closer ratio is not evidence and a wider one is not")
    print("   refutation.")
    print("-" * 78)

    # ---------------------------------------------------------- per diagnostic
    if per_diag:
        print("\n" + "=" * 78)
        print("PER-DIAGNOSTIC MOG EXCEEDANCE — what the ensemble mean hides")
        print("=" * 78)
        print("   The 21-member mean can look right while one member is far off;")
        print("   a badly wrong diagnostic moves it by only 1/21. Prosser's")
        print("   Figure S5 is the published version of this breakdown.\n")
        a, b = years[0], years[-1]
        print(f"   {'diagnostic':<23}{f'DJF {a}':>11}{f'DJF {b}':>11}{'change':>9}")
        rows = sorted(per_diag[a], key=lambda n: per_diag[a][n], reverse=True)
        for n in rows:
            r0, r1 = per_diag[a][n], per_diag[b][n]
            ch = (r1 - r0) / r0 if r0 > 0 else np.nan
            print(f"   {n:<23}{r0:>10.3%}{r1:>11.3%}{ch:>9.0%}")
        vals = np.array([per_diag[a][n] for n in rows])
        print(f"\n   spread across the 21 at DJF {a}: "
              f"min {vals.min():.3%}, median {np.median(vals):.3%}, "
              f"max {vals.max():.3%}")
        print(f"   ratio max/min = {vals.max() / max(vals.min(), 1e-12):,.0f}x")
        print("   A diagnostic flagging orders more or fewer cells than its")
        print("   siblings is the signal to chase, not the ensemble mean.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
