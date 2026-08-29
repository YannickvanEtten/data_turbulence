#!/usr/bin/env python
"""
ada/trend_check.py
==================
The TREND check: five DJF seasons across the record, fitted, against
Prosser (2023) Table 1's fitted 1979 -> 2020 change.

    pixi run python ada/trend_check.py [--thresholds PATH] [--per-diagnostic]

WHY THIS AND NOT ANOTHER LEVEL CHECK
------------------------------------
The calibration check (ada/calibration_check.py) verified a turbulence
frequency LEVEL in one season of one year, and landed within 24 % of Prosser.
But the paper's claim -- and the reason this project exists -- is the CHANGE:
DJF MOG rising 22.3 -> 30.6 h, +37 %, over 1979-2020. Layers 5 and 6 (annual
aggregation, per-gridpoint regression) have never touched real data.

Two things make this the right next test rather than more of the same:

1. A CONSTANT LEVEL OFFSET CANCELS IN A RELATIVE CHANGE. Our levels came in
   at 0.76x Prosser's, and §11.3 of CALIBRATION_REFERENCE attributes that to
   the wider vertical stencil. If that reading is right, the offset divides
   out of a ratio and THE TREND SHOULD REPLICATE BETTER THAN THE LEVEL DID.
   That is a falsifiable prediction, not a hope.

2. IT SEPARATES STENCIL FROM NOISE. Prosser's Table 1 numbers are FITTED
   endpoints from a 42-year regression -- his words, "a guide to the
   underlying turbulence statistics in the absence of interannual
   variability". The calibration check compared a single RAW DJF 1979 against
   a FITTED 1979, so some unknown part of that 24 % is simply 1979 being one
   draw. Five seasons show the scatter and let the comparison be
   fitted-against-fitted.

FIVE SEASONS, NOT TWO
---------------------
1979, 1990, 2000, 2010, 2020. Two points cannot distinguish a trend from two
noisy draws; five can, and they are spread evenly enough to be a fair sample
of a 42-year line. This is a small-sample estimate of a slope Prosser fitted
from 42 -- so read the SIGN and rough MAGNITUDE, not the third digit.

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

YEARS = [1979, 1990, 2000, 2010, 2020]
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


def ols(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Slope, intercept, R^2 of a straight line through five points."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return float(slope), float(intercept), (1 - ss_res / ss_tot if ss_tot else np.nan)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--thresholds", default=None,
                    help="thresholds JSON (default: newest in <base>/calibration)")
    ap.add_argument("--per-diagnostic", action="store_true",
                    help="also break the first and last season down by diagnostic")
    args = ap.parse_args()

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
    print(f"    severities {provenance['severities_percentiles']}\n")

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
    for year in YEARS:
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

        if args.per_diagnostic and year in (YEARS[0], YEARS[-1]):
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
    print("   Prosser fits 42 years; this fits 5. Read sign and magnitude,")
    print("   not the third digit.\n")
    print(f"   {'':<6}{'our 1979':>10}{'our 2020':>10}{'ours %':>9}"
          f"{'Prosser %':>11}{'ratio':>8}{'R^2':>7}")

    x = np.array(YEARS, float)
    summary = {}
    for s in ORDER:
        y = np.array(rates[s])
        slope, intercept, r2 = ols(x, y)
        fit79 = slope * 1979 + intercept
        fit20 = slope * 2020 + intercept
        rel = (fit20 - fit79) / fit79 if fit79 > 0 else np.nan
        p_rel = PROSSER_DJF[s][2]
        summary[s] = (fit79, fit20, rel, p_rel, r2)
        print(f"   {LABEL[s]:<6}{fit79:>9.3%}{fit20:>10.3%}{rel:>8.0%}"
              f"{p_rel:>10.0%}{rel / p_rel:>8.2f}{r2:>7.2f}")

    print(f"\n   {'':<6}{'ours (h)':>11}{'':>10}{'Prosser (h)':>13}")
    for s in ORDER:
        fit79, fit20, *_ = summary[s]
        print(f"   {LABEL[s]:<6}{fit79 * djf_days(1979) * 24:>11.1f}"
              f" ->{fit20 * djf_days(2020) * 24:>7.1f}"
              f"{PROSSER_DJF[s][0]:>10.1f} ->{PROSSER_DJF[s][1]:>6.1f}")

    # ------------------------------------------------------------------ verdict
    mog_rel, mog_p = summary["moderate"][2], PROSSER_DJF["moderate"][2]
    signs_ok = all(summary[s][2] > 0 for s in ORDER)
    monotone = all(summary[ORDER[i]][2] <= summary[ORDER[i + 1]][2] + 0.05
                   for i in range(len(ORDER) - 1))
    level_ratio = summary["moderate"][0] / (PROSSER_DJF["moderate"][0]
                                            / (djf_days(1979) * 24))

    print("\n" + "-" * 78)
    print(f"   1. ALL TRENDS POSITIVE       {'PASS' if signs_ok else 'FAIL'}")
    print(f"   2. MOG trend {mog_rel:.0%} vs Prosser {mog_p:.0%}   "
          f"ratio {mog_rel / mog_p:.2f}   "
          f"{'PASS' if 0.5 <= mog_rel / mog_p <= 2.0 else 'CHECK'}")
    print(f"   3. STRONGER AT HIGHER SEVERITY   {'PASS' if monotone else 'CHECK'}"
          "   (Prosser: 21/30/37/43/49 %)")
    print(f"\n   Level ratio at fitted 1979: {level_ratio:.2f}"
          f"   (single raw season gave 0.76)")
    print(f"   THE PREDICTION: the trend ratio should be CLOSER to 1.00 than")
    print(f"   the level ratio, because a constant offset cancels in a change.")
    print(f"   trend {abs(mog_rel / mog_p - 1):.0%} off vs level "
          f"{abs(level_ratio - 1):.0%} off -- "
          f"{'CONFIRMED' if abs(mog_rel / mog_p - 1) < abs(level_ratio - 1) else 'NOT confirmed'}")
    print("-" * 78)

    # ---------------------------------------------------------- per diagnostic
    if per_diag:
        print("\n" + "=" * 78)
        print("PER-DIAGNOSTIC MOG EXCEEDANCE — what the ensemble mean hides")
        print("=" * 78)
        print("   The 21-member mean can look right while one member is far off;")
        print("   a badly wrong diagnostic moves it by only 1/21. Prosser's")
        print("   Figure S5 is the published version of this breakdown.\n")
        a, b = YEARS[0], YEARS[-1]
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
