#!/usr/bin/env python
"""
ada/check_shape_ratios.py
=========================
The p97/median shape test across all 21. FORMULA_AUDIT.md §7.1.

    pixi run python ada/check_shape_ratios.py                     # defaults
    pixi run python ada/check_shape_ratios.py --box prosser --no-daily-mean

WHY THIS EXISTS
---------------
CALIBRATION_REFERENCE.md §4.1 concludes that the two Williams tables cannot
anchor MAGNITUDE, because a 2 deg GCM and 0.25 deg ERA5 differ in every
gradient-based quantity -- and then sets them aside. That gives up too much.

The two tables come from the SAME model, box, season and sampling:

    Williams (2017) Table 2   p97 onset thresholds, 21 diagnostics
    Williams & Joshi (2013)   medians,              21 diagnostics
    both: GFDL-CM2.1, DJF, 200 hPa, 50-75N/10-60W, DAILY-MEAN fields

so their RATIO is a pure distribution-shape statistic. It is dimensionless,
which is the point: it survives the unit conversions that make magnitude
comparisons fragile, and it is far more transferable across resolution than
either number alone. No scale factors appear anywhere in this script, because
a ratio of two quantiles of the same field does not need any.

It also targets the right thing. STATUS.md §1 says the econometric direction
is TAIL-BEHAVIOUR modelling; a test of distribution shape is closer to that
than any test of level.

WHAT THE PUBLISHED TABLE LOOKS LIKE, AND WHY IT IS INTERNALLY CONSISTENT
-----------------------------------------------------------------------
The published ratios are not scattered -- they track algebraic degree:

    degree 1  wind_speed, deformation, temperature_gradient, vws     2.5-2.8
    degree 2  ti1, ti2, ngm1, ngm2, nva, rva  (products of two)      6.2-6.6
              and 2.7^2 = 7.3
    square    vorticity_squared                                      11.1
    cubic     ncsu1 (also clipped at zero, which crushes the median) 108

That structure is what makes the test diagnostic rather than decorative: a
diagnostic whose ratio is off by 10x while its algebraic siblings are off by
1.5x is flagged, whatever the absolute values do.

PREDICTIONS, STATED BEFORE RUNNING SO THE TEST CAN FAIL
-------------------------------------------------------
  1. Most diagnostics land within a factor of ~2 of the published ratio.
  2. The SHEAR-BASED ones come in systematically LOW, because the 50 hPa
     stencil damps the tail (CALIBRATION_REFERENCE.md §10.3). That would be a
     second, independent measurement of the stencil effect -- from
     distribution shape rather than from exceedance frequency.
  3. `f2d` under variant A comes in FAR TOO HIGH -- hundreds, not 13.6 --
     because a signed material derivative is centred on zero and its median
     collapses. See FORMULA_AUDIT.md §4.3. If it does, variant A is refuted
     independently of anything in the paper.

CAVEATS THAT BELONG IN ANY WRITE-UP
-----------------------------------
  - Williams averages the WIND AND TEMPERATURE to daily means and then
    computes the diagnostic; --daily-mean (the default) averages the
    DIAGNOSTIC. For a non-linear diagnostic those differ. It is the closer of
    the two available approximations, not an equivalence. Run
    --no-daily-mean too: the gap between the two columns bounds the size of
    the approximation.
  - Resolution does affect distribution shape, so this is a RELATIVE screen
    across the 21, not an absolute test of any one of them.
  - Two diagnostics are negative-valued (colson_panofsky, negative_richardson)
    and a quantile ratio is not a shape statistic for them. They get a
    standardised spread instead, and no published comparison, because the
    published IQR was never tabulated.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _checkutil import (BASE, PROSSER_BOX, WILLIAMS_BOX, cos_phi_weights,  # noqa: E402
                        daily_mean, load_module, peak_rss_gb, subset_box)
from calib_weighted_percentile import weighted_percentile  # noqa: E402

# Williams (2017) Table 2 "Light" column = p97 onset threshold.
# Williams & Joshi (2013) Table 1 = median.
# Both in the SAME units per diagnostic, so only the ratio is kept here --
# see the module docstring on why no scale factors are needed.
#
# NVA and RVA: Williams prints 10^-9 s^-2 where W&J print 10^-10 s^-2, so the
# medians below are converted to 10^-9 (2.33 -> 0.233, 2.05 -> 0.205). Getting
# this backwards inflated those two by 100x once already -- STATUS.md §11.5 #3.
PUBLISHED = {
    #                        p97 (W2017)  median (W&J)
    "magnitude_pv":          (8.33,       6.84),
    "brown1":                (99.0,       77.1),
    "temperature_gradient":  (14.7,       5.75),
    "horizontal_divergence": (11.9,       2.82),
    "vertical_wind_shear":   (5.3,        1.88),
    "endlich":               (3.21,       0.952),
    "deformation":           (50.9,       18.6),
    "wind_speed":            (40.9,       14.9),
    "ngm2":                  (53.0,       8.17),
    "rva_magnitude":         (1.44,       0.233),
    "ubf":                   (1230.0,     161.0),
    "nva":                   (1.33,       0.205),
    "brown2":                (870.0,      116.0),
    "vorticity_squared":     (2.46,       0.221),
    "ti1":                   (195.0,      31.5),
    "ngm1":                  (1.65,       0.251),
    "ti2":                   (184.0,      28.8),
    "f2d":                   (770.0,      56.6),
    "ncsu1":                 (1200.0,     11.1),
}

# Negative-valued: a quantile ratio is meaningless, report spread instead.
SIGNED = ["negative_richardson", "colson_panofsky"]

# Algebraic degree, for the family-structure column. Purely descriptive.
DEGREE = {
    "wind_speed": 1, "deformation": 1, "temperature_gradient": 1,
    "vertical_wind_shear": 1, "horizontal_divergence": 1, "magnitude_pv": 1,
    "brown1": 1, "endlich": 1,
    "ti1": 2, "ti2": 2, "ngm1": 2, "ngm2": 2, "nva": 2, "rva_magnitude": 2,
    "vorticity_squared": 2, "ubf": 2, "f2d": 3, "brown2": 3, "ncsu1": 3,
}


def month_paths(base: Path, domain: str, year: int, months) -> list[Path]:
    d = base / "derived" / domain
    if domain == "global":
        # Calibration months carry no day suffix in the derived name.
        return [d / f"diagnostics_glob_{year}-{m}.zarr" for m in months]
    return [d / f"diagnostics_na_{year}-{m}.zarr" for m in months]


def stats_for(values: np.ndarray, weights: np.ndarray) -> dict:
    ok = np.isfinite(values)
    v, w = values[ok], weights[ok]
    if v.size < 100:
        return dict(n=int(v.size), p50=np.nan, p97=np.nan,
                    p25=np.nan, p75=np.nan)
    q = weighted_percentile(v, w, np.array([25.0, 50.0, 75.0, 97.0]))
    return dict(n=int(v.size), p25=float(q[0]), p50=float(q[1]),
                p75=float(q[2]), p97=float(q[3]))


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--domain", default="global", choices=["global", "north_atlantic"],
                    help="global (default) is the only one that covers "
                         "Williams' 50-75N box; the NA download stops at 60N.")
    ap.add_argument("--year", type=int, default=2000)
    ap.add_argument("--months", default="01,02,12",
                    help="DJF by default, matching Williams' season")
    ap.add_argument("--box", default="williams", choices=["williams", "prosser"])
    ap.add_argument("--no-daily-mean", action="store_true",
                    help="use the 3-hourly instantaneous fields instead of "
                         "daily means. Run BOTH: the difference between them "
                         "bounds the sampling approximation.")
    args = ap.parse_args()

    base = Path(args.base)
    months = [m.strip() for m in args.months.split(",")]
    box = WILLIAMS_BOX if args.box == "williams" else PROSSER_BOX

    paths = month_paths(base, args.domain, args.year, months)
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        print(f"!! missing derived months: {missing}")
        print(f"   looked in {base / 'derived' / args.domain}")
        return 1

    print(f">>> domain {args.domain}, {args.year} months {months}")
    print(f">>> box {args.box}: lat {box['lat']}, lon {box['lon']}")
    print(f">>> sampling: {'3-hourly instantaneous' if args.no_daily_mean else 'daily mean'}")

    parts = []
    for p in paths:
        ds = subset_box(xr.open_zarr(p), **box)
        if not args.no_daily_mean:
            ds = daily_mean(ds)
        parts.append(ds.load())
        print(f"    {p.name}: {dict(parts[-1].sizes)}   [rss {peak_rss_gb():.1f} GB]")
    ds = xr.concat(parts, dim="time")
    del parts

    variant = ds.attrs.get("f2d_variant", "unrecorded (pre-2026-08-29)")
    conv = ds.attrs.get("deformation_convention", "unrecorded — assume DEF^2")
    print(f"    f2d variant: {variant}")
    print(f"    deformation: {conv}\n")

    # Canonical W&J ordering when 2_diagnostics can be imported; otherwise
    # alphabetical. Only the ORDER depends on it, so a missing rojak (running
    # this on a laptop against a copied zarr, say) degrades to a cosmetic
    # difference rather than to a failure.
    try:
        dmod = load_module("diagnostics", "2_diagnostics.py")
        names = [k for k in dmod.REFERENCE_TABLE if k in ds.data_vars]
    except Exception as exc:                       # noqa: BLE001
        print(f"    (2_diagnostics unavailable -- {type(exc).__name__}; "
              f"falling back to alphabetical order)")
        names = sorted(ds.data_vars)

    results = {}
    for k in names:
        da = ds[k]
        w = cos_phi_weights(da)
        v = np.asarray(da.values, dtype=np.float64).ravel()
        results[k] = stats_for(v, w)

    # ------------------------------------------------------------------ report
    print("=" * 84)
    print("DISTRIBUTION SHAPE — p97/median, ours vs Williams (2017)/W&J (2013)")
    print("=" * 84)
    print(f"   {'diagnostic':<23}{'deg':>4}{'our p97/med':>13}"
          f"{'published':>11}{'ours/pub':>10}   note")

    flagged = []
    rows = sorted((k for k in names if k in PUBLISHED),
                  key=lambda k: -(PUBLISHED[k][0] / PUBLISHED[k][1]))
    for k in rows:
        s = results[k]
        pub = PUBLISHED[k][0] / PUBLISHED[k][1]
        if not np.isfinite(s["p50"]) or s["p50"] == 0:
            print(f"   {k:<23}{DEGREE.get(k, 0):>4}{'—':>13}{pub:>11.2f}{'—':>10}"
                  f"   median is zero: ratio undefined")
            flagged.append(k)
            continue
        ours = s["p97"] / s["p50"]
        rel = ours / pub
        note = ""
        if not np.isfinite(rel):
            note = "undefined"
        elif rel > 5 or rel < 0.2:
            note = "<<< OFF BY MORE THAN 5x"
            flagged.append(k)
        elif rel > 2 or rel < 0.5:
            note = "< outside 2x"
        print(f"   {k:<23}{DEGREE.get(k, 0):>4}{ours:>13.2f}{pub:>11.2f}"
              f"{rel:>10.2f}   {note}")

    print("\n" + "-" * 84)
    print("NEGATIVE-VALUED — no published shape statistic exists for these")
    print(f"   {'diagnostic':<23}{'p25':>13}{'median':>13}{'p75':>13}{'p97':>13}")
    for k in SIGNED:
        if k not in results:
            continue
        s = results[k]
        print(f"   {k:<23}{s['p25']:>13.4g}{s['p50']:>13.4g}"
              f"{s['p75']:>13.4g}{s['p97']:>13.4g}")
    print("   Reported for the record and for comparison across runs — a")
    print("   quantile RATIO across a sign change carries no information")
    print("   (CALIBRATION_REFERENCE.md §10.4 records that mistake being made).")

    # ------------------------------------------------------------------ reading
    shear_family = ["vertical_wind_shear", "ti1", "ti2", "ngm2", "brown2"]
    have = [k for k in shear_family if k in results and
            np.isfinite(results[k]["p50"]) and results[k]["p50"] != 0]
    if have:
        rels = [(results[k]["p97"] / results[k]["p50"]) /
                (PUBLISHED[k][0] / PUBLISHED[k][1]) for k in have]
        print("\n" + "-" * 84)
        print("PREDICTION 2 — shear-based diagnostics should be systematically LOW")
        print(f"   {', '.join(have)}")
        print(f"   median ours/published across that family: {np.median(rels):.2f}")
        print("   Below 1.0 supports the stencil-damping reading of")
        print("   CALIBRATION_REFERENCE.md §10.3, measured from shape rather")
        print("   than from exceedance frequency. Around 1.0 does not refute it")
        print("   — it says shape is less sensitive to the stencil than the")
        print("   exceedance tail is, which is itself worth knowing.")

    if "f2d" in results and np.isfinite(results["f2d"]["p50"]):
        s = results["f2d"]
        ratio = s["p97"] / s["p50"] if s["p50"] else float("inf")
        print("\n" + "-" * 84)
        print("PREDICTION 3 — f2d")
        print(f"   variant in this data: {variant}")
        print(f"   our median {s['p50']:.4g}, our p97 {s['p97']:.4g}, "
              f"ratio {ratio:.1f} vs published 13.6")
        print("   A ratio in the hundreds means the field is centred on zero,")
        print("   i.e. a SIGNED material derivative — which is not what the")
        print("   published pair describes. See FORMULA_AUDIT.md §4.3 and run")
        print("   ada/check_f2d_variants.py to compare all four readings.")

    print("\n" + "=" * 84)
    if flagged:
        print(f"   FLAGGED: {', '.join(flagged)}")
        print("   Chase these individually; the ensemble mean of 21 hides each")
        print("   of them at about 1/21 of its size.")
    else:
        print("   No diagnostic outside 5x of its published shape.")
    print(f"   PEAK RSS {peak_rss_gb():.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
