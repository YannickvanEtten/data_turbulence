"""
ada/per_diagnostic_trend.py
===========================
The per-diagnostic trend, FITTED over all 42 years, with an interval.

WHY, WHEN full_trend_check.py ALREADY HAS --per-diagnostic
-----------------------------------------------------------
That flag reports the first and last year's exceedance and the percentage
between them. Two points. STATUS.md §12.8 retired exactly that construction
for the ensemble and wrote down why:

    "a comparison of point estimates from a small sample is not a test. Every
     future check in this project that compares a fitted quantity against a
     published one should carry an interval, or it is decoration."

A two-point change is the most extreme version of the thing §12.8 warns about,
and §12.4 measured the damage: across nine winters the MOG frequency spans a
factor of 1.65 while the fitted 41-year change is 43%. Interannual variability
is larger than the signal, so 1979 and 2020 are two draws from a noisy
distribution and their ratio is not a trend. This script fits all 42.

WHAT IT IS FOR
--------------
STATUS.md §12.6 is the project's per-diagnostic table and it is a NINE-SEASON,
DJF-ONLY, endpoint number. It is the basis for the top open item in §7:

    ubf is 15th of 21 on level but LAST on trend, at +11% where every sibling
    is +38% to +226%. It breaks the pattern that explains all the others.

That claim has never been checked against the full record — and the full
record has already overturned §12.5, which was built on the same nine seasons
and which predicted a "trend excess" that does not exist at n=42. So §12.6 is
under exactly the same suspicion, and `ubf` may not be anomalous at all.

This needs NO NEW DATA. The 504 North Atlantic stores already exist.

REUSES full_trend_check.py RATHER THAN REIMPLEMENTING
------------------------------------------------------
The box, the cos(phi) weighting, the season slicing, the year loader and the
OLS all come from `ada/full_trend_check.py` by import. Not for brevity: two
implementations of "Prosser's box" that disagree by one gridpoint would
produce two different per-diagnostic tables and no way to tell which was
right. STATUS §4g's lesson runs the other way too -- agreement between two
implementations is only as good as their shared assumptions, so where one
implementation will do, use one.

CROSS-CHECK BUILT IN
--------------------
The mean of the 21 per-diagnostic rates is printed next to the ensemble figure
that `aggregate.exceedance_mean_all_severities` produces. They should agree
closely but not exactly: the ensemble averages the 21 binary FIELDS per cell
before weighting, and f2d carries a different NaN mask from the other twenty
(STATUS §4d), so the two orderings of the average differ slightly. A large
gap would mean something is wrong with this script, not with the data.

USAGE
-----
    pixi run python ada/per_diagnostic_trend.py                    # annual, MOG
    pixi run python ada/per_diagnostic_trend.py --season djf       # compare to §12.6
    pixi run python ada/per_diagnostic_trend.py --severity severe
    pixi run python ada/per_diagnostic_trend.py --csv out.csv
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
sys.path.insert(0, str(REPO))

import aggregate                                     # noqa: E402
import calibration                                   # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# STATUS.md §12.6 — the nine-season, DJF-only, 1979-vs-2020 endpoint table this
# script exists to re-examine. Kept verbatim so the comparison is in the output
# rather than in someone's memory of a document.
STATUS_12_6 = {
    "vertical_wind_shear": (2.321, 49), "brown2": (1.993, 52), "ti1": (1.639, 51),
    "ti2": (1.566, 46), "temperature_gradient": (1.349, 68), "ngm1": (1.234, 82),
    "rva_magnitude": (0.857, 91), "nva": (0.827, 77), "endlich": (0.785, 95),
    "wind_speed": (0.713, 105), "deformation": (0.642, 58), "ncsu1": (0.610, 54),
    "brown1": (0.592, 67), "horizontal_divergence": (0.572, 38),
    "ubf": (0.249, 11), "vorticity_squared": (0.236, 99), "ngm2": (0.115, 56),
    "magnitude_pv": (0.075, 132), "colson_panofsky": (0.056, 148),
    "f2d": (0.026, 226), "negative_richardson": (0.022, 176),
}


def main() -> int:
    ap = argparse.ArgumentParser(description="Fitted per-diagnostic trend, n=42.")
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--thresholds", default=None,
                    help="default: newest calibration/thresholds_*.json")
    ap.add_argument("--season", default="annual",
                    choices=["djf", "mam", "jja", "son", "annual"])
    ap.add_argument("--severity", default="moderate",
                    choices=["light", "light_to_moderate", "moderate",
                             "moderate_to_severe", "severe"])
    ap.add_argument("--start-year", type=int, default=1979)
    ap.add_argument("--end-year", type=int, default=2020)
    ap.add_argument("--csv", default=None)
    args = ap.parse_args()

    base = Path(args.base)
    ftc = _load("full_trend_check", REPO / "ada" / "full_trend_check.py")

    tpath = (Path(args.thresholds) if args.thresholds
             else sorted((base / "calibration").glob("thresholds_*.json"))[-1])
    if not tpath.is_absolute():
        tpath = base / tpath
    thresholds, signs, prov = calibration.load_thresholds(tpath)

    dmod = _load("diagnostics", REPO / "2_diagnostics.py")
    names = [k for k in dmod.REFERENCE_TABLE if k in thresholds]
    years = list(range(args.start_year, args.end_year + 1))
    sev = args.severity

    print(f">>> thresholds: {tpath.name}")
    print(f"    domain {prov['calibration_domain']}, period {prov['period']}")
    print(f">>> {len(names)} diagnostics, {len(years)} years "
          f"({years[0]}-{years[-1]}), season {ftc.SEASON_TITLE[args.season]}, "
          f"severity {sev}")
    print(f">>> box {ftc.PROSSER_BOX}\n")

    rates: dict[str, list[float]] = {n: [] for n in names}
    ens: list[float] = []
    t0 = time.time()

    for year in years:
        ds_box = ftc.subset_box(ftc.load_year(base, year), **ftc.PROSSER_BOX)
        w = ftc.lat_weights_for(ds_box)
        mask = ds_box["time"].dt.month.isin(ftc.SEASON_MONTHS[args.season])
        ds_s = ds_box.isel(time=mask.values)

        expected = ftc.season_days(year, args.season) * 8
        if ds_s.sizes["time"] != expected:
            print(f"   !! {year}: {ds_s.sizes['time']} timesteps, expected "
                  f"{expected} — using actual, treat this year with caution")

        fields = {n: ds_s[n] for n in names}
        for n in names:
            rates[n].append(ftc.weighted_rate(
                aggregate.exceedance_field(fields[n], thresholds[n][sev],
                                           signs[n]), w))
        per_sev = aggregate.exceedance_mean_all_severities(
            fields, thresholds, [sev], signs)
        ens.append(ftc.weighted_rate(per_sev[sev]["exceedance_mean"], w))
        print(f"   {year}  done")
        del ds_box, ds_s, fields

    print(f"\n   loaded in {(time.time() - t0) / 60:.1f} min\n")

    x = np.asarray(years, float)
    out = {}
    for n in names:
        y = np.asarray(rates[n], float)
        f = ftc.ols(x, y)
        y0 = f["slope"] * x[0] + f["intercept"]
        y1 = f["slope"] * x[-1] + f["intercept"]
        span = x[-1] - x[0]
        out[n] = dict(fit=f, y0=y0, y1=y1,
                      change=(y1 - y0) / y0 if y0 > 0 else np.nan,
                      lo=(f["slope"] - f["half"]) * span / y0 if y0 > 0 else np.nan,
                      hi=(f["slope"] + f["half"]) * span / y0 if y0 > 0 else np.nan)

    title = f"PER-DIAGNOSTIC FITTED TREND — {ftc.SEASON_TITLE[args.season]}, {sev}, n={len(years)}"
    print("=" * 96)
    print(title)
    print("=" * 96)
    hdr = (f"   {'diagnostic':<23}{'1979':>9}{'2020':>9}{'change':>9}"
           f"{'95% CI':>18}{'t':>7}{'R2':>7}{'sig':>5}   §12.6")
    print(hdr)
    print("-" * 96)
    tc = ftc.tcrit(len(years) - 2)
    for n in sorted(out, key=lambda k: -out[k]["y0"]):
        o, f = out[n], out[n]["fit"]
        sig = "yes" if abs(f["t"]) > tc else "NO"
        ref = STATUS_12_6.get(n)
        refs = "{:+d}%".format(ref[1]) if ref else "--"
        ci = "[{:+.0%}, {:+.0%}]".format(o["lo"], o["hi"])
        print("   {:<23}{:>8.3%}{:>9.3%}{:>+9.0%}{:>18}{:>7.2f}{:>7.2f}{:>5}   {:>6}"
              .format(n, o["y0"], o["y1"], o["change"], ci,
                      f["t"], f["r2"], sig, refs))

    n_sig = sum(1 for n in out if abs(out[n]["fit"]["t"]) > tc)
    print("-" * 96)
    print(f"   significant at 5% (|t| > {tc:.2f}, df={len(years)-2}): "
          f"{n_sig}/{len(out)}")

    mean_of_21 = float(np.mean([out[n]["y0"] for n in out]))
    ens_fit = ftc.ols(x, np.asarray(ens, float))
    ens_y0 = ens_fit["slope"] * x[0] + ens_fit["intercept"]
    print(f"\n   CROSS-CHECK  mean of the 21 fitted 1979 levels : {mean_of_21:.4%}")
    print(f"                ensemble fitted 1979 level        : {ens_y0:.4%}")
    print(f"                difference                        : "
          f"{abs(mean_of_21 - ens_y0) / ens_y0:.2%}  "
          f"(small is expected — f2d's NaN mask differs, STATUS §4d)")

    # -------------------------------------------------------------- the point
    print("\n" + "=" * 96)
    print("ubf AGAINST STATUS.md §12.6")
    print("=" * 96)
    u = out.get("ubf")
    if u is None:
        print("   ubf absent from the threshold set — cannot judge.")
    else:
        changes = {n: out[n]["change"] for n in out if np.isfinite(out[n]["change"])}
        rank = sorted(changes, key=lambda k: changes[k]).index("ubf") + 1
        lvl_rank = sorted(out, key=lambda k: -out[k]["y0"]).index("ubf") + 1
        print(f"   §12.6 (n=9, DJF, endpoints): level 15/21, trend LAST at +11%,")
        print(f"                                 siblings +38% to +226%")
        print(f"   here  (n={len(years)}, {ftc.SEASON_TITLE[args.season]}, fitted): "
              f"level {lvl_rank}/21, trend rank {rank}/21 at "
              f"{u['change']:+.0%}")
        print(f"   95% CI [{u['lo']:+.0%}, {u['hi']:+.0%}], t={u['fit']['t']:.2f}, "
              f"{'significant' if abs(u['fit']['t']) > tc else 'NOT significant'}")
        vals = np.array(sorted(changes.values()), float)
        med = float(np.median(vals))
        mad = float(np.median(np.abs(vals - med)))
        print(f"   spread of the 21 fitted changes: {vals[0]:+.0%} to {vals[-1]:+.0%}, "
              f"median {med:+.0%}, MAD {mad:+.0%}")
        print()

        # NO AUTOMATED VERDICT HERE, DELIBERATELY.
        #
        # Two were tried and both failed on real data:
        #   1. "bottom-three rank AND change < median/2" -- decided by a
        #      knife-edge 0.18 < 0.18 on the DJF run.
        #   2. "does the 95% CI exclude the ensemble median" -- CONTAINS it on
        #      DJF (+36% median, CI [-5%, +41%]) and EXCLUDES it by ONE
        #      percentage point on Annual (+38% median, CI [+3%, +37%]).
        #      The same diagnostic on the same data flips verdict between two
        #      seasons. A test that does that is measuring rounding.
        #
        # STATUS §12.8's lesson generalises: report the statistic and let a
        # human read it. A robust z, (change - median)/MAD, is printed for
        # every diagnostic so `ubf` is judged against its peers rather than
        # against a threshold someone invented.
        print("   ROBUST OUTLIER SCORE  z = (change - median) / MAD")
        print("   (a diagnostic is worth investigating when |z| is large AND")
        print("    the same sign appears in more than one season)")
        ranked = sorted(changes, key=lambda k: (changes[k] - med) / mad if mad else 0)
        for n in ranked[:4]:
            z = (changes[n] - med) / mad if mad else float("nan")
            star = "  <-- " + n if n == "ubf" else ""
            print(f"     {n:<24}{changes[n]:+7.0%}   z={z:+5.2f}"
                  f"   t={out[n]['fit']['t']:5.2f}{star}")
        print(f"     ...")
        for n in ranked[-2:]:
            z = (changes[n] - med) / mad if mad else float("nan")
            print(f"     {n:<24}{changes[n]:+7.0%}   z={z:+5.2f}"
                  f"   t={out[n]['fit']['t']:5.2f}")
        print()
        print("   §12.6's TWO SPECIFIC CLAIMS, checked against this run:")
        print(f"     'ubf is LAST on trend'          -> it is rank {rank}/21 "
              f"({'TRUE' if rank == 1 else 'FALSE'})")
        print(f"     'siblings run +38% to +226%'    -> observed range "
              f"{vals[0]:+.0%} to {vals[-1]:+.0%} "
              f"({'TRUE' if vals[0] >= 0.38 else 'FALSE'})")

    if args.csv:
        p = Path(args.csv)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w") as fh:
            fh.write("diagnostic,level_1979,level_2020,change,ci_lo,ci_hi,t,r2,"
                     "status_12_6_level,status_12_6_change\n")
            for n in out:
                o, f = out[n], out[n]["fit"]
                r = STATUS_12_6.get(n, ("", ""))
                fh.write(f"{n},{o['y0']:.8f},{o['y1']:.8f},{o['change']:.6f},"
                         f"{o['lo']:.6f},{o['hi']:.6f},{f['t']:.4f},{f['r2']:.4f},"
                         f"{r[0]},{r[1]}\n")
        print(f"\nwrote {p}")
        # per-year series too: the table above is a fit, and anyone checking it
        # will want the points it was fitted to.
        ps = p.with_name(p.stem + "_series" + p.suffix)
        with ps.open("w") as fh:
            fh.write("year," + ",".join(names) + ",ensemble\n")
            for i, yr in enumerate(years):
                fh.write(f"{yr}," + ",".join(f"{rates[n][i]:.8f}" for n in names)
                         + f",{ens[i]:.8f}\n")
        print(f"wrote {ps}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
