#!/usr/bin/env python
"""
ada/per_diagnostic_all_severities.py
====================================
(FINAL-CHECKS-2026-09-16)
The fitted per-diagnostic trend (annual, n = 42) at ALL FIVE severities in ONE
pass over the North Atlantic series -- plus a finite-value census of every
diagnostic in every year, which is the silent-failure sweep the final
selection asked for.

    pixi run python ada/per_diagnostic_all_severities.py \
        --derived-subdir derived/north_atlantic_audit \
        --thresholds $BASE/calibration/thresholds_2026-09-11.json

WHY A NEW SCRIPT RATHER THAN THREE RUNS OF per_diagnostic_trend.py
-----------------------------------------------------------------
`ada/per_diagnostic_trend.py` does one severity per run and reads every store
lazily twice per variable per year (once per diagnostic rate, once for the
ensemble): 72 min per severity on the annual season. Three severities would be
~3.6 h, over the 3 h budget for the final checks. This script loads each
year's box ONCE into memory and derives all five severities from that copy.

WHAT IS REUSED AND WHAT IS NOT
------------------------------
Reused by import from ada/full_trend_check.py, exactly as per_diagnostic_trend
does: the box, the cos(phi) weights, the year loader, the expected timestep
count, the OLS and the critical t. Reused from per_diagnostic_trend.py: the
STATUS 12.6 reference columns, so the CSVs have the identical format and
`ada/prosser_scorecard.py` / `ada/prosser_figures.py` read them unchanged.

NOT reused: `ftc.weighted_rate` and `aggregate.exceedance_*`. They are
re-expressed in numpy because the weights depend on latitude only, so

    sum(w * exceed) / sum(w over finite cells)
      = sum_lat w_lat * count_exceed_lat / sum_lat w_lat * count_finite_lat

which is the same arithmetic with the per-cell weight broadcast removed. The
ensemble keeps aggregate.py's populated-only semantics (Q-INTEG-3): per cell,
the mean over the diagnostics that are finite there, and the cell counts in
the denominator only if at least one is.

THE BUILT-IN PROOF THAT THIS IS THE SAME COMPUTATION
----------------------------------------------------
The MOG table this writes must equal the existing
cat_outputs/per_diagnostic_annual_moderate_audit.csv (job 1108148, written by
per_diagnostic_trend.py) to the printed precision. `ada/final_scorecard.py`
checks exactly that before it prints anything else. It was also checked
locally against per_diagnostic_trend.py on synthetic stores with NaN holes
(max difference 0 at printed precision) before this file was committed.

CHECKPOINTING
-------------
A TIMEOUT is this project's most expensive failure (STATUS 15.7 item 2): full
cost, zero output. Every finished year is appended to
<outdir>/_checkpoint_annual<tag>.csv and flushed. Resubmitting the same command
resumes after the last finished year, and refuses to resume if the series,
thresholds file or severity list differ from the checkpoint's.

OUTPUT (default outdir cat_outputs/final/, so nothing existing is overwritten)
    per_diagnostic_annual_<severity><tag>.csv          one per severity
    per_diagnostic_annual_<severity><tag>_series.csv   the 42 points
    finite_fraction_annual<tag>.csv                    year x diagnostic
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import resource
import socket
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
sys.path.insert(0, str(REPO))

import calibration                                   # noqa: E402

SEVERITIES = ["light", "light_to_moderate", "moderate",
              "moderate_to_severe", "severe"]


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def peak_rss_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def year_rates(ftc, base: Path, subdir: str, year: int, names: list[str],
               thresholds: dict, signs: dict, sevs: list[str]):
    """One year: per-diagnostic and ensemble weighted exceedance rates at every
    severity, plus the finite fraction of every diagnostic in the box."""
    ds = ftc.subset_box(ftc.load_year(base, year, subdir), **ftc.PROSSER_BOX)
    ds = ds[names].load()
    w_lat = np.asarray(ftc.lat_weights_for(ds).values, dtype=np.float64)

    n_time = ds.sizes["time"]
    expected = ftc.season_days(year, "annual") * 8
    if n_time != expected:
        print(f"   !! {year}: {n_time} timesteps, expected {expected}")

    rates = {s: {} for s in sevs}
    finite = {}
    k_sum = {s: None for s in sevs}     # per-cell count of exceeding diagnostics
    m_sum = None                        # per-cell count of finite diagnostics
    for n in names:
        da = ds[n]
        ilat = da.get_axis_num("latitude")
        arr = np.moveaxis(da.values, ilat, 0)          # (lat, ...)
        ok = np.isfinite(arr)
        other = tuple(range(1, arr.ndim))
        n_ok_lat = ok.sum(axis=other, dtype=np.int64)
        denom = float((w_lat * n_ok_lat).sum())
        finite[n] = float(ok.mean())
        m_sum = ok.astype(np.uint8) if m_sum is None else m_sum + ok
        for s in sevs:
            thr = thresholds[n][s]
            if signs[n] == "+":
                ex = arr >= thr                     # NaN compares False
            elif signs[n] == "-":
                ex = arr <= thr
            else:
                raise ValueError(f"{n}: sign {signs[n]!r} not supported here "
                                 f"(all 21 are '+', FORMULAS_AND_DECISIONS 3)")
            num = float((w_lat * ex.sum(axis=other, dtype=np.int64)).sum())
            rates[s][n] = num / denom if denom > 0 else float("nan")
            k_sum[s] = ex.astype(np.uint8) if k_sum[s] is None else k_sum[s] + ex
        del arr, ok

    ens = {}
    populated = m_sum > 0
    w_cell = np.broadcast_to(
        w_lat.reshape((-1,) + (1,) * (m_sum.ndim - 1)), m_sum.shape)
    denom = float(w_cell[populated].sum())
    for s in sevs:
        mean_cell = np.zeros(m_sum.shape, dtype=np.float64)
        np.divide(k_sum[s], m_sum, out=mean_cell, where=populated)
        ens[s] = float((mean_cell * w_cell).sum() / denom)
        del mean_cell
    del ds, k_sum, m_sum
    return rates, ens, finite, n_time


def main() -> int:
    ap = argparse.ArgumentParser(description="Per-diagnostic trend, all severities, one pass.")
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--derived-subdir", required=True)
    ap.add_argument("--thresholds", required=True)
    ap.add_argument("--severities", nargs="+", default=SEVERITIES, choices=SEVERITIES)
    ap.add_argument("--start-year", type=int, default=1979)
    ap.add_argument("--end-year", type=int, default=2020)
    ap.add_argument("--outdir", default=str(REPO / "cat_outputs" / "final"))
    ap.add_argument("--tag", default=None,
                    help="filename tag; default '_audit' for an audit series, '' otherwise")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    base = Path(args.base)
    ftc = _load("full_trend_check", REPO / "ada" / "full_trend_check.py")
    pdt = _load("per_diagnostic_trend", REPO / "ada" / "per_diagnostic_trend.py")
    dmod = _load("diagnostics", REPO / "2_diagnostics.py")

    tpath = Path(args.thresholds)
    if not tpath.is_absolute():
        tpath = base / tpath
    thresholds, signs, prov = calibration.load_thresholds(tpath)
    names = [k for k in dmod.REFERENCE_TABLE if k in thresholds]
    sevs = [s for s in SEVERITIES if s in args.severities]
    tag = args.tag if args.tag is not None else (
        "_audit" if "audit" in args.derived_subdir else "")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    years = list(range(args.start_year, args.end_year + 1))

    # ------------------------------------------------ configuration, first lines
    print(f">>> host        : {socket.gethostname()}  pid {os.getpid()}  "
          f"{time.strftime('%Y-%m-%dT%H:%M:%S')}")
    print(f">>> series      : {base / args.derived_subdir}")
    print(f">>> thresholds  : {tpath}")
    cal_dir = str(prov.get("derived_subdir", prov.get("period", "")))
    print(f"    calibrated on {cal_dir!r}, domain {prov.get('calibration_domain')}")
    print(f">>> severities  : {', '.join(sevs)}")
    print(f">>> years       : {years[0]}-{years[-1]} (n={len(years)}), annual, "
          f"box {ftc.PROSSER_BOX}")
    print(f">>> diagnostics : {len(names)}")
    print(f">>> output      : {outdir}  tag {tag!r}")

    # ------------------------------------------------ pairing guard, two ways
    if "derived/" not in cal_dir:
        print("!! thresholds provenance names no derived directory -- cannot verify "
              "the convention pairing. Refusing: the final checks must be paired.")
        return 2
    if ("audit" in args.derived_subdir) != ("audit" in cal_dir):
        print(f"!! CONVENTION MISMATCH -- series {args.derived_subdir!r} vs "
              f"thresholds calibrated on {cal_dir!r}. Refusing.")
        return 2
    first = base / args.derived_subdir / f"diagnostics_na_{years[0]}-01.zarr" / ".zattrs"
    try:
        fixes = json.loads(first.read_text()).get("audit_fixes", "baseline")
    except OSError as exc:
        print(f"!! cannot read {first}: {exc}")
        return 2
    print(f">>> store attrs : audit_fixes={fixes!r} ({first.parent.name})")
    if ("audit" in args.derived_subdir) != (fixes != "baseline"):
        print("!! the directory name and the stores' own audit_fixes disagree. Refusing.")
        return 2
    print("    pairing OK\n")

    # ------------------------------------------------ checkpoint
    ckpt = outdir / f"_checkpoint_annual{tag}.csv"
    meta_path = ckpt.with_suffix(".json")
    meta = dict(series=str(base / args.derived_subdir), thresholds=str(tpath),
                thresholds_mtime=tpath.stat().st_mtime, severities=sevs,
                names=names, box=str(ftc.PROSSER_BOX), years=years)
    done: dict[int, dict] = {}
    if ckpt.exists() and not args.no_resume:
        old = json.loads(meta_path.read_text()) if meta_path.exists() else None
        if old != meta:
            print(f"!! checkpoint {ckpt.name} was written for a different "
                  f"configuration. Move it aside or pass --no-resume.")
            return 2
        with ckpt.open() as fh:
            for r in csv.DictReader(fh):
                y = int(r["year"])
                d = done.setdefault(y, {"rates": {s: {} for s in sevs},
                                        "ens": {}, "finite": {}})
                if r["kind"] == "finite":
                    d["finite"] = {n: float(r[n]) for n in names}
                    d["n_time"] = int(float(r["ensemble"]))
                else:
                    d["rates"][r["severity"]] = {n: float(r[n]) for n in names}
                    d["ens"][r["severity"]] = float(r["ensemble"])
        print(f">>> resuming: {len(done)} year(s) already in {ckpt.name}\n")
    else:
        meta_path.write_text(json.dumps(meta, indent=1))
        with ckpt.open("w", newline="") as fh:
            csv.writer(fh).writerow(["year", "kind", "severity", *names, "ensemble"])

    t_start = time.time()
    todo = [y for y in years if y not in done]
    for i, year in enumerate(todo, 1):
        t0 = time.time()
        rates, ens, finite, n_time = year_rates(
            ftc, base, args.derived_subdir, year, names, thresholds, signs, sevs)
        with ckpt.open("a", newline="") as fh:
            wr = csv.writer(fh)
            for s in sevs:
                wr.writerow([year, "rate", s, *[f"{rates[s][n]:.12g}" for n in names],
                             f"{ens[s]:.12g}"])
            wr.writerow([year, "finite", "", *[f"{finite[n]:.12g}" for n in names],
                         n_time])
            fh.flush()
            os.fsync(fh.fileno())
        done[year] = dict(rates=rates, ens=ens, finite=finite, n_time=n_time)
        el = time.time() - t_start
        worst = min(finite, key=finite.get)
        print(f"   {year} done in {time.time() - t0:5.1f} s   "
              f"min finite {finite[worst]:.6f} ({worst})   rss {peak_rss_gb():.1f} G   "
              f"~{el / i * (len(todo) - i) / 60:.0f} min to go")

    print(f"\n   loaded in {(time.time() - t_start) / 60:.1f} min\n")

    # ------------------------------------------------ fits, identical to per_diagnostic_trend
    x = np.asarray(years, float)
    span = x[-1] - x[0]
    tc = ftc.tcrit(len(years) - 2)
    for s in sevs:
        out = {}
        for n in names:
            y = np.asarray([done[yr]["rates"][s][n] for yr in years], float)
            f = ftc.ols(x, y)
            y0 = f["slope"] * x[0] + f["intercept"]
            y1 = f["slope"] * x[-1] + f["intercept"]
            out[n] = dict(fit=f, y0=y0, y1=y1,
                          change=(y1 - y0) / y0 if y0 > 0 else np.nan,
                          lo=(f["slope"] - f["half"]) * span / y0 if y0 > 0 else np.nan,
                          hi=(f["slope"] + f["half"]) * span / y0 if y0 > 0 else np.nan)
        p = outdir / f"per_diagnostic_annual_{s}{tag}.csv"
        with p.open("w") as fh:
            fh.write("diagnostic,level_1979,level_2020,change,ci_lo,ci_hi,t,r2,"
                     "status_12_6_level,status_12_6_change\n")
            for n in out:
                o, f = out[n], out[n]["fit"]
                r = pdt.STATUS_12_6.get(n, ("", ""))
                fh.write(f"{n},{o['y0']:.8f},{o['y1']:.8f},{o['change']:.6f},"
                         f"{o['lo']:.6f},{o['hi']:.6f},{f['t']:.4f},{f['r2']:.4f},"
                         f"{r[0]},{r[1]}\n")
        ps = p.with_name(p.stem + "_series" + p.suffix)
        with ps.open("w") as fh:
            fh.write("year," + ",".join(names) + ",ensemble\n")
            for yr in years:
                fh.write(f"{yr}," + ",".join(f"{done[yr]['rates'][s][n]:.8f}" for n in names)
                         + f",{done[yr]['ens'][s]:.8f}\n")

        n_sig = sum(1 for n in out if np.isfinite(out[n]["fit"]["t"])
                    and abs(out[n]["fit"]["t"]) > tc)
        dgn = [n for n in out if not np.isfinite(out[n]["fit"]["t"])
               or not np.isfinite(out[n]["change"])]
        ens_fit = ftc.ols(x, np.asarray([done[yr]["ens"][s] for yr in years]))
        ens_y0 = ens_fit["slope"] * x[0] + ens_fit["intercept"]
        mean21 = float(np.mean([out[n]["y0"] for n in out]))
        print(f"   {s:<20} significant {n_sig}/21   degenerate {dgn or 'none'}   "
              f"mean-of-21 1979 {mean21:.4%} vs ensemble {ens_y0:.4%}   -> {p.name}")

    # ------------------------------------------------ finite census
    fp = outdir / f"finite_fraction_annual{tag}.csv"
    with fp.open("w") as fh:
        fh.write("year,n_time," + ",".join(names) + "\n")
        for yr in years:
            fh.write(f"{yr},{done[yr]['n_time']}," +
                     ",".join(f"{done[yr]['finite'][n]:.10f}" for n in names) + "\n")
    below = [(yr, n, done[yr]["finite"][n]) for yr in years for n in names
             if done[yr]["finite"][n] < 1.0]
    print(f"\nFINITE CENSUS — {len(years)} years x {len(names)} diagnostics in the box -> {fp.name}")
    if not below:
        print("   every diagnostic is finite in every cell of every year")
    else:
        by_diag: dict[str, list] = {}
        for yr, n, v in below:
            by_diag.setdefault(n, []).append((yr, v))
        for n, lst in by_diag.items():
            vals = [v for _, v in lst]
            print(f"   {n:<24} {len(lst):2d} year(s) below 1.0, min {min(vals):.6f}, "
                  f"max {max(vals):.6f}")
        print("   (a whole missing month shows as ~0.917; 2 NaN steps per month as ~0.992)")
    print(f"\n   wall {(time.time() - t_start) / 60:.1f} min, peak RSS {peak_rss_gb():.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
