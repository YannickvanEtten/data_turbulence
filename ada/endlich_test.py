#!/usr/bin/env python
"""
ada/endlich_test.py
===================
(ENDLICH-READING-2026-11-20)
Which reading of Sharman (2006) A25 should `endlich` use? Two readings of the
same formula, side by side, on data already on the share. Only endlich is
computed; nothing is written outside --out.

    pixi run python ada/endlich_test.py            # paths default to the share

THE QUESTION
    Endlich index, Sharman A25:  E = |V| * |d psi / dz|,  psi = wind direction.
    Sharman says only "centred second-order finite differences". On a 50 hPa
    stencil that admits two estimators of d psi/dz at 200 hPa:

    E0  "angle"       difference of the direction ANGLE between 175 and 225 hPa,
                      wrapped to the shorter arc, over dz (rojak; BOTH production
                      series use this -- DECISION_RECORD section 2).
    E1  "components"  the closed form d psi/dz = (u dv/dz - v du/dz) / |V|^2,
                      with du/dz, dv/dz centred and u, v at 200 hPa (audit F10).

    They agree as the stencil shrinks and disagree on 41-66 % of the exceedance
    set on ours (job 1102936). endlich is bit-identical across the baseline and
    audit conventions, so this choice is independent of everything else.

WHY IT MATTERS NOW
    endlich is class C on: level ratio vs Prosser 0.67 / 0.51 / 0.46
    (LOG / MOG / SOG, 1979), worse than its stencil-group median 0.93 / 0.85 /
    0.83, and at SOG his +36 % trend outside our interval (+58 % [+37, +78]).
    If it is to be used as a separate series, the reading has to be settled by
    measurement rather than inherited from rojak.

    Note: endlich has NO horizontal derivative, so the polar curvature defect
    (REFEREE_PASS F5) cannot reach it, whatever that pass's ranking suggests.

STAGES
    1  CONTROL: E0 on the first global file against
       derived/global_audit/diagnostics_glob_2000-01.zarr. Exit 2 on mismatch.
    2  global reference year, 48-day sub-sample: thresholds, zonal tail
       profile, NH/SH, Prosser-box density, flip rate E1 vs E0.
    3  North Atlantic 1979: annual level in Prosser's box vs his 1979 level.

PRE-REGISTERED, BEFORE THE NUMBERS
    Q1  E0 reproduces production: Spearman >= 0.99999.
    Q2  IF the angle reading causes endlich's excess deficit, then under E1:
          the NA 1979 level ratio moves toward its stencil-group median
          (0.93 / 0.85 / 0.83) at all three severities, and box density rises.
    Q3  FALSIFIED if E1/E0 box density is within 0.80-1.20 at BOTH MOG and SOG.
        Then the reading is not the cause; keep the angle reading (literal,
        and what both series use) and endlich stays class C.
    Q4  If E1 moves the level toward Prosser, the trend still has to be
        checked: ubf/endlich-only 42-year re-derive, scored against Prosser S4
        at three severities, then the class rule. Adopting E1 would make
        endlich a documented one-diagnostic departure from both series.

    Statistics, not verdicts.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "essentials"))
sys.path.insert(0, str(HERE))
import config as C                       # noqa: E402
import diagnostics as E                  # noqa: E402
from ubf_stencil_test import (BANDS, PROSSER_BOX, REPORT, SEVERITIES,  # noqa: E402
                              file_sha, lat_counts, peak_rss_gb, per_day,
                              spearman_sub, tail_thresholds)

VERSION = "ENDLICH-READING-2026-11-20"
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
# Prosser (2023) Fig. S4-LOG, Fig. 4, Fig. S4-SOG, endlich 1979 level (%); our
# production fitted level ratio (FINAL_DATASET.md section 5).
PROSSER_1979 = {"light": 5.25, "moderate": 0.92, "severe": 0.23}
PRODUCTION_RATIO = {"light": 0.67, "moderate": 0.51, "severe": 0.46}
STENCIL_GROUP_MEDIAN = {"light": 0.93, "moderate": 0.85, "severe": 0.83}
VERSIONS = ("E0", "E1")


def endlich_versions(ds: xr.Dataset) -> dict:
    """endlich at 200 hPa under both readings, via essentials.endlich itself."""
    u, v, z = ds["eastward_wind"], ds["northward_wind"], ds["geopotential"]
    s = dict(u=u, v=v, z=z, du_dz=E.d_dz(u, z), dv_dz=E.d_dz(v, z))
    out = {}
    for name, reading in (("E0", "angle"), ("E1", "components")):
        s["conv"] = {"endlich": reading}
        out[name] = E.endlich(s).sel(pressure_level=C.TARGET_LEVEL_HPA)
    return out


def load(path) -> xr.Dataset:
    ds = E.prepare(E.open_era5(path))
    return ds[["eastward_wind", "northward_wind", "geopotential"]]


def stage1(global_file: Path, control_zarr: Path) -> dict:
    print("=" * 78 + "\nSTAGE 1  CONTROL: E0 against the production store\n" + "=" * 78)
    ds = load(global_file)
    z = xr.open_zarr(control_zarr)["endlich"]
    tdim = next(d for d in ("time", "valid_time") if d in z.dims)
    if tdim != "time":
        z = z.rename({tdim: "time"})
    if float(z["longitude"].max()) > 180.0:
        z = z.assign_coords(longitude=(((z["longitude"] + 180) % 360) - 180)).sortby("longitude")
    times = np.intersect1d(ds["time"].values, z["time"].values)[:8]
    if times.size == 0:
        print("!! no common timesteps"); return {"ok": False}
    mine = endlich_versions(ds.sel(time=times).load())["E0"]
    ref = z.sel(time=times).load()
    mine, ref = xr.align(mine, ref.transpose(*mine.dims), join="inner")
    a = mine.values.astype(np.float64); b = ref.values.astype(np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    rel = np.abs(a[ok] - b[ok]) / np.maximum(np.abs(b[ok]), 1e-30)
    rho = spearman_sub(a, b)
    out = dict(ok=bool(rho >= 0.99999), n_steps=int(times.size), spearman=rho,
               rel_median=float(np.median(rel)), rel_p99=float(np.percentile(rel, 99)))
    print(f"  {times.size} timesteps, {ok.sum():,} cells   Spearman {rho:.7f}   "
          f"rel diff median {out['rel_median']:.2e}  p99 {out['rel_p99']:.2e}")
    print("  Q1 requires Spearman >= 0.99999" + ("" if out["ok"] else "   -- NOT MET, stopping"))
    return out


def stage2(files):
    print("\n" + "=" * 78 + "\nSTAGE 2  global reference year (48-day sub-sample)\n" + "=" * 78)
    chunks = {n: [] for n in VERSIONS}
    lat = lon = None
    t0 = time.time()
    for k, path in enumerate(files):
        for day in per_day(load(path)):
            if lat is None:
                lat = day["latitude"].values.astype(np.float64)
                lon = day["longitude"].values.astype(np.float64)
            out = endlich_versions(day)
            for n in VERSIONS:
                chunks[n].append(out[n].transpose("latitude", "longitude", "time").values.astype(np.float32))
        print(f"  [{k + 1:>2}/{len(files)}] {Path(path).name}  "
              f"{(time.time() - t0) / 60:5.1f} min  rss {peak_rss_gb():.1f} GB")
    data = {n: np.concatenate(chunks[n], axis=2) for n in VERSIONS}
    del chunks
    w_lat = np.cos(np.deg2rad(lat)); w_lat[np.abs(lat) >= 90.0] = 0.0
    thr = {n: dict(zip(SEVERITIES, map(float, tail_thresholds(data[n], w_lat, list(SEVERITIES.values())))))
           for n in VERSIONS}

    lon_in = (lon >= PROSSER_BOX["lon"][0]) & (lon <= PROSSER_BOX["lon"][1])
    lat_in = (lat >= PROSSER_BOX["lat"][0]) & (lat <= PROSSER_BOX["lat"][1])
    area_band = np.array([w_lat[(lat >= lo) & (lat < hi)].sum() for _, lo, hi in BANDS]) / w_lat.sum()
    nt = data["E0"].shape[2]
    w_all = w_lat.sum() * lon.size * nt
    w_box = w_lat[lat_in].sum() * lon_in.sum() * nt
    stats, exceed = {}, {}
    for n in VERSIONS:
        for sev in REPORT:
            m = data[n] >= thr[n][sev]
            exceed[(n, sev)] = m
            mass_lat = lat_counts(m) * w_lat
            gfrac = mass_lat.sum() / w_all
            band = np.array([mass_lat[(lat >= lo) & (lat < hi)].sum() for _, lo, hi in BANDS]) / mass_lat.sum()
            box = float((m[lat_in][:, lon_in].reshape(lat_in.sum(), -1).sum(axis=1) * w_lat[lat_in]).sum())
            stats[(n, sev)] = dict(
                version=n, severity=sev, threshold=thr[n][sev], global_pct=100 * gfrac,
                **{f"density_{b}": float(d) for (b, _, _), d in zip(BANDS, band / area_band)},
                nh_sh=float(mass_lat[lat > 0].sum() / max(mass_lat[lat < 0].sum(), 1e-30)),
                box_pct=100 * box / w_box, box_density=float(box / w_box / gfrac))
    rho = spearman_sub(data["E1"], data["E0"])
    for sev in REPORT:
        a, b = exceed[("E1", sev)], exceed[("E0", sev)]
        flip = (lat_counts(a ^ b) * w_lat).sum() / max((lat_counts(a | b) * w_lat).sum(), 1e-30)
        b0 = stats[("E0", sev)]["box_density"]
        stats[("E1", sev)].update(spearman_vs_E0=rho, flip_vs_E0=float(flip),
                                  box_density_vs_E0=stats[("E1", sev)]["box_density"] / b0 if b0 > 0 else float("nan"))

    print(f"\n  {nt} timesteps x {lat.size} x {lon.size};  rss {peak_rss_gb():.1f} GB")
    print("  tail DENSITY by band (1.00 = even), NH/SH mass ratio, Prosser-box density")
    print(f"  {'ver':<4}{'sev':<10}{'global%':>8}" + "".join(f"{b:>9}" for b, _, _ in BANDS)
          + f"{'NH/SH':>7}{'box':>7}{'rho':>9}{'flip':>7}")
    for sev in REPORT:
        for n in VERSIONS:
            s = stats[(n, sev)]
            print(f"  {n:<4}{sev:<10}{s['global_pct']:>8.3f}"
                  + "".join(f"{s['density_' + b]:>9.2f}" for b, _, _ in BANDS)
                  + f"{s['nh_sh']:>7.2f}{s['box_density']:>7.2f}"
                  + (f"{s['spearman_vs_E0']:>9.5f}{100 * s['flip_vs_E0']:>6.1f}%" if n == "E1" else ""))
        print()
    print("  box density, E1 / E0:  " + "   ".join(
        f"{sev} {stats[('E1', sev)]['box_density_vs_E0']:.2f}" for sev in REPORT))
    print("  Q3: falsified if E1/E0 box density within 0.80-1.20 at BOTH moderate and severe")
    return thr, stats


def stage3(files, thr) -> dict:
    print("\n" + "=" * 78 + "\nSTAGE 3  North Atlantic 1979, Prosser's box, each reading's own thresholds\n" + "=" * 78)
    num = {(n, s): 0.0 for n in VERSIONS for s in REPORT}
    den = 0.0
    t0 = time.time()
    lat_in = None
    for k, path in enumerate(files):
        for day in per_day(load(path)):
            if lat_in is None:
                lat = day["latitude"].values; lon = day["longitude"].values
                lat_in = (lat >= PROSSER_BOX["lat"][0]) & (lat <= PROSSER_BOX["lat"][1])
                lon_in = (lon >= PROSSER_BOX["lon"][0]) & (lon <= PROSSER_BOX["lon"][1])
                w = np.cos(np.deg2rad(lat[lat_in]))
            out = endlich_versions(day)
            den += w.sum() * lon_in.sum() * day.sizes["time"]
            for n in VERSIONS:
                x = out[n].transpose("latitude", "longitude", "time").values[lat_in][:, lon_in]
                for s in REPORT:
                    num[(n, s)] += float(((x >= thr[n][s]).reshape(x.shape[0], -1).sum(axis=1) * w).sum())
        print(f"  [{k + 1:>2}/{len(files)}] {Path(path).name}  {(time.time() - t0) / 60:5.1f} min")
    res = {}
    print(f"\n  {'ver':<4}" + "".join(f"{s + ' %':>16}" for s in REPORT)
          + "".join(f"{'ratio ' + s[:3]:>12}" for s in REPORT))
    for n in VERSIONS:
        lv = {s: 100 * num[(n, s)] / den for s in REPORT}
        rt = {s: lv[s] / PROSSER_1979[s] for s in REPORT}
        res[n] = dict(level=lv, ratio_vs_prosser=rt)
        print(f"  {n:<4}" + "".join(f"{lv[s]:>16.4f}" for s in REPORT) + "".join(f"{rt[s]:>12.3f}" for s in REPORT))
    print("  Prosser 1979  " + "  ".join(f"{s} {PROSSER_1979[s]} %" for s in REPORT))
    imp = {s: (PRODUCTION_RATIO[s] * res["E1"]["level"][s] / res["E0"]["level"][s]
               if res["E0"]["level"][s] > 0 else float("nan")) for s in REPORT}
    res["E1"]["implied_fitted_ratio"] = imp
    print("\n  implied FITTED level ratio = production fitted ratio x (E1 / E0, raw 1979):")
    print("  E1   " + "   ".join(f"{s} {PRODUCTION_RATIO[s]:.2f} -> {imp[s]:.2f}" for s in REPORT))
    print("  stencil-group median       " + "   ".join(f"{s} {STENCIL_GROUP_MEDIAN[s]:.2f}" for s in REPORT))
    print("  (class rule band for a level ratio: 0.75-1.33)")
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--global-glob", default="raw/global_sub48/*.grib")
    ap.add_argument("--na-glob", default="raw/north_atlantic/era5_na_1979-*.grib")
    ap.add_argument("--control-zarr", default="derived/global_audit/diagnostics_glob_2000-01.zarr")
    ap.add_argument("--out", default=str(REPO / "cat_outputs" / "endlich_reading"))
    ap.add_argument("--stages", default="123")
    args = ap.parse_args()

    base = Path(args.base)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    print(f">>> {VERSION}   script sha256 {file_sha(Path(__file__))}")
    gfiles = sorted(glob.glob(str(base / args.global_glob)))
    nfiles = sorted(glob.glob(str(base / args.na_glob)))
    print(f">>> global files {len(gfiles)}  NA files {len(nfiles)}\n")
    summary = {"version": VERSION}
    if "1" in args.stages:
        if not gfiles:
            print("!! no global files"); return 2
        summary["stage1"] = stage1(Path(gfiles[0]), base / args.control_zarr)
        if not summary["stage1"].get("ok"):
            (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str)); return 2
    if "2" in args.stages:
        if len(gfiles) != 12:
            print(f"!! expected 12 global files, found {len(gfiles)}"); return 2
        thr, stats = stage2(gfiles)
        summary["thresholds"] = thr
        rows = list(stats.values())
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(out / "global_tail.csv", "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader()
            for r in rows:
                wr.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
        if "3" in args.stages:
            if len(nfiles) != 12:
                print(f"!! expected 12 NA files, found {len(nfiles)}"); return 2
            summary["stage3"] = stage3(nfiles, thr)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\n>>> wrote {out}/summary.json and global_tail.csv   peak rss {peak_rss_gb():.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
