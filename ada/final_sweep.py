#!/usr/bin/env python
"""
ada/final_sweep.py
==================
(FINAL-CHECKS-2026-09-16)
Provenance, structure and silent-failure sweep of the FINAL dataset: every
store of `derived/north_atlantic_audit` (504) and `derived/global_audit` (12),
the paired thresholds file, and the job logs that wrote the stores.

    pixi run python ada/final_sweep.py              # defaults = the final selection

WHY THIS AND NOT jobs/15
    ada/check_production_complete.py hardcodes derived/north_atlantic and
    derived/global and checks only f2d_variant and deformation_convention. It
    has never been run on the audit directories, and it does not look at
    `audit_fixes` -- the attribute that distinguishes the two series. As of
    2026-09-16 the record verified the count, the size and ONE store's
    attributes of the audit series. This checks all of them.

WHAT IS CHECKED, PER STORE (metadata, plus one small read per variable)
    1. present, with consolidated metadata (.zmetadata) -- a store killed
       mid-write has none
    2. .zattrs: audit_fixes, f2d_variant, deformation_convention,
       target_level_hPa, and source_file naming the right month
    3. exactly the 21 expected variables, float32, on the expected grid
    4. time axis: 8 steps x days-in-month, starting at 00 UTC on the 1st,
       uniformly 3-hourly
    5. chunk census: chunk files present versus the number the array shape
       implies. zarr-python 3 does not write a chunk that is entirely fill
       value (NaN), so a shortfall means an all-NaN day-block or an unwritten
       one. If this zarr build does write empty chunks, the census is simply
       uninformative -- it cannot produce a false alarm.
    6. one timestep read per variable (the middle of the month), which must be
       finite everywhere; for f2d also the first and last step, to record how
       the month edges are handled (finite = one-sided time difference,
       NaN = dropped). The 42-year finite census of every value in the
       Prosser box is done by ada/per_diagnostic_all_severities.py.

THRESHOLDS FILE
    provenance names the paired reference directory; 21 diagnostics; five
    finite, strictly increasing values each (the one-tailed-upper property of
    FORMULAS_AND_DECISIONS 3, now checked on OUR ladders); sample sizes.

LOGS
    For every store, the most recent job log that WROTE it (a skip writes
    nothing and is ignored):
    exit status, '21 of 21', the fix label printed before computing, and any
    '!!' line. A task that failed would have exited 2 and been visible in
    sacct; this confirms it store by store rather than trusting the job-level
    summary.

Exit 0 when nothing is listed, 1 otherwise. Every listed item is printed with
its store name, so the output is the repair list.
"""
from __future__ import annotations

import argparse
import calendar
import csv
import importlib.util
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
EXPECT_FIXES = "meridional_metric+four_variable_provenance"
SEVERITIES = ["light", "light_to_moderate", "moderate",
              "moderate_to_severe", "severe"]


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def expected_stores(kind: str):
    if kind == "na":
        return [(f"diagnostics_na_{y}-{m:02d}.zarr", y, m, f"era5_na_{y}-{m:02d}.grib")
                for y in range(1979, 2021) for m in range(1, 13)]
    return [(f"diagnostics_glob_2000-{m:02d}.zarr", 2000, m, f"2000-{m:02d}")
            for m in range(1, 13)]


def chunk_census(store: Path, zmeta: dict, var: str) -> tuple[int, int]:
    za = zmeta["metadata"].get(f"{var}/.zarray")
    if za is None:
        return -1, -1
    expected = math.prod(math.ceil(s / c) for s, c in zip(za["shape"], za["chunks"]))
    vdir = store / var
    present = sum(1 for p in vdir.rglob("*") if p.is_file() and not p.name.startswith("."))
    return present, expected


def check_store(path: Path, year: int, month: int, src_expect: str, names: list[str],
                grid: tuple[int, int], fixes: str) -> tuple[list[str], dict]:
    problems: list[str] = []
    info: dict = {}
    if not path.exists():
        return ["MISSING"], info
    zm = path / ".zmetadata"
    if not zm.exists():
        return ["no .zmetadata (incomplete write)"], info
    zmeta = json.loads(zm.read_text())
    attrs = zmeta["metadata"].get(".zattrs", {})
    info["audit_fixes"] = attrs.get("audit_fixes", "<absent>")
    if attrs.get("audit_fixes", "baseline") != fixes:
        problems.append(f"audit_fixes={attrs.get('audit_fixes')!r}")
    if attrs.get("f2d_variant") != "A":
        problems.append(f"f2d_variant={attrs.get('f2d_variant')!r}")
    if not str(attrs.get("deformation_convention", "")).startswith("DEF (un-squared"):
        problems.append(f"deformation_convention={attrs.get('deformation_convention')!r}")
    if int(attrs.get("target_level_hPa", -1)) != 200:
        problems.append(f"target_level_hPa={attrs.get('target_level_hPa')!r}")
    src = str(attrs.get("source_file", ""))
    info["source_file"] = src
    if src_expect not in src:
        problems.append(f"source_file={src!r} (expected to contain {src_expect!r})")

    ds = xr.open_zarr(path, consolidated=True)
    have = set(ds.data_vars)
    if have != set(names):
        problems.append(f"variables: missing {sorted(set(names) - have)}, "
                        f"extra {sorted(have - set(names))}")
    nt = calendar.monthrange(year, month)[1] * 8
    if ds.sizes.get("time") != nt:
        problems.append(f"time steps {ds.sizes.get('time')} != {nt}")
    if (ds.sizes.get("latitude"), ds.sizes.get("longitude")) not in {grid, (grid[0], grid[1] + 1)}:
        problems.append(f"grid {ds.sizes.get('latitude')}x{ds.sizes.get('longitude')} != {grid}")
    t = ds["time"].values
    if len(t) and t[0] != np.datetime64(f"{year}-{month:02d}-01T00:00"):
        problems.append(f"first time {t[0]}")
    if len(t) > 1 and not np.all(np.diff(t) == np.timedelta64(3, "h")):
        problems.append("time axis not uniformly 3-hourly")

    short = []
    for v in sorted(have & set(names)):
        if ds[v].dtype != np.float32:
            problems.append(f"{v} dtype {ds[v].dtype}")
        present, expected = chunk_census(path, zmeta, v)
        if present < expected:
            short.append(f"{v} {present}/{expected}")
        mid = ds[v].isel(time=nt // 2).values
        fin = float(np.isfinite(mid).mean())
        if fin < 1.0:
            problems.append(f"{v} mid-month step only {fin:.4%} finite")
    if short:
        problems.append("chunk files short: " + ", ".join(short))
    if "f2d" in have:
        info["f2d_edge_finite"] = (float(np.isfinite(ds["f2d"].isel(time=0).values).mean()),
                                   float(np.isfinite(ds["f2d"].isel(time=-1).values).mean()))
    return problems, info


LOG_OUT = re.compile(r"^\s*output\s*:\s*(\S+\.zarr)")
LOG_VARS = re.compile(r"^\s*variables\s*:\s*(\d+) of 21")
LOG_FIX = re.compile(r"^\s*audit fixes\s+(\S+)")
LOG_EXIT = re.compile(r"exit[= ](\d+)")


def scan_logs(dirs: list[Path]) -> dict[str, list[dict]]:
    """'<parent dir>/<store name>' -> compute records, oldest log first.

    Only COMPUTES are collected: a skip writes nothing, so it cannot change a
    store. Keyed by parent directory as well as name because the baseline and
    audit series share store names."""
    recs: dict[str, list[dict]] = {}
    files = sorted({p for d in dirs if d.is_dir() for p in d.glob("*.out")},
                   key=lambda p: p.stat().st_mtime)
    for f in files:
        if f.stat().st_size > 5e6:
            continue
        text = f.read_text(errors="replace")
        if "variables :" not in text or ".zarr" not in text:
            continue
        lines = text.splitlines()
        fix = next((m.group(1) for l in lines if (m := LOG_FIX.match(l))), None)
        nvars = next((int(m.group(1)) for l in lines if (m := LOG_VARS.match(l))), None)
        exits = [int(m.group(1)) for l in lines if "===" in l and (m := LOG_EXIT.search(l))]
        bangs = [l.strip() for l in lines if l.strip().startswith("!!")]
        for l in lines:
            m = LOG_OUT.match(l)
            if m:
                p = Path(m.group(1))
                recs.setdefault(f"{p.parent.name}/{p.name}", []).append(dict(
                    log=f.name, fixes=fix, nvars=nvars,
                    exit=exits[-1] if exits else None, bangs=bangs))
    return recs


def judge_logs(recs: list[dict] | None, fixes: str) -> list[str]:
    """The LAST compute of a store is what is on disk. It must have run under
    the expected fix set, written 21 of 21, exited 0 and printed no '!!'."""
    if not recs:
        return ["no job log found that wrote this store"]
    last = recs[-1]
    out = []
    if last["fixes"] != fixes:
        out.append(f"last writer ({last['log']}) ran under {last['fixes']!r}")
    if last["nvars"] != 21:
        out.append(f"last writer ({last['log']}) wrote {last['nvars']} of 21")
    if last["exit"] not in (0, None):
        out.append(f"last writer ({last['log']}) exit {last['exit']}")
    if last["bangs"]:
        out.append(f"last writer ({last['log']}) printed: {last['bangs'][:3]}")
    return out


def check_thresholds(path: Path, pair_dir: str, names: list[str]) -> list[str]:
    payload = json.loads(path.read_text())
    prov, thr = payload["provenance"], payload["thresholds"]
    problems = []
    period = str(prov.get("period", ""))
    print(f"   file      : {path.name}  created {payload.get('created')}")
    print(f"   period    : {period!r}")
    if not period.startswith(pair_dir + " "):
        problems.append(f"period {period!r} does not name {pair_dir!r}")
    if set(thr) != set(names):
        problems.append(f"diagnostics differ: {sorted(set(names) ^ set(thr))}")
    for n in names:
        ladder = [thr.get(n, {}).get(s, float("nan")) for s in SEVERITIES]
        if not all(math.isfinite(v) for v in ladder):
            problems.append(f"{n}: non-finite threshold {ladder}")
        elif not all(b > a for a, b in zip(ladder, ladder[1:])):
            problems.append(f"{n}: ladder not strictly increasing {ladder}")
    sizes = prov.get("sample_sizes", {})
    uniq = sorted(set(sizes.values()))
    print(f"   n per diagnostic: {uniq}  (full leap year 2000 on 721x1440 = "
          f"{366 * 8 * 721 * 1440:,})")
    signs = set(payload.get("signs", {}).values())
    print(f"   signs     : {sorted(signs)}")
    if signs != {"+"}:
        problems.append(f"signs {signs}")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--na-subdir", default="derived/north_atlantic_audit")
    ap.add_argument("--global-subdir", default="derived/global_audit")
    ap.add_argument("--thresholds", default="calibration/thresholds_2026-09-11.json")
    ap.add_argument("--fixes", default=EXPECT_FIXES)
    ap.add_argument("--logs", nargs="*", default=None,
                    help="log directories (default: <repo>/logs and <base>/logs)")
    ap.add_argument("--no-logs", action="store_true")
    ap.add_argument("--csv", default=str(REPO / "cat_outputs" / "final" / "store_sweep_audit.csv"))
    args = ap.parse_args()

    t0 = time.time()
    base = Path(args.base)
    dmod = _load("diagnostics", REPO / "2_diagnostics.py")
    names = list(dmod.REFERENCE_TABLE)
    print(f">>> final sweep  {time.strftime('%Y-%m-%dT%H:%M:%S')}")
    print(f">>> base {base}\n>>> expected audit_fixes {args.fixes!r}, f2d A, DEF, 21 variables")

    logs = {} if args.no_logs else scan_logs(
        [Path(p) for p in (args.logs or [REPO / "logs", base / "logs"])])
    if not args.no_logs:
        print(f">>> logs: {sum(len(v) for v in logs.values())} store writes found, "
              f"covering {len(logs)} stores (all series)\n")

    total = 0
    rows = []
    for kind, sub, grid in (("na", args.na_subdir, (121, 301)),
                            ("global", args.global_subdir, (721, 1440))):
        d = base / sub
        exp = expected_stores(kind)
        on_disk = sorted(p.name for p in d.glob("*.zarr")) if d.is_dir() else []
        extra = sorted(set(on_disk) - {e[0] for e in exp})
        print("=" * 90)
        print(f"{sub}: expected {len(exp)}, on disk {len(on_disk)}"
              f"{', UNEXPECTED ' + str(extra) if extra else ''}")
        print("=" * 90)
        n_bad = 0
        edges = set()
        fixes_seen: dict[str, int] = {}
        for name, y, m, src in exp:
            try:
                probs, info = check_store(d / name, y, m, src, names, grid, args.fixes)
            except Exception as exc:                      # report, never crash the sweep
                probs, info = [f"could not be checked: {type(exc).__name__}: {exc}"], {}
            if not args.no_logs:
                probs += [f"log: {p}" for p in judge_logs(logs.get(f"{Path(sub).name}/{name}"), args.fixes)]
            fixes_seen[info.get("audit_fixes", "<none>")] = fixes_seen.get(
                info.get("audit_fixes", "<none>"), 0) + 1
            if "f2d_edge_finite" in info:
                edges.add(tuple(round(v, 4) for v in info["f2d_edge_finite"]))
            rows.append(dict(subdir=sub, store=name, n_problems=len(probs),
                             problems=" | ".join(probs)))
            if probs:
                n_bad += 1
                print(f"   {name}: " + " | ".join(probs))
        print(f"   audit_fixes values: {fixes_seen}")
        print(f"   f2d finite fraction at (first, last) step of the month: {sorted(edges)}")
        print(f"   stores with problems: {n_bad} of {len(exp)}\n")
        total += n_bad + len(extra)

    print("=" * 90)
    print("THRESHOLDS")
    print("=" * 90)
    tp = Path(args.thresholds)
    tp = tp if tp.is_absolute() else base / tp
    try:
        tprob = check_thresholds(tp, args.global_subdir, names)
    except Exception as exc:
        tprob = [f"could not be checked: {type(exc).__name__}: {exc}"]
    for p in tprob:
        print(f"   !! {p}")
    if not tprob:
        print("   21 ladders finite and strictly increasing; pairing names the audit reference year")
    total += len(tprob)

    Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.csv, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print(f"\nwrote {args.csv}")
    print(f"\nRESULT: {total} item(s) listed above  (wall {time.time() - t0:.0f} s)")
    return 0 if total == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
