#!/usr/bin/env python
"""
ada/calibration_check.py
========================
Stage 2 -> Stage 3: calibrate on the global year 2000, apply to North Atlantic
DJF 1979, and compare against Prosser et al. (2023) Table 1.

    pixi run python ada/calibration_check.py [--out-dir DIR] [--skip-split-half]

THE ACCEPTANCE CRITERION
------------------------
Prosser Table 1, DJF 1979, for an average point in his North Atlantic box:

    LOG  128.9 h    LMOG 45.6 h    MOG 22.3 h    MSOG 12.1 h    SOG 6.4 h

over 90 days = 2160 hours, i.e. 5.97 / 2.11 / 1.03 / 0.56 / 0.30 %.

Read it at three strengths:
  1. ORDER OF MAGNITUDE -- MOG near 1%, not 10% or 0.01%. Failing this means
     stop; something is wrong that no amount of interpretation fixes.
  2. THE LADDER SHAPE -- the five frequencies in the ratios
     5.97 : 2.11 : 1.03 : 0.56 : 0.30. This is the strong test. It constrains
     the whole tail shape of the North Atlantic distribution relative to the
     global one, using five numbers at once. No accidental agreement fakes it.
  3. ABSOLUTE agreement within ~30% would be a strong result given the
     sub-sampled calibration and the 175/200/225 vs 188/197/206 hPa stencil.

WHY THIS IS THE RIGHT CHECK, AND MAGNITUDE COMPARISONS ARE NOT
--------------------------------------------------------------
Every published per-diagnostic table (Williams & Joshi 2013 Table 1, Williams
2017 Table 2) comes from GFDL-CM2.1 at ~2 deg. We compute on ERA5 at 0.25 deg,
and Williams (2017) says directly beneath his own table that "the thresholds
are dependent on the grid resolution of the atmospheric model". So absolute
magnitudes CANNOT match and a magnitude comparison is a smell test, not a
criterion.

The exceedance frequency is different: it is percentile-calibrated on our own
data, so a resolution offset moves the data and the threshold by the same
factor and cancels exactly. That is why this number is comparable across
studies when the raw diagnostic values are not.

WHAT ELSE IT CHECKS, IN ORDER
-----------------------------
  A. thresholds from the global year 2000, cos(phi)-weighted (Prosser s2)
  B. SPLIT-HALF: days 1+17 vs 9+25 -- was 48 days enough sampling? Raw point
     counts overstate the effective sample size because gridpoints within a
     timestep are spatially correlated, so this measures it rather than
     assuming it.
  C. IDENTITY ON REAL DATA: apply the thresholds back to the calibration data.
     Must return 3.0 / 0.9 / 0.4 / 0.2 / 0.1% exactly, by construction. The
     synthetic version of this passes (tests/test_calibration_roundtrip.py);
     this additionally exercises real NaN patterns.
  D. Williams (2017) Table 2: SIGN and order of magnitude only.
  E. The Prosser comparison itself.
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
sys.path.insert(0, str(REPO))

import aggregate                                    # noqa: E402
import calibration                                  # noqa: E402
from calib_weighted_percentile import weighted_percentile  # noqa: E402


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, REPO / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Williams (2017) Table 1, adopted verbatim by Prosser (2023) s2.
SEVERITIES = {
    "light":              97.0,
    "light_to_moderate":  99.1,
    "moderate":           99.6,
    "moderate_to_severe": 99.8,
    "severe":             99.9,
}
ORDER = list(SEVERITIES)
LABEL = {"light": "LOG", "light_to_moderate": "LMOG", "moderate": "MOG",
         "moderate_to_severe": "MSOG", "severe": "SOG"}

# Prosser (2023) Table 1, DJF 1979, hours per season at an average point in
# his North Atlantic box.
PROSSER_DJF_1979_HOURS = {"light": 128.9, "light_to_moderate": 45.6,
                          "moderate": 22.3, "moderate_to_severe": 12.1,
                          "severe": 6.4}
DJF_HOURS = 90 * 24                      # Jan 31 + Feb 28 + Dec 31 = 90 days

# Prosser (2023), Figure 2 caption: "the North Atlantic (36-60N and 55-10W)".
# A strict subset of our 30-60N / 75-0W download, so this is exact, not
# approximate. Comparing on our larger box would give systematically different
# frequencies and look like an error.
PROSSER_BOX = dict(lat=(36.0, 60.0), lon=(-55.0, -10.0))

CALIB_DAYS = [1, 9, 17, 25]
SPLIT_HALVES = {"days 1+17": [1, 17], "days 9+25": [9, 25]}

# Williams (2017) Table 2 -- GFDL-CM2.1, so magnitudes are NOT comparable.
# Kept for SIGN and order-of-magnitude only. In the table's own units.
W2017_TABLE2 = {
    "negative_richardson":   [-15.4, -9.8, -7.9, -6.7, -5.9],
    "vertical_wind_shear":   [5.3, 6.6, 7.4, 7.9, 8.4],
    "colson_panofsky":       [-29.3, -27.0, -25.2, -23.7, -22.2],
    "f2d":                   [770, 1280, 1660, 1980, 2340],
    "brown1":                [99, 106, 110, 113, 118],
    "brown2":                [870, 1370, 1730, 2030, 2330],
    "ti1":                   [195, 292, 360, 419, 472],
    "ti2":                   [184, 282, 356, 419, 477],
    "deformation":           [50.9, 60.9, 66.9, 71.8, 76.3],
    "magnitude_pv":          [8.33, 8.73, 8.98, 9.19, 9.41],
    "vorticity_squared":     [2.46, 3.74, 4.70, 5.50, 6.24],
    "temperature_gradient":  [14.7, 17.6, 19.4, 20.8, 22.0],
    "wind_speed":            [40.9, 48.4, 52.4, 55.3, 58.5],
    "endlich":               [3.21, 3.94, 4.39, 4.72, 5.08],
    "ngm1":                  [1.65, 2.29, 2.76, 3.17, 3.54],
    "ngm2":                  [53, 84, 106, 127, 151],
    "ubf":                   [1230, 1840, 2270, 2610, 2960],
    "horizontal_divergence": [11.9, 15.7, 18.2, 20.4, 22.5],
    "ncsu1":                 [1200, 3600, 6300, 9300, 13000],
    "nva":                   [1.33, 1.86, 2.23, 2.56, 2.93],
    "rva_magnitude":         [1.44, 1.99, 2.34, 2.66, 3.00],
}
# Williams 2017 prints NVA/RVA in 10^-9 s^-2 while REFERENCE_TABLE (following
# W&J 2013) uses 10^-10. To express a value given in 10^-9 units in 10^-10
# units you MULTIPLY by ten: V x 10^-9 / 10^-10 = 10V.
#
# The first version of this divided, which made the published values ten times
# too small and inflated the printed ratios by a hundred -- nva read 1081 when
# it is 10.8, rva 1327 when it is 13.3. Those two then looked like wild
# outliers rather than sitting in the same resolution family as everything
# else. A reminder that a "sanity check" is only as sane as its own arithmetic.
W2017_UNIT_MULTIPLIER = {"nva": 10.0, "rva_magnitude": 10.0}


# ---------------------------------------------------------------------------
def open_months(paths: list[Path], days: list[int] | None = None) -> xr.Dataset:
    """Open zarr stores lazily and concatenate along time.

    Lazy is not an optimisation here, it is what makes this run at all: 12
    global months x 21 diagnostics is ~33 GB. calibration.compute_thresholds
    loops over diagnostics and materialises one at a time, so the peak is set
    by ONE diagnostic (~1.6 GB float32, ~20 GB once weighted_percentile casts
    to float64 and sorts), not by the whole set.
    """
    dss = [xr.open_zarr(p) for p in paths]
    ds = xr.concat(dss, dim="time") if len(dss) > 1 else dss[0]
    if days is not None:
        ds = ds.sel(time=ds["time.day"].isin(days))
    return ds


def subset_box(ds: xr.Dataset, lat: tuple, lon: tuple) -> xr.Dataset:
    """Subset by bounds regardless of coordinate order.

    ERA5 latitudes descend (90 -> -90), so a naive slice(36, 60) silently
    returns an EMPTY selection rather than an error. That would make the
    comparison quietly meaningless.
    """
    la = ds["latitude"].values
    lo = ds["longitude"].values
    lat_sl = slice(lat[1], lat[0]) if la[0] > la[-1] else slice(lat[0], lat[1])
    lon_sl = slice(lon[1], lon[0]) if lo[0] > lo[-1] else slice(lon[0], lon[1])
    out = ds.sel(latitude=lat_sl, longitude=lon_sl)
    if out.sizes["latitude"] == 0 or out.sizes["longitude"] == 0:
        raise ValueError(
            f"box {lat} x {lon} selected nothing from latitude "
            f"{la.min()}..{la.max()}, longitude {lo.min()}..{lo.max()}"
        )
    return out


def lat_weights_for(ds: xr.Dataset) -> xr.DataArray:
    return xr.DataArray(np.cos(np.deg2rad(ds["latitude"].values)),
                        coords={"latitude": ds["latitude"]}, dims=("latitude",),
                        name="lat_weights")


def progress_percentile(names: list[str], stage: str):
    """Wrap weighted_percentile so a long stage reports where it is.

    calibration.compute_thresholds takes weighted_percentile as an INJECTED
    argument, so progress reporting needs no change to the library.

    Sorting 21 arrays of 4e8 float64 values takes about an hour. Without this
    the log shows a stage header and then nothing, which is indistinguishable
    from a hang -- the same problem PYTHONUNBUFFERED solved at the job level,
    one layer down.

    TWO CALL SHAPES
    ---------------
    Since the 2026-08-29 fix, compute_thresholds asks for all five severities
    in ONE call, passing an array of percentiles. Older copies asked one
    severity at a time. This handles both rather than silently mislabelling
    progress if the two files ever drift:

      * array `pct`  -> one call per diagnostic; report on every call
      * scalar `pct` -> len(SEVERITIES) calls per diagnostic; report on the
                        first of each group, and memoise the other four onto
                        the same sort so the 5x cost is not paid either way
    """
    import time
    pcts = list(SEVERITIES.values())
    state = {"n": 0, "i": 0, "t0": time.time(), "cache": {}}

    def _report():
        state["i"] += 1
        i = state["i"]
        elapsed = time.time() - state["t0"]
        eta = elapsed / (i - 1) * (len(names) - i + 1) if i > 1 else 0.0
        print(f"   [{stage} {i:>2}/{len(names)}] {names[i - 1]:<23} "
              f"{elapsed / 60:>5.1f} min elapsed"
              + (f", ~{eta / 60:>4.1f} min left" if i > 1 else ""), flush=True)

    def wrapped(values, weights, pct):
        if np.ndim(pct) > 0:                       # one call, all severities
            result = weighted_percentile(values, weights, pct)
            _report()
            return result
        if state["n"] % len(pcts) == 0:            # legacy: first of a group
            allp = np.atleast_1d(weighted_percentile(
                values, weights, np.asarray(pcts, dtype=float)))
            state["cache"] = dict(zip(pcts, [float(v) for v in allp]))
            _report()
        state["n"] += 1
        return state["cache"][pct]

    return wrapped


def weighted_rate(exceed: xr.DataArray, weights: xr.DataArray) -> float:
    """cos(phi)-weighted mean of a 0/1/NaN exceedance field.

    NaN cells are excluded from numerator AND denominator, matching
    aggregate.exceedance_field's deliberate NaN propagation.
    """
    w = weights.broadcast_like(exceed)
    populated = exceed.notnull()
    return float((exceed.fillna(0.0) * w).sum() / w.where(populated).sum())


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--out-dir", default=None,
                    help="where thresholds_<date>.json goes "
                         "(default <base>/calibration)")
    ap.add_argument("--skip-split-half", action="store_true",
                    help="skip check B; it roughly doubles the runtime")
    args = ap.parse_args()

    base = Path(args.base)
    out_dir = Path(args.out_dir) if args.out_dir else base / "calibration"

    diag = _load("diagnostics", "2_diagnostics.py")
    pipeline = _load("pipeline", "3_pipeline.py")
    reference_table = diag.REFERENCE_TABLE

    glob_paths = sorted((base / "derived/global").glob("diagnostics_glob_2000-*.zarr"))
    na_paths = [base / "derived/north_atlantic" / f"diagnostics_na_{m}.zarr"
                for m in ("1979-01", "1979-02", "1979-12")]

    missing = [p for p in glob_paths + na_paths if not p.exists()]
    if missing or len(glob_paths) != 12:
        print(f"!! expected 12 global + 3 NA stores; found "
              f"{len(glob_paths)} global, missing: {[p.name for p in missing]}")
        return 1

    # =====================================================================
    print("=" * 78)
    print("A. CALIBRATION -- global year 2000, cos(phi)-weighted percentiles")
    print("=" * 78)
    calib = open_months(glob_paths)
    print(f"   {len(glob_paths)} months, dims {dict(calib.sizes)}")
    n_points = int(np.prod([calib.sizes[d] for d in ("latitude", "longitude", "time")]))
    print(f"   {n_points:,} points per diagnostic "
          f"({n_points * 0.004:,.0f} above the MOG threshold)")

    fields = {k: calib[k] for k in reference_table if k in calib}
    if len(fields) != 21:
        print(f"!! only {len(fields)} of 21 diagnostics present in the stores")
        return 1

    print("   sorting 21 arrays of 4e8 values; progress below\n")
    thresholds, signs, sample_sizes = calibration.compute_thresholds(
        calibration_fields=fields,
        lat_weights=lat_weights_for(calib),
        severities=SEVERITIES,
        reference_table=reference_table,
        weighted_percentile=progress_percentile(list(fields), "A"),
    )

    print(f"\n{'diagnostic':<23}" + "".join(f"{LABEL[s]:>13}" for s in ORDER))
    for name in sorted(thresholds, key=lambda k: reference_table[k]["num"]):
        print(f"{name:<23}" + "".join(f"{thresholds[name][s]:>13.4g}" for s in ORDER))

    out_dir.mkdir(parents=True, exist_ok=True)
    stamp = dt.date.today().isoformat()
    path = calibration.save_thresholds(
        out_dir / f"thresholds_{stamp}.json",
        thresholds, signs, SEVERITIES,
        domain="global",
        period="2000, days 1/9/17/25 of each month, 3-hourly",
        pressure_levels=[175, 200, 225],
        sample_sizes=sample_sizes,
        notes="Sub-sampled global year 2000 (48 days). Evaluated at 200 hPa "
              "from a 175/225 stencil; Prosser uses 197 hPa from 188/206. "
              "ERA5, not ERA5.1 (Prosser used ERA5.1 for 2000-2006).",
    )
    print(f"\n   saved -> {path}")

    # =====================================================================
    if not args.skip_split_half:
        print("\n" + "=" * 78)
        print("B. SPLIT-HALF -- was 48 days enough sampling?")
        print("=" * 78)
        halves = {}
        for label, days in SPLIT_HALVES.items():
            sub = open_months(glob_paths, days=days)
            print(f"   {label}: {sub.sizes['time']} timesteps")
            halves[label], _, _ = calibration.compute_thresholds(
                calibration_fields={k: sub[k] for k in fields},
                lat_weights=lat_weights_for(sub),
                severities=SEVERITIES,
                reference_table=reference_table,
                weighted_percentile=progress_percentile(
                    list(fields), f"B {label}"),
            )

        a, b = halves.values()
        rows = []
        for name in fields:
            for sev in ORDER:
                va, vb = a[name][sev], b[name][sev]
                denom = max(abs(va), abs(vb), 1e-300)
                rows.append({"diagnostic": name, "severity": sev,
                             "rel_diff": abs(va - vb) / denom})
        df = pd.DataFrame(rows)
        worst = df.loc[df["rel_diff"].idxmax()]
        print(f"\n   median disagreement between halves : {df['rel_diff'].median():.2%}")
        print(f"   90th percentile                    : {df['rel_diff'].quantile(0.9):.2%}")
        print(f"   worst                              : {worst['rel_diff']:.2%} "
              f"({worst['diagnostic']} @ {worst['severity']})")
        if df["rel_diff"].quantile(0.9) < 0.05:
            print("   -> 48 days is ENOUGH. Sub-sampling is not a limitation here.")
        elif df["rel_diff"].quantile(0.9) < 0.15:
            print("   -> marginal. Usable for this check; add days 5/13/21/29 "
                  "before any production calibration.")
        else:
            print("   -> NOT ENOUGH. Re-download with more days per month.")

    # =====================================================================
    print("\n" + "=" * 78)
    print("C. IDENTITY ON REAL DATA -- thresholds applied back to year 2000")
    print("=" * 78)
    print("   by construction the p-th percentile must be exceeded by (100-p)%")
    cw = lat_weights_for(calib)
    # Diagnostic OUTER, severity INNER, and materialise once: the naive
    # ordering re-reads each diagnostic from zarr five times, turning 21 reads
    # of 1.6 GB into 105.
    rates: dict[str, list[float]] = {sev: [] for sev in ORDER}
    import time as _time
    _t0 = _time.time()
    for i, n in enumerate(fields, 1):
        vals = fields[n].compute()
        for sev in ORDER:
            rates[sev].append(weighted_rate(
                aggregate.exceedance_field(vals, thresholds[n][sev], signs[n]), cw))
        del vals
        el = _time.time() - _t0
        print(f"   [C {i:>2}/{len(fields)}] {n:<23} {el / 60:>5.1f} min elapsed, "
              f"~{el / i * (len(fields) - i) / 60:>4.1f} min left", flush=True)

    print(f"\n   {'severity':<20}{'expected':>10}{'observed':>10}{'rel err':>10}")
    identity_worst = 0.0
    for sev in ORDER:
        expected = (100.0 - SEVERITIES[sev]) / 100.0
        observed = float(np.mean(rates[sev]))
        rel = abs(observed - expected) / expected
        identity_worst = max(identity_worst, rel)
        flag = "" if rel < 0.02 else "   <-- OFF"
        print(f"   {LABEL[sev]:<20}{expected:>9.4%}{observed:>10.4%}{rel:>10.2%}{flag}")
    print(f"\n   {'PASS' if identity_worst < 0.02 else 'FAIL'} "
          f"(worst {identity_worst:.2%})")

    # =====================================================================
    print("\n" + "=" * 78)
    print("D. vs WILLIAMS (2017) TABLE 2 -- SIGN and ORDER only")
    print("=" * 78)
    print("   GFDL-CM2.1 at ~2 deg vs ERA5 at 0.25 deg. Williams states the")
    print("   thresholds are resolution-dependent, so magnitudes CANNOT match.")
    print("   A sign mismatch, however, would be a real defect.\n")
    print(f"   {'diagnostic':<23}{'ours (MOG)':>14}{'W2017 (MOG)':>14}"
          f"{'ratio':>9}  sign")
    sign_problems = []
    for name in sorted(fields, key=lambda k: reference_table[k]["num"]):
        if name not in W2017_TABLE2:
            continue
        native = thresholds[name]["moderate"]
        val = native
        if name in pipeline.PRETRANSFORM:
            val = float(pipeline.PRETRANSFORM[name](np.asarray(val)))
        ours = val * pipeline.SCALE_TO_TABLE.get(name, 1.0)
        theirs = W2017_TABLE2[name][2] * W2017_UNIT_MULTIPLIER.get(name, 1.0)
        ok = np.sign(ours) == np.sign(theirs)
        if not ok:
            sign_problems.append(name)
        print(f"   {name:<23}{ours:>14.4g}{theirs:>14.4g}"
              f"{ours / theirs:>9.2f}  {'ok' if ok else 'MISMATCH'}")
    print(f"\n   {len(sign_problems)} sign mismatch(es)"
          + (f": {sign_problems}" if sign_problems else ""))

    # =====================================================================
    print("\n" + "=" * 78)
    print("E. THE CHECK -- North Atlantic DJF 1979 vs Prosser (2023) Table 1")
    print("=" * 78)
    na = open_months(na_paths)
    print(f"   loaded {na.sizes['time']} timesteps "
          f"(expect 720 = 248+224+248, minus F2D's 2 per month)")
    na_box = subset_box(na, **PROSSER_BOX)
    print(f"   Prosser box 36-60N / 55-10W -> {na_box.sizes['latitude']} x "
          f"{na_box.sizes['longitude']} gridpoints")

    na_fields = {k: na_box[k] for k in fields}
    per_sev = aggregate.exceedance_mean_all_severities(
        na_fields, thresholds, ORDER, signs)
    nw = lat_weights_for(na_box)

    print(f"\n   {'':<8}{'ours':>10}{'Prosser':>10}{'ratio':>9}"
          f"{'ours (h)':>11}{'Prosser (h)':>13}")
    results = {}
    for sev in ORDER:
        obs = weighted_rate(per_sev[sev]["exceedance_mean"], nw)
        target = PROSSER_DJF_1979_HOURS[sev] / DJF_HOURS
        results[sev] = obs
        print(f"   {LABEL[sev]:<8}{obs:>9.3%}{target:>10.3%}"
              f"{obs / target:>9.2f}{obs * DJF_HOURS:>11.1f}"
              f"{PROSSER_DJF_1979_HOURS[sev]:>13.1f}")

    mog = results["moderate"] / (PROSSER_DJF_1979_HOURS["moderate"] / DJF_HOURS)
    ours = np.array([results[s] for s in ORDER])
    theirs = np.array([PROSSER_DJF_1979_HOURS[s] / DJF_HOURS for s in ORDER])
    # Ladder shape: normalise both by MOG and compare the remaining four.
    shape_err = float(np.max(np.abs((ours / ours[2]) / (theirs / theirs[2]) - 1)))

    print("\n" + "-" * 78)
    print(f"   1. ORDER OF MAGNITUDE   MOG ratio {mog:.2f}   "
          f"{'PASS' if 0.2 <= mog <= 5 else 'FAIL -- STOP'}")
    print(f"   2. LADDER SHAPE         worst deviation {shape_err:.1%}   "
          f"{'PASS' if shape_err < 0.35 else 'CHECK'}")
    print(f"   3. ABSOLUTE (<30%)      {abs(mog - 1):.1%} from Prosser   "
          f"{'PASS' if abs(mog - 1) < 0.30 else 'not met (not required)'}")
    print("-" * 78)
    print("\n   Criterion 2 is the strong one: it constrains the whole tail")
    print("   shape with five numbers at once, and it is what no accidental")
    print("   agreement can fake.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
