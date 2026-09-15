"""
thresholds.py -- Stage 3: severity thresholds from the global year 2000.

What is computed
----------------
For each diagnostic D and each severity s with percentile q_s (config.SEVERITIES:
97, 99.1, 99.6, 99.8, 99.9), the threshold is the q_s-th percentile of D over
every grid cell and every 3-hourly time step of the global year 2000, with each
cell weighted by the area it represents, w = cos(latitude):

    threshold_s = Q_w(q_s)   where   Q_w = weighted quantile function of D.

The weighted quantile is the Hazen-type one. Sort the values, give value i the
weight w_i, let C_i be the cumulative weight up to and including i and W the
total weight. Value i sits at the plotting position

    P_i = (C_i - w_i / 2) / W

and Q_w(q) is found by linear interpolation of the values against P_i (equal
values are merged first). With equal weights this is exactly the classic
(i - 1/2)/n Hazen rule. `weighted_percentile` below does this on a whole array.

Why a second, streaming version
-------------------------------
A full global year is 721 x 1440 x 2928 = 3.0e9 values per diagnostic: far too
much to sort in memory. But the formula above only needs the values near the
top. Everything below a cut enters P_i through ONE number, the total weight
below the cut, W_below:

    P_i = (W_below + C_i^tail - w_i / 2) / W      for every value i above the cut.

So one pass over the year that keeps the values above a cut (about the top 12 %),
and adds up W and W_below, gives exactly the same thresholds as sorting
everything (tests/test_operators.py checks this to 1e-9). The cut is placed at
p88 of a thinned sample, well below p97; a guard stops the run if the kept tail
does not reach far enough down, because then np.interp would silently return
the cut itself as the "threshold".

Output: a JSON file in the same format as the original calibration.py (so both
code bases can read each other's files), plus, optionally, the kept tails as
.npz files (values sorted ascending, with their latitude row) -- the
peaks-over-threshold archive for extreme-value work.
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import numpy as np
import xarray as xr

import config as C


# ---------------------------------------------------------------------------
# The weighted Hazen percentile on a whole array (small data, tests, notebook)
# ---------------------------------------------------------------------------
def weighted_percentile(values, weights, percentiles):
    """Weighted Hazen percentiles (percentiles in 0..100). NaNs are dropped."""
    values = np.asarray(values, dtype=np.float64).ravel()
    weights = np.asarray(weights, dtype=np.float64).ravel()
    keep = np.isfinite(values)
    values, weights = values[keep], weights[keep]
    order = np.argsort(values, kind="stable")
    v_unique, inverse = np.unique(values[order], return_inverse=True)
    w_unique = np.zeros(v_unique.shape)
    np.add.at(w_unique, inverse, weights[order])
    position = (np.cumsum(w_unique) - 0.5 * w_unique) / weights.sum()
    q = np.atleast_1d(np.asarray(percentiles, dtype=np.float64) / 100.0)
    return np.interp(q, position, v_unique, left=v_unique[0], right=v_unique[-1])


def latitude_weights(latitude) -> np.ndarray:
    """w = cos(phi), float64."""
    return np.cos(np.deg2rad(np.asarray(latitude, dtype=np.float64)))


# ---------------------------------------------------------------------------
# The streaming (tail) version for the full year
# ---------------------------------------------------------------------------
def time_blocks(da: xr.DataArray, block: int):
    """Yield the field as (lat, lon, time) float32 numpy blocks of `block` steps."""
    da = da.transpose("latitude", "longitude", "time")
    n = da.sizes["time"]
    for a in range(0, n, block):
        yield np.asarray(da.isel(time=slice(a, min(a + block, n))).values, dtype=np.float32)


def place_cut(name, stores, w_lat, stride, cut_percentile) -> float:
    """p88 (by default) of every `stride`-th time step: where to cut the tail.

    The stride must not be a multiple of 8, or the sample would contain a
    single time of day (3-hourly data has 8 steps per day).
    """
    if stride % C.STEPS_PER_DAY == 0:
        raise ValueError(f"stride {stride} samples one time of day; use e.g. 7")
    values, weights = [], []
    for store in stores:
        block = np.asarray(xr.open_zarr(store)[name].isel(time=slice(None, None, stride))
                           .transpose("latitude", "longitude", "time").values, dtype=np.float32)
        finite = np.isfinite(block)
        values.append(block[finite].astype(np.float64))
        weights.append(np.broadcast_to(w_lat[:, None, None], block.shape)[finite])
    return float(weighted_percentile(np.concatenate(values), np.concatenate(weights),
                                     [cut_percentile])[0])


def scan_tail(name, stores, w_lat, cut, block):
    """One pass over the whole year: keep values >= cut, add up W and W_below."""
    kept_values, kept_rows = [], []
    w_total = w_below = 0.0
    n_finite = 0
    rows = np.arange(w_lat.size, dtype=np.int16)
    for store in stores:
        for x in time_blocks(xr.open_zarr(store)[name], block):
            finite = np.isfinite(x)
            keep = finite & (x >= cut)                    # '>=': a tie is never split
            n_fin_lat = finite.sum(axis=(1, 2))           # per latitude row
            n_keep_lat = keep.sum(axis=(1, 2))
            w_total += float(np.dot(w_lat, n_fin_lat))
            w_below += float(np.dot(w_lat, n_fin_lat - n_keep_lat))
            n_finite += int(n_fin_lat.sum())
            if n_keep_lat.sum():
                kept_values.append(x[keep])
                kept_rows.append(np.broadcast_to(rows[:, None, None], x.shape)[keep])
    if not kept_values:
        raise ValueError(f"{name}: nothing above the cut {cut}")
    return np.concatenate(kept_values), np.concatenate(kept_rows), w_below, w_total, n_finite


def thresholds_from_tail(values, rows, w_lat, w_below, w_total, percentiles, name, margin):
    """The weighted Hazen percentiles, from the kept tail plus W_below and W."""
    v = values.astype(np.float64)
    w = w_lat[rows]
    order = np.argsort(v, kind="stable")
    v, w = v[order], w[order]
    v_unique, inverse = np.unique(v, return_inverse=True)
    w_unique = np.zeros(v_unique.shape)
    np.add.at(w_unique, inverse, w)
    position = (w_below + np.cumsum(w_unique) - 0.5 * w_unique) / w_total

    q = np.asarray(percentiles, dtype=np.float64) / 100.0
    if position[0] > q.min() - margin:                   # THE GUARD
        raise SystemExit(f"GUARD FAILED for {name}: the kept tail starts at plotting position "
                         f"{position[0]:.5f}, not below {q.min() - margin:.5f}. "
                         f"Lower the cut percentile and re-run.")
    out = np.interp(q, position, v_unique, left=v_unique[0], right=v_unique[-1])
    return out, v.astype(np.float32), rows[order], float(position[0])


def calibrate(derived: Path, out_json: Path, convention_name: str,
              archive_dir: Path | None = None, cut_percentile=C.TAIL_CUT_PERCENTILE,
              stride=C.TAIL_CUT_STRIDE, block=C.TAIL_TIME_CHUNK) -> dict:
    """Thresholds for all 21 diagnostics from the 12 global stores in `derived`."""
    stores = sorted(Path(derived).glob("diagnostics_glob_*.zarr"))
    if len(stores) != 12:
        raise SystemExit(f"expected 12 monthly global stores in {derived}, found {len(stores)}")
    check_store_conventions(stores, convention_name)

    first = xr.open_zarr(stores[0])
    names = [k for k in C.DIAGNOSTIC_KEYS if k in first.data_vars]
    w_lat = latitude_weights(first["latitude"].values)
    percentiles = np.asarray(list(C.SEVERITIES.values()), dtype=float)
    if archive_dir:
        Path(archive_dir).mkdir(parents=True, exist_ok=True)

    thresholds, sample_sizes = {}, {}
    for i, name in enumerate(names, 1):
        cut = place_cut(name, stores, w_lat, stride, cut_percentile)
        values, rows, w_below, w_total, n = scan_tail(name, stores, w_lat, cut, block)
        out, v_sorted, rows_sorted, p_first = thresholds_from_tail(
            values, rows, w_lat, w_below, w_total, percentiles, name, C.TAIL_GUARD_MARGIN)
        thresholds[name] = {s: float(x) for s, x in zip(C.SEVERITIES, out)}
        sample_sizes[name] = n
        if archive_dir:
            np.savez(Path(archive_dir) / f"{name}.npz", values=v_sorted, lat_index=rows_sorted,
                     w_below=np.float64(w_below), w_total=np.float64(w_total),
                     cut=np.float64(cut), n_finite=np.int64(n))
        print(f"  [{i:2d}/21] {name:24s} cut={cut: .6g}  kept={1 - w_below / w_total:6.2%}  "
              f"n={n:,}  first position={p_first:.5f}", flush=True)

    save_thresholds(out_json, thresholds, sample_sizes, convention_name,
                    period=f"{derived} ({len(stores)} stores)",
                    notes=f"essentials/thresholds.py, tail cut p{cut_percentile}, stride {stride}")
    return thresholds


def check_store_conventions(stores, convention_name):
    """Refuse to calibrate on stores written under another convention."""
    for store in stores:
        attrs = json.loads((Path(store) / ".zattrs").read_text())
        have = attrs.get("convention")
        if have is None:            # a store written by the original pipeline
            have = original_store_convention(attrs)
        if have != convention_name:
            raise SystemExit(f"{Path(store).name} was built under '{have}', not '{convention_name}'")


def original_store_convention(attrs: dict) -> str:
    """Name the convention of a store written by the ORIGINAL pipeline."""
    fixes = attrs.get("audit_fixes", "baseline")
    variant = attrs.get("f2d_variant")
    for name, conv in C.CONVENTIONS.items():
        want_fixes = "+".join(k for k, on in (("meridional_metric", conv["dy"] == "map"),
                                              ("endlich_component_shear", conv["endlich"] == "components"),
                                              ("four_variable_provenance", conv["fields"] == "computed")) if on)
        if set(fixes.split("+")) == set((want_fixes or "baseline").split("+")) and variant == conv["f2d_variant"]:
            return name
    return f"unknown(audit_fixes={fixes}, f2d_variant={variant})"


# ---------------------------------------------------------------------------
# File format (identical schema to the original calibration.py, version 1)
# ---------------------------------------------------------------------------
def save_thresholds(path, thresholds, sample_sizes, convention_name, period, notes=""):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "created": dt.datetime.now().isoformat(timespec="seconds"),
        "provenance": {
            "calibration_domain": "global",
            "period": period,
            "pressure_levels": [str(p) for p in C.PRESSURE_LEVELS_HPA],
            "severities_percentiles": dict(C.SEVERITIES),
            "rojak_version": None,
            "sample_sizes": sample_sizes,
            "notes": notes,
            "convention": convention_name,
            "convention_settings": C.CONVENTIONS[convention_name],
        },
        "signs": {k: C.SIGN[k] for k in thresholds},
        "thresholds": thresholds,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(f"    saved -> {path}")
    return path


def load_thresholds(path):
    """Return (thresholds, convention_name or None, provenance)."""
    payload = json.loads(Path(path).read_text())
    if payload.get("schema_version") != 1:
        raise ValueError(f"{path}: unknown schema_version {payload.get('schema_version')}")
    prov = payload["provenance"]
    return payload["thresholds"], prov.get("convention"), prov


def compare_threshold_files(a, b, tolerance=1e-6) -> bool:
    """Largest relative difference between two threshold files."""
    ta, _, _ = load_thresholds(a)
    tb, _, _ = load_thresholds(b)
    worst, where = 0.0, ""
    for name in C.DIAGNOSTIC_KEYS:
        for sev in C.SEVERITIES:
            x, y = ta[name][sev], tb[name][sev]
            rel = abs(x - y) / max(abs(y), 1e-300)
            if rel > worst:
                worst, where = rel, f"{name}/{sev}"
    print(f"largest relative difference {worst:.3e} at {where} (tolerance {tolerance:.0e})")
    return worst <= tolerance
