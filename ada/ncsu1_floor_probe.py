"""
ada/ncsu1_floor_probe.py
========================
Audit fix F5. Measure, do not change, the Richardson floor inside NCSU1.

WHY
---
`ncsu1` (#19 in Williams 2017 order) is the largest single disagreement with
Prosser in the whole replication: a 1979 exceedance level ratio of **6.07**, and
the only one of the 21 whose error has the opposite sign to everything else.
AUDIT_diagnostics_vs_literature.md §5.6 lists three candidate causes and cannot
separate them from the published figures alone. This measures the first one.

Sharman (2006) A36, p. 284 -- reproduced verbatim by Williams & Storer (2022)
Eq. (7), p. 1428, so it is definitively the form being replicated -- reads

    NCSU1 = [1 / MAX(Ri, 1e-5)] * MAX(u du/dx + v dv/dy, 0) * |grad zeta|

`MAX(Ri, 1e-5)` is a LOWER bound, not `|Ri|`. Wherever `Ri <= 1e-5` -- i.e.
essentially wherever `N^2 <= 0`, since `Ri = N^2/Sv^2` and `Sv^2 > 0` -- the
floor binds and `1/Ri` returns **1e5**, against roughly 1e-2 for a typical
`Ri ~ 100` at 200 hPa. That is a seven-order-of-magnitude amplification, applied
to a set of cells whose geography is decided by the sign of the static stability
on a 50 hPa stencil.

If that set is large enough to reach the 97th percentile, then NCSU1's entire
turbulent tail is a map of `N^2 <= 0` rather than of anything the other twenty
diagnostics respond to -- and since Prosser resolves `N^2` on an 18 hPa stencil,
his map of `N^2 <= 0` is a different map. That would explain the 6.07 and it
would explain the sign.

THIS IS A MEASUREMENT, NOT A FIX. A36 is transcribed correctly and must not be
"corrected": the floor is in the published formula, twice confirmed. The output
here decides which of §5.6's hypotheses survives, nothing more.

HOW TO READ THE OUTPUT
----------------------
The decisive number is the LAST one: the cos(phi)-weighted fraction of NCSU1's
own p97 exceedance set that sits on the floor.

  >= ~50%   The tail IS the floor. §5.6 candidate 1 is confirmed and the 6.07 is
            a stencil effect operating through the sign of N^2, not through its
            magnitude. Report it; do not change A36.
  ~3-50%    The floor contributes but does not dominate. Both candidate 1 and
            candidate 2 (|grad zeta| from ERA5's archived vorticity) are live.
  << 1%     Candidate 1 is dead. The floor never binds hard enough to matter and
            the 6.07 must come from |grad zeta| or from the MAX(.,0) clip.

USAGE
-----
    python ada/ncsu1_floor_probe.py <in.grib|in.nc> [--target-level 200]
                                    [--out data_prosser/ab/]

Run it on a GLOBAL file. The thresholds NCSU1 is calibrated against are global
percentiles, so the fraction that matters is the global one; a North Atlantic
file cannot tell you whether the floor reaches the global p97.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from calib_weighted_percentile import weighted_percentile  # noqa: E402

RI_FLOOR = 1e-5          # Sharman A36, p. 284, verbatim
LIGHT_PCT = 97.0         # Williams (2017) Table 1, p. 580
LAT_BAND = 10.0


def _load(module_name: str, filename: str):
    spec = importlib.util.spec_from_file_location(module_name, REPO / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def run(path: Path, target_level: int, out_dir: Path | None,
        max_timesteps: int | None = None) -> dict:
    diag = _load("diagnostics", "2_diagnostics.py")

    print(f"input        {path}")
    print(f"target level {target_level} hPa")
    print(f"floor        MAX(Ri, {RI_FLOOR:g})  -- Sharman A36 p. 284 / W&S Eq. 7 p. 1428\n")

    ds_raw = diag.load_era5(path)
    for _t in ("time", "valid_time"):
        if _t in ds_raw.dims and max_timesteps and ds_raw.sizes[_t] > max_timesteps:
            print(f"  capping {_t}: {ds_raw.sizes[_t]} -> {max_timesteps} timesteps")
            ds_raw = ds_raw.isel({_t: slice(0, max_timesteps)})
            break
    catdata = diag.prepare_for_rojak(ds_raw)
    ds = catdata._dataset

    ri = diag.richardson(ds, negative=False).sel(pressure_level=target_level)
    ncsu1 = diag.ncsu1(ds, target_level=target_level)

    lat = ri["latitude"]
    w = np.cos(np.deg2rad(lat)).clip(min=0.0).broadcast_like(ri)

    ri_v = np.asarray(ri.values, dtype=np.float64).ravel()
    n_v = np.asarray(ncsu1.values, dtype=np.float64).ravel()
    w_v = np.asarray(w.values, dtype=np.float64).ravel()
    lat_v = np.asarray(lat.broadcast_like(ri).values, dtype=np.float64).ravel()

    finite = np.isfinite(ri_v) & np.isfinite(n_v)
    floored = finite & (ri_v <= RI_FLOOR)
    negative_n2 = finite & (ri_v < 0.0)

    wt = w_v[finite].sum()
    frac_floored = float(w_v[floored].sum() / wt)
    frac_negative = float(w_v[negative_n2].sum() / wt)

    thr = float(weighted_percentile(n_v, w_v, LIGHT_PCT))
    tail = finite & (n_v >= thr)
    w_tail = float(w_v[tail].sum())
    frac_tail_floored = float(w_v[tail & floored].sum() / w_tail) if w_tail > 0 else float("nan")

    print(f"{'quantity':<52}{'weighted fraction':>20}")
    print("-" * 72)
    print(f"{'cells where the floor binds (Ri <= 1e-5)':<52}{100 * frac_floored:>19.4f}%")
    print(f"{'  of which convectively unstable (Ri < 0, i.e. N2 < 0)':<52}{100 * frac_negative:>19.4f}%")
    print(f"{'light-CAT tail of NCSU1 (p97), by construction':<52}{100 * (1 - LIGHT_PCT / 100):>19.4f}%")
    print()
    print(f"  NCSU1 p97 threshold        {thr:.6e} s-3")
    print(f"  median Ri                  {np.median(ri_v[finite]):.4e}")
    print(f"  1st percentile of Ri       {np.percentile(ri_v[finite], 1):.4e}")
    print()
    print(f"*** DECISIVE: {100 * frac_tail_floored:.2f} % of NCSU1's p97 exceedance set "
          f"sits on the floor ***")
    if frac_tail_floored >= 0.50:
        verdict = ("THE TAIL IS THE FLOOR. Audit 5.6 candidate 1 confirmed: NCSU1's "
                   "turbulent tail is a map of N2 <= 0 on a 50 hPa stencil, and "
                   "Prosser's 18 hPa stencil maps a different set. Explains both the "
                   "6.07 and its sign. Report it; do NOT change A36.")
    elif frac_tail_floored >= 0.03:
        verdict = ("THE FLOOR CONTRIBUTES BUT DOES NOT DOMINATE. Candidates 1 and 2 "
                   "(|grad zeta| from archived vorticity) are both still live; test "
                   "candidate 2 by recomputing zeta from vector_derivatives.")
    else:
        verdict = ("CANDIDATE 1 IS DEAD. The floor never reaches the tail. The 6.07 "
                   "must come from |grad zeta| or from the MAX(.,0) clip -- audit 5.6 "
                   "candidates 2 and 3.")
    print(f"\n{verdict}\n")

    edges = np.arange(-90.0, 90.0 + LAT_BAND, LAT_BAND)
    print("where the floor binds, by latitude band (weighted % of the band):")
    bands = {}
    for lo, hi in zip(edges[:-1], edges[1:]):
        band = finite & (lat_v >= lo) & (lat_v < hi)
        bw = float(w_v[band].sum())
        if bw <= 0:
            continue
        f = float(w_v[band & floored].sum() / bw)
        bands[f"{lo:+.0f}..{hi:+.0f}"] = 100 * f
        bar = "#" * int(round(min(f, 0.2) * 250))
        print(f"   {lo:+4.0f}..{hi:+4.0f}   {100 * f:7.3f}%  {bar}")

    payload = {
        "input": str(path), "target_level": target_level, "ri_floor": RI_FLOOR,
        "frac_floored": frac_floored, "frac_negative_n2": frac_negative,
        "ncsu1_p97_threshold": thr,
        "frac_of_p97_tail_floored": frac_tail_floored,
        "verdict": verdict, "floor_by_lat_band_pct": bands,
    }
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        f = out_dir / f"ncsu1_floor_{path.stem}.json"
        f.write_text(json.dumps(payload, indent=2))
        print(f"\nwritten: {f}")
    return payload


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path)
    p.add_argument("--target-level", type=int, default=200)
    p.add_argument("--max-timesteps", type=int, default=None,
                   help="use only the first N timesteps; needed to keep a "
                        "global month inside a modest memory request")
    p.add_argument("--out", type=Path, default=None)
    a = p.parse_args()
    run(a.input, a.target_level, a.out, a.max_timesteps)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
