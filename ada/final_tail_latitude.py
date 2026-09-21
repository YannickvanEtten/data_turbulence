#!/usr/bin/env python
"""
ada/final_tail_latitude.py
==========================
(FINAL-CHECKS-2026-09-16)
WHERE each diagnostic's global tail lives, read from the tail archive that the
reference-year calibration already wrote. No diagnostic data is opened.

    pixi run python ada/final_tail_latitude.py

WHY
    The global maps (Prosser Figs 1, 2, S5) were never computed and will not
    be: a 42-year global run is out of scope. But the question S5 answers for
    a single diagnostic -- is its tail where the literature says it is? --
    does not need 42 years. calibration/tails_2026-09-11/<name>.npz holds
    every retained reference-year value above the p88 cut WITH ITS LATITUDE
    ROW (ada/tail_thresholds.py). A zonal profile of the LOG, MOG and SOG
    exceedance sets therefore costs one read of the archive.

    The comparison target is qualitative and printed in Williams & Storer
    (2022, QJRMS 148, p. 1428-1429), same lineage as Prosser, 200 hPa, light
    CAT, ERA-Interim and HadGEM2-ES:
      - negative Richardson number and Colson-Panofsky: a single TROPICAL peak
      - Ellrod TI1: two MIDLATITUDE peaks
      - residual of the nonlinear balance equation: STRONGER IN THE NORTHERN
        HEMISPHERE in both seasons
      - frontogenesis function: no clear midlatitude peaks

    That last-but-one line is the one this exists for: our `ubf` puts 0.35x
    (MOG) and 0.10x (SOG) of an even share of its tail into Prosser's North
    Atlantic box, where Prosser has 1.8x and 1.35x. If our global UBF tail is
    not Northern-Hemisphere weighted either, the difference is structural (an
    implementation difference from the Williams lineage), not regional.

OUTPUT
    For every diagnostic and severity: the tail's DENSITY in six latitude
    bands, i.e. (band share of the cos-weighted exceedance) / (band share of
    the cos-weighted area). 1.00 everywhere = an evenly spread tail. Plus the
    NH/SH mass ratio. And the achieved global fraction, which must sit at 3 /
    0.4 / 0.1 % (the identity check, recomputed from the archive).
    Written to cat_outputs/final/tail_latitude_2000_audit.csv. No verdicts.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
SEVS = {"light": 3.0, "moderate": 0.4, "severe": 0.1}
BANDS = ["90S-60S", "60S-30S", "30S-0", "0-30N", "30N-60N", "60N-90N"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--archive", default="calibration/tails_2026-09-11")
    ap.add_argument("--thresholds", default="calibration/thresholds_2026-09-11.json")
    ap.add_argument("--latitude-from", default="derived/global_audit/diagnostics_glob_2000-01.zarr")
    ap.add_argument("--csv", default=str(REPO / "cat_outputs" / "final" / "tail_latitude_2000_audit.csv"))
    args = ap.parse_args()

    base = Path(args.base)
    arch = base / args.archive
    thr = json.loads((base / args.thresholds).read_text())
    period = thr["provenance"]["period"]
    thresholds = thr["thresholds"]
    lat = xr.open_zarr(base / args.latitude_from)["latitude"].values.astype(np.float64)
    w = np.cos(np.deg2rad(lat))
    print(f">>> archive    {arch}")
    print(f">>> thresholds {args.thresholds}  (period {period!r})")
    print(f">>> latitude   {lat.size} rows, {lat[0]} .. {lat[-1]}\n")
    if "global_audit" not in period:
        print("!! thresholds file is not the audit reference year -- refusing")
        return 2

    band_rows = [lat < -60, (lat >= -60) & (lat < -30), (lat >= -30) & (lat < 0),
                 (lat >= 0) & (lat < 30), (lat >= 30) & (lat < 60), lat >= 60]
    area = np.array([w[m].sum() for m in band_rows]) / w.sum()
    nh, sh = lat > 0, lat < 0

    out = []
    print("   " + f"{'diagnostic':<22}{'sev':<5}{'global %':>9}  "
          + "".join(f"{b:>9}" for b in BANDS) + f"{'NH/SH':>8}")
    t0 = time.time()
    for f in sorted(arch.glob("*.npz")):
        name = f.stem
        if name not in thresholds:
            continue
        z = np.load(f)
        v, li = z["values"], z["lat_index"]
        w_total = float(z["w_total"])
        for sev, target in SEVS.items():
            i0 = int(np.searchsorted(v, thresholds[name][sev], side="left"))
            cnt = np.bincount(li[i0:].astype(np.int64), minlength=lat.size)
            mass = cnt * w
            frac = mass.sum() / w_total * 100
            share = np.array([mass[m].sum() for m in band_rows]) / mass.sum()
            dens = share / area
            ratio = mass[nh].sum() / mass[sh].sum() if mass[sh].sum() else float("inf")
            out.append(dict(diagnostic=name, severity=sev, global_pct=frac,
                            **{f"density_{b}": d for b, d in zip(BANDS, dens)},
                            nh_sh_ratio=ratio))
            print(f"   {name:<22}{sev[:4]:<5}{frac:>9.4f}  "
                  + "".join(f"{d:>9.2f}" for d in dens) + f"{ratio:>8.2f}")
        del v, li, z
    print(f"\n   {len(out) // 3} diagnostics in {(time.time() - t0) / 60:.1f} min")
    print("   density = band share of the exceedance / band share of the area (1.00 = even)")
    print("   W&S (2022) p. 1428-9: -Ri and CP tropical; TI1 two midlatitude peaks; "
          "UBF stronger in the NH; F2D no clear midlatitude peaks")
    Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.csv, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out[0]))
        wr.writeheader()
        for r in out:
            wr.writerow({k: (f"{x:.6g}" if isinstance(x, float) else x) for k, x in r.items()})
    print(f"   wrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
