"""
cat_month_frames.py -- STEP 1 of 2 (runs on ADA): reduce one month of the 21
stored diagnostics to a small "frames" file that can be animated anywhere.

Why two steps: the ADA pixi environment has numpy/xarray/zarr but no
matplotlib, cartopy or ffmpeg (checked against pixi.lock). Rather than change
the reproducibility record, ADA does the heavy reading here and writes a file
of ~10-40 MB; cat_month_render.py draws it on any machine with matplotlib.

What is stored, per 3-hourly step and 0.25 deg cell (time, lat, lon):
  count_<severity>   uint8   how many diagnostics are >= their threshold, for all
                             five severities (light .. severe)
  n_defined          uint8   how many diagnostics are defined (not NaN) there
                             -> Ebar = count / n_defined, exactly as in trends.py
  wind               uint8   200 hPa wind speed in m/s (rounded) for jet contours
  bits_<severity>    uint8   per-diagnostic exceedance at --per-diagnostic
                             severity, bit-packed along longitude (21, T, Y, ceil(X/8))

Usage on ADA (repo root, pixi env), one month or a whole year:
    pixi run python ada/cat_month_frames.py --year 2020 --month 2
    pixi run python ada/cat_month_frames.py --year 2020            # all 12 months
Locally on the one-day file in the repo (preview):
    python ada/cat_month_frames.py --nc cat_outputs/diagnostics.nc \
        --thresholds calibration/thresholds_2026-09-11.json --out preview
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import xarray as xr

BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
DEFAULT_DERIVED = BASE / "derived" / "north_atlantic_audit"      # THE DATASET (FINAL_DATASET.md)
DEFAULT_THRESHOLDS = BASE / "calibration" / "thresholds_2026-09-11.json"
DEFAULT_OUT = BASE / "explore" / "frames"

SEVERITIES = ["light", "light_to_moderate", "moderate", "moderate_to_severe", "severe"]
DIAGNOSTICS = [
    "magnitude_pv", "colson_panofsky", "brown1", "temperature_gradient",
    "horizontal_divergence", "vertical_wind_shear", "endlich", "deformation",
    "wind_speed", "ngm2", "negative_richardson", "rva_magnitude", "ubf", "nva",
    "brown2", "vorticity_squared", "ti1", "ngm1", "ti2", "f2d", "ncsu1",
]


def open_month(derived: Path, year: int, month: int, nc: str | None) -> xr.Dataset:
    if nc:
        ds = xr.open_dataset(nc)
    else:
        store = derived / f"diagnostics_na_{year}-{month:02d}.zarr"
        if not store.exists():
            raise FileNotFoundError(store)
        ds = xr.open_zarr(store)
    return ds.transpose("time", "latitude", "longitude")


def reduce_month(ds: xr.Dataset, thr: dict, per_diag_sev: str | None) -> dict:
    T, Y, X = ds.sizes["time"], ds.sizes["latitude"], ds.sizes["longitude"]
    counts = {s: np.zeros((T, Y, X), np.uint8) for s in SEVERITIES}
    n_def = np.zeros((T, Y, X), np.uint8)
    bits = []
    present = [k for k in DIAGNOSTICS if k in ds]
    for key in present:                      # one diagnostic in memory at a time
        x = ds[key].values.astype(np.float32)
        ok = np.isfinite(x)
        n_def += ok
        for s in SEVERITIES:
            counts[s] += ok & (x >= thr[key][s])
        if per_diag_sev:
            bits.append(np.packbits(ok & (x >= thr[key][per_diag_sev]), axis=-1))
        if key == "wind_speed":
            wind = np.clip(np.rint(np.nan_to_num(x)), 0, 255).astype(np.uint8)
    out = {f"count_{s}": counts[s] for s in SEVERITIES}
    out.update(n_defined=n_def, diagnostics=np.array(present),
               latitude=ds.latitude.values.astype(np.float32),
               longitude=ds.longitude.values.astype(np.float32),
               time=ds.time.values.astype("datetime64[ns]").astype(np.int64))
    if "wind_speed" in present:
        out["wind"] = wind
    if per_diag_sev:
        out[f"bits_{per_diag_sev}"] = np.stack(bits)
        out["bits_n_longitude"] = np.array(X)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, default=2020)
    ap.add_argument("--month", type=int, default=None, help="omit for all 12 months of --year")
    ap.add_argument("--derived", default=str(DEFAULT_DERIVED))
    ap.add_argument("--thresholds", default=str(DEFAULT_THRESHOLDS))
    ap.add_argument("--per-diagnostic", default="moderate", choices=SEVERITIES + ["none"],
                    help="severity at which the 21 single-diagnostic fields are kept")
    ap.add_argument("--nc", default=None, help="read one netCDF file instead of a zarr month")
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    thr_doc = json.loads(Path(args.thresholds).read_text())
    thr = thr_doc["thresholds"]
    per_diag = None if args.per_diagnostic == "none" else args.per_diagnostic
    months = [args.month] if (args.month or args.nc) else list(range(1, 13))
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    for month in months:
        ds = open_month(Path(args.derived), args.year, month or 1, args.nc)
        red = reduce_month(ds, thr, per_diag)
        tag = np.datetime_as_string(ds.time.values[0], unit="M")
        path = out_dir / f"cat_frames_{tag}.npz"
        np.savez_compressed(
            path, **red,
            thresholds_file=np.array(Path(args.thresholds).name),
            thresholds_created=np.array(thr_doc.get("created", "")),
            source=np.array(args.nc or str(Path(args.derived) / f"diagnostics_na_{tag}.zarr")),
        )
        mean_light = (red["count_light"] / np.maximum(red["n_defined"], 1)).mean()
        print(f"{tag}: {red['n_defined'].shape[0]} steps, {len(red['diagnostics'])} diagnostics, "
              f"mean Ebar(light) {mean_light:.4f} -> {path} ({path.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
