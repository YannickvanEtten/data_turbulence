"""
main.py -- the one entry point. Run from the repository root inside the pixi env:

    pixi run python essentials/main.py <command> [options]

Stages, in order (each writes only under config.OUT_ROOT, except download):

  1  download-na      --index N | --year Y --month M     one North Atlantic month
     download-global  --index N (0-35)                     one global day-block of 2000
     merge-global                                          join the blocks into 12 months
  2  diagnostics      --domain north_atlantic|global --index N | --year Y --month M
                      --convention baseline|audit
  3  thresholds       --convention ...                     from the 12 global stores
  4  trends           --convention ...                     1979-2020, all seasons

Checks against the ORIGINAL pipeline (read-only, nothing is overwritten):

     compare-store       NEW.zarr OLD.zarr        bit-for-bit, all 21 diagnostics
     compare-thresholds  NEW.json OLD.json        largest relative difference
     thresholds --derived <old derived/global_audit> ...   re-calibrate old stores
     trends     --derived <old derived/north_atlantic_audit> --thresholds <old json>

Small test (no ADA needed):

     demo --input some_day.grib --convention audit [--thresholds file.json]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))

import config as C          # noqa: E402
import diagnostics as D     # noqa: E402


def month_from_index(index: int) -> tuple[int, int]:
    """SLURM array index -> (year, month); 0 = 1979-01, 503 = 2020-12."""
    return C.TREND_YEARS[0] + index // 12, index % 12 + 1


def cmd_download_na(a):
    import download
    year, month = (a.year, a.month) if a.index is None else month_from_index(a.index)
    print(f"north_atlantic {year}-{month:02d}: {download.download_month(year, month, a.force)}")


def cmd_download_global(a):
    import download
    month, block = a.index // 3 + 1, a.index % 3
    result = download.download_global_block(C.CALIBRATION_YEAR, month, block, a.force)
    print(f"global {C.CALIBRATION_YEAR}-{month:02d} block {block}: {result}")


def cmd_merge_global(a):
    import download
    for month in range(1, 13):
        print(f"global {C.CALIBRATION_YEAR}-{month:02d}: {download.merge_blocks(C.CALIBRATION_YEAR, month)}")


def cmd_diagnostics(a):
    if a.input:
        src, dst = Path(a.input), Path(a.output)
    else:
        if a.domain == "global":
            year, month = C.CALIBRATION_YEAR, (a.index + 1 if a.index is not None else a.month)
        else:
            year, month = (a.year, a.month) if a.index is None else month_from_index(a.index)
        src = C.raw_month_path(a.domain, year, month)
        dst = C.derived_dir(a.domain, a.convention) / C.month_store_name(a.domain, year, month)
    return D.run_month(src, dst, a.convention, a.chunk_days, skip_if_current=not a.no_skip)


def cmd_thresholds(a):
    import datetime as dt
    import thresholds as T
    derived = Path(a.derived) if a.derived else C.derived_dir("global", a.convention)
    stamp = dt.date.today().isoformat()
    out = Path(a.out) if a.out else C.OUT_ROOT / "calibration" / f"thresholds_{a.convention}_{stamp}.json"
    archive = (C.OUT_ROOT / "calibration" / f"tails_{a.convention}_{stamp}") if a.archive else None
    T.calibrate(derived, out, a.convention, archive)
    if a.compare_with:
        return 0 if T.compare_threshold_files(out, a.compare_with) else 1


def cmd_trends(a):
    import trends
    derived = Path(a.derived) if a.derived else C.derived_dir("north_atlantic", a.convention)
    if a.thresholds:
        tfile = Path(a.thresholds)
    else:   # newest file written by `thresholds`, else the original pipeline's file
        mine = sorted((C.OUT_ROOT / "calibration").glob(f"thresholds_{a.convention}_*.json"))
        tfile = mine[-1] if mine else C.OLD_THRESHOLDS[a.convention]
    years = range(a.years[0], a.years[1] + 1) if a.years else None
    trends.run(derived, tfile, a.convention, C.OUT_ROOT / "results", a.seasons, years)


def cmd_compare_store(a):
    new, old = xr.open_zarr(a.new), xr.open_zarr(a.old)
    all_same = True
    for k in C.DIAGNOSTIC_KEYS:
        x = new[k].transpose("latitude", "longitude", "time").values
        y = old[k].transpose("latitude", "longitude", "time").values
        same = x.shape == y.shape and np.array_equal(x, y, equal_nan=True)
        all_same &= same
        n_diff = int((~((x == y) | (np.isnan(x) & np.isnan(y)))).sum()) if x.shape == y.shape else -1
        print(f"  {k:24s} {'identical' if same else f'DIFFERS in {n_diff} cells'}")
    print("ALL 21 IDENTICAL" if all_same else "DIFFERENCES FOUND")
    return 0 if all_same else 1


def cmd_compare_thresholds(a):
    import thresholds as T
    return 0 if T.compare_threshold_files(a.new, a.old, a.tolerance) else 1


def cmd_demo(a):
    """21 diagnostics on a small file, and (optionally) exceedance frequencies."""
    import thresholds as T
    ds = D.compute_month(D.open_era5(a.input), C.CONVENTIONS[a.convention], a.chunk_days)
    w = xr.DataArray(T.latitude_weights(ds["latitude"]), dims="latitude")
    print(f"\n{'diagnostic':24s}{'median':>12s}{'p99':>12s}")
    for k in C.DIAGNOSTIC_KEYS:
        q = T.weighted_percentile(ds[k].values, w.broadcast_like(ds[k]).values, [50, 99])
        print(f"{k:24s}{q[0]:>12.4g}{q[1]:>12.4g}")
    if a.thresholds:
        thr, _, _ = T.load_thresholds(a.thresholds)
        print(f"\nexceedance frequency over the whole file (cos-latitude weighted)")
        print(f"{'diagnostic':24s}" + "".join(f"{C.SEVERITY_LABEL[s]:>9s}" for s in C.SEVERITIES))
        for k in C.DIAGNOSTIC_KEYS:
            x = ds[k]
            wb = w.broadcast_like(x).where(x.notnull())
            row = [float(((x >= thr[k][s]) * wb).sum() / wb.sum()) for s in C.SEVERITIES]
            print(f"{k:24s}" + "".join(f"{r:>9.3%}" for r in row))
    if a.output:
        ds.attrs.update(D.convention_attrs(a.convention))
        ds.to_zarr(a.output, mode="w", zarr_format=2, consolidated=True)
        print(f"wrote {a.output}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("download-na")
    p.add_argument("--index", type=int)
    p.add_argument("--year", type=int)
    p.add_argument("--month", type=int)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_download_na)

    p = sub.add_parser("download-global")
    p.add_argument("--index", type=int, required=True, help="0-35: month*3 + block")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_download_global)

    p = sub.add_parser("merge-global")
    p.set_defaults(func=cmd_merge_global)

    p = sub.add_parser("diagnostics")
    p.add_argument("--domain", choices=list(C.DOMAINS), default="north_atlantic")
    p.add_argument("--index", type=int, help="NA: 0-503 (1979-01..2020-12); global: 0-11")
    p.add_argument("--year", type=int)
    p.add_argument("--month", type=int)
    p.add_argument("--input", help="explicit GRIB path (overrides --domain/--index)")
    p.add_argument("--output", help="explicit .zarr path (with --input)")
    p.add_argument("--convention", choices=list(C.CONVENTIONS), default=C.DEFAULT_CONVENTION)
    p.add_argument("--chunk-days", type=int, default=1)
    p.add_argument("--no-skip", action="store_true", help="recompute even if the store is current")
    p.set_defaults(func=cmd_diagnostics)

    p = sub.add_parser("thresholds")
    p.add_argument("--convention", choices=list(C.CONVENTIONS), default=C.DEFAULT_CONVENTION)
    p.add_argument("--derived", help="folder with the 12 diagnostics_glob_2000-MM.zarr stores")
    p.add_argument("--out", help="output JSON")
    p.add_argument("--archive", action="store_true", help="also write the tail .npz archive (~41 GB)")
    p.add_argument("--compare-with", help="an existing thresholds JSON to compare against")
    p.set_defaults(func=cmd_thresholds)

    p = sub.add_parser("trends")
    p.add_argument("--convention", choices=list(C.CONVENTIONS), default=C.DEFAULT_CONVENTION)
    p.add_argument("--derived", help="folder with the 504 diagnostics_na_YYYY-MM.zarr stores")
    p.add_argument("--thresholds", help="thresholds JSON (default: newest from `thresholds`, "
                                        "else the original file for this convention)")
    p.add_argument("--seasons", nargs="+", default=list(C.SEASON_MONTHS), choices=list(C.SEASON_MONTHS))
    p.add_argument("--years", nargs=2, type=int, metavar=("FIRST", "LAST"))
    p.set_defaults(func=cmd_trends)

    p = sub.add_parser("compare-store")
    p.add_argument("new")
    p.add_argument("old")
    p.set_defaults(func=cmd_compare_store)

    p = sub.add_parser("compare-thresholds")
    p.add_argument("new")
    p.add_argument("old")
    p.add_argument("--tolerance", type=float, default=1e-6)
    p.set_defaults(func=cmd_compare_thresholds)

    p = sub.add_parser("demo")
    p.add_argument("--input", required=True)
    p.add_argument("--convention", choices=list(C.CONVENTIONS), default=C.DEFAULT_CONVENTION)
    p.add_argument("--thresholds")
    p.add_argument("--output")
    p.add_argument("--chunk-days", type=int, default=1)
    p.set_defaults(func=cmd_demo)

    args = ap.parse_args()
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
