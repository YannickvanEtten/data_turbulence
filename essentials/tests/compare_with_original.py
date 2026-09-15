"""
Compare essentials/diagnostics.py with the ORIGINAL pipeline on one ERA5 file.

Needs the original repo (the parent folder of essentials/) and rojak @ 25b8685,
i.e. run it inside the project's pixi environment:

    pixi run python essentials/tests/compare_with_original.py <file.grib> [--chunk-days 1]

For each convention it runs the original ada/diagnostics_global.compute_chunked
with the matching FixSet / f2d variant, runs essentials' compute_month, and
prints per diagnostic: exact bitwise match of the float32 output, the largest
relative difference, and whether the NaN pattern is identical.
Exit code 0 only if every diagnostic matches bit for bit under both conventions.
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
ESSENTIALS = HERE.parents[1]
ORIGINAL = ESSENTIALS.parent
sys.path.insert(0, str(ESSENTIALS))

import config as C                     # noqa: E402
import diagnostics as D                # noqa: E402


def load_original():
    spec = importlib.util.spec_from_file_location("dg", ORIGINAL / "ada" / "diagnostics_global.py")
    dg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dg)
    diag = dg._load("diagnostics_original", "2_diagnostics.py")
    return dg, diag


def original_settings(diag, name):
    conv = C.CONVENTIONS[name]
    fixes = diag.FixSet(meridional_metric=(conv["dy"] == "map"),
                        four_variable_provenance=(conv["fields"] == "computed"),
                        endlich_component_shear=(conv["endlich"] == "components"))
    return fixes, conv["f2d_variant"]


def compare(old, new):
    rows, all_equal = [], True
    for k in C.DIAGNOSTIC_KEYS:
        a = old[k].transpose("latitude", "longitude", "time").values
        b = new[k].transpose("latitude", "longitude", "time").values
        same_nan = np.array_equal(np.isnan(a), np.isnan(b))
        equal = np.array_equal(a, b, equal_nan=True)
        fin = np.isfinite(a) & np.isfinite(b)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.abs(a[fin].astype(float) - b[fin]) / np.maximum(np.abs(a[fin].astype(float)), 1e-300)
        rows.append((k, equal, same_nan, float(rel.max()) if rel.size else 0.0,
                     int((a[fin] != b[fin]).sum()), int(fin.sum())))
        all_equal &= equal
    return rows, all_equal


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("grib")
    ap.add_argument("--chunk-days", type=int, default=1)
    ap.add_argument("--conventions", nargs="+", default=list(C.CONVENTIONS) + ["baseline_endlich_components"])
    args = ap.parse_args()

    C.CONVENTIONS.setdefault("baseline_endlich_components",
                             {**C.CONVENTIONS["baseline"], "endlich": "components"})
    dg, diag = load_original()
    raw = D.open_era5(args.grib)
    ok = True
    for name in args.conventions:
        fixes, variant = original_settings(diag, name)
        old, failures = dg.compute_chunked(diag, raw, args.chunk_days, verbose=False,
                                           f2d_variant=variant, fixes=fixes)
        assert not failures, failures
        new = D.compute_month(raw, C.CONVENTIONS[name], args.chunk_days, verbose=False)
        rows, equal = compare(old, new)
        ok &= equal
        print(f"\n=== {name}: {C.CONVENTIONS[name]}   original FixSet={fixes.label()} f2d={variant}")
        print(f"{'diagnostic':24s}{'bitwise':>9s}{'NaNs same':>11s}{'max rel diff':>14s}{'cells differ':>14s}")
        for k, eq, nan_eq, rel, ndiff, nfin in rows:
            print(f"{k:24s}{str(eq):>9s}{str(nan_eq):>11s}{rel:>14.2e}{ndiff:>8d}/{nfin}")
        print(f"--> {'ALL 21 BIT-IDENTICAL' if equal else 'DIFFERENCES FOUND'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
