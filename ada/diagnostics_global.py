#!/usr/bin/env python
"""
ada/diagnostics_global.py
=========================
Compute the 21 diagnostics on ONE global calibration file and save them in the
form the calibration actually needs — nothing more.

    pixi run python ada/diagnostics_global.py <in.grib> <out.zarr>

WHY NOT JUST USE 3_pipeline.py?
-------------------------------
`3_pipeline.run()` is built for the analysis path and does three things that
are wrong for a calibration input:

  1. `save_outputs()` writes the SAME data as NetCDF *and* Zarr. For the
     analysis run that is a known waste (STATUS 10.12); for twelve global
     months it is ~134 GB instead of ~33 GB.
  2. It saves float64. The diagnostics are computed from 16-bit-packed GRIB
     and feed a percentile — float32 is well past sufficient and halves both
     the file and the peak memory.
  3. It builds the W&J comparison table, which is a per-month statement about
     the *analysis* box. On a global month it is not meaningful.

Layer 2 needs exactly one thing from each month: every diagnostic value, with
its latitude. That is what this writes.

MEMORY
------
This is the first time the diagnostics see a global grid, and the shape is
very different from the North Atlantic month that produced the only MaxRSS
measurement on record (14.8 GB on 121 x 301 x 248).

    NA month      121 x  301 x 248  =  9.0e6 points
    global file   721 x 1440 x  32  =  3.3e7 points   (3.7x)

A linear extrapolation says ~55 GB, but that is an extrapolation across a
different array *shape*, and rojak's gradient operators need not scale with
point count alone. Do not trust it. Run ONE month with a generous request,
read the real MaxRSS off sacct, and size the rest from that — the same way
03_diagnostics.sbatch was sized after the fact rather than guessed.

The script prints its own peak RSS at the end so the number is in the log
even if you lose the job id.
"""
from __future__ import annotations

import importlib.util
import resource
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[1]


def _load(name: str, filename: str):
    """Load a module whose filename is not a legal Python identifier.

    Registering in sys.modules BEFORE exec_module is mandatory: 2_diagnostics.py
    uses `from __future__ import annotations`, so @dataclass resolves its field
    annotations through sys.modules[cls.__module__]. Leaving that None is the
    Phase-0 crash in STATUS.md.
    """
    spec = importlib.util.spec_from_file_location(name, REPO / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def peak_rss_gb() -> float:
    """Peak RSS of this process, in GB. Linux reports ru_maxrss in kB."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        print("usage: diagnostics_global.py <in.grib> <out.zarr>")
        return 1

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    if not in_path.exists():
        print(f"!! input not found: {in_path}")
        return 1

    diag = _load("diagnostics", "2_diagnostics.py")

    t0 = time.time()
    print(f">>> loading {in_path.name}  ({in_path.stat().st_size / 1e9:.2f} GB)")
    ds_raw = diag.load_era5(in_path)
    print(f"    raw vars : {sorted(ds_raw.data_vars)}")

    catdata = diag.prepare_for_rojak(ds_raw)
    sizes = dict(catdata._dataset.sizes)
    print(f"    dims     : {sizes}")
    n_points = int(np.prod([v for k, v in sizes.items()
                            if k in ("latitude", "longitude", "time")]))
    print(f"    points   : {n_points:,} per diagnostic per level")
    print(f"    [rss after load: {peak_rss_gb():.1f} GB]")

    print("\n>>> computing 21 diagnostics at 200 hPa")
    diagnostics, failures = diag.compute_all_21(catdata, target_level=200)
    print(f"    [rss after compute: {peak_rss_gb():.1f} GB]")

    if failures:
        # A diagnostic that failed here would silently become an all-NaN
        # placeholder, and calibration.compute_thresholds REFUSES to build a
        # threshold from an all-NaN field. Better to see it now than to have
        # Layer 2 die twelve months later.
        print(f"\n!! {len(failures)} diagnostic(s) FAILED on this file:")
        for f in failures:
            print(f"     #{f.wj_number:>2} {f.key:<23} {f.exception_type}: {f.message}")
        print("   These become all-NaN placeholders and CANNOT be calibrated.")

    # float32: the inputs are 16-bit-packed GRIB and the output feeds a
    # percentile, so float64 buys nothing and costs 2x memory and 2x disk.
    print("\n>>> casting to float32 and writing zarr")
    ds_out = xr.Dataset({k: v.astype(np.float32) for k, v in diagnostics.items()})

    for key, da in ds_out.data_vars.items():
        finite = int(np.isfinite(da.values).sum())
        if finite == 0:
            print(f"    !! {key}: NO finite values -- will break Layer 2")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds_out.to_zarr(out_path, mode="w", zarr_format=2)

    size_gb = sum(f.stat().st_size for f in out_path.rglob("*") if f.is_file()) / 1e9
    elapsed = time.time() - t0

    print(f"\n=== RESULT ===")
    print(f"  output    : {out_path}  ({size_gb:.2f} GB)")
    print(f"  variables : {len(ds_out.data_vars)} of 21")
    print(f"  wall time : {int(elapsed // 60)} min {int(elapsed % 60)} s")
    print(f"  PEAK RSS  : {peak_rss_gb():.1f} GB   <-- size the remaining months from this")
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
