#!/usr/bin/env python
"""
ada/diagnostics_global.py
=========================
Compute the 21 diagnostics on ONE global calibration file and save them in the
form the calibration actually needs -- nothing more.

    pixi run python ada/diagnostics_global.py <in.grib> <out.zarr> [--chunk-days N]

WHY NOT JUST USE 3_pipeline.py?
-------------------------------
`3_pipeline.run()` is built for the analysis path and does three things that
are wrong for a calibration input:

  1. `save_outputs()` writes the SAME data as NetCDF *and* Zarr. For the
     analysis run that is a known waste (STATUS 10.12); for twelve global
     months it is ~134 GB instead of ~33 GB.
  2. It saves float64. The diagnostics are computed from 16-bit-packed GRIB
     and feed a percentile -- float32 is well past sufficient and halves both
     the file and the peak memory.
  3. It builds the W&J comparison table, which is a per-month statement about
     the *analysis* box. On a global month it is not meaningful.

Layer 2 needs exactly one thing from each month: every diagnostic value, with
its latitude. That is what this writes.

MEMORY, AND WHY --chunk-days EXISTS
-----------------------------------
This is the first time the diagnostics see a global grid, and the shape is
very different from the North Atlantic month that produced the only MaxRSS
measurement on record (14.8 GB on 121 x 301 x 248).

    NA month      121 x  301 x 248  =  9.0e6 points
    global file   721 x 1440 x  32  =  3.3e7 points   (3.7x)

A linear extrapolation says ~55 GB. That is an extrapolation across a
different array SHAPE, so do not trust it -- but more importantly, a job
asking for that much waits for a large contiguous free block. On a busy defq
(47 jobs queued, every node in 'mix') a 120 GB request was scheduled 22 hours
out, while the 4-8 GB download jobs started instantly. Memory footprint is the
single biggest lever on time-to-start here.

--chunk-days processes the file in time slices, so the peak is set by the
chunk rather than by the file. --chunk-days 1 should land near the NA month's
~15 GB, which fits ANY defq node including the 59 GB ones.

THE OVERLAP BUFFER -- WHY CHUNKING COSTS NOTHING
------------------------------------------------
Chunking naively would corrupt F2D. `frontogenesis_isentropic` takes a CENTRED
material derivative, so it drops the first and last timestep of whatever it is
given: chop 32 steps into four blocks of 8 and F2D loses 8 steps instead of 2.

So each chunk is computed with ONE EXTRA TIMESTEP on each side, and the result
is trimmed back to the target range afterwards. Every diagnostic then sees the
same neighbours it would have seen in an unchunked run, and the concatenated
output is identical to processing all 32 steps at once -- F2D included, losing
only the two steps at the file's genuine ends.

This is the same device `chunk_stitch.py` applies at month boundaries. Here it
is applied within a file; the reason is identical.

Peak RSS is printed per chunk and at the end, so the number lands in the job
log without depending on SLURM step accounting -- which on this cluster has
not been returning MaxRSS.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[1]

STEPS_PER_DAY = 8          # 3-hourly
OVERLAP = 1                # timesteps of buffer each side; F2D needs exactly 1


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
    """Peak RSS of this process so far, in GB. Linux reports ru_maxrss in kB."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def existing_output_matches(out_path: Path, f2d_variant: str) -> tuple[bool, str]:
    """Is there already a COMPLETE output here built with THIS configuration?

    The production array is submitted repeatedly -- while the download is still
    running, and again after it finishes -- so it has to distinguish three
    states, not two:

        no output                  -> compute
        output from an older config -> RECOMPUTE, do not skip
        output from this config     -> skip in a second

    The middle case is the one a plain `if [ -d "$OUT" ]` gets wrong, and it is
    not hypothetical: every zarr written before 2026-08-29 holds DEF^2 rather
    than DEF, and every zarr written before 2026-08-30 holds f2d variant A
    rather than C. A directory-exists check would silently preserve both and
    hand back a dataset that is half one convention and half the other -- the
    same shape of failure as Q-INTEG-3's silent level substitution.

    Completeness is judged by the consolidated-metadata marker, so a store
    killed part-way through a write is recomputed rather than trusted.
    """
    if not out_path.exists():
        return False, "no output yet"

    complete = (out_path / ".zmetadata").exists() or (out_path / "zarr.json").exists()
    if not complete:
        return False, "incomplete store (no consolidated metadata) — will rewrite"

    zattrs = out_path / ".zattrs"
    if not zattrs.exists():
        return False, "no .zattrs — written before 2026-08-29, config unknown"
    try:
        attrs = json.loads(zattrs.read_text())
    except (OSError, ValueError) as exc:
        return False, f"unreadable .zattrs ({type(exc).__name__}) — will rewrite"

    have_variant = attrs.get("f2d_variant")
    have_def = attrs.get("deformation_convention")
    if have_variant is None or have_def is None:
        return False, "provenance attributes absent — pre-2026-08-29 output"
    if have_variant != f2d_variant:
        return False, f"f2d_variant {have_variant!r} != requested {f2d_variant!r}"
    return True, f"f2d_variant={have_variant}, deformation={have_def}"


def _time_dim(da: xr.DataArray) -> str | None:
    for candidate in ("time", "valid_time", "step"):
        if candidate in da.dims:
            return candidate
    return None


def compute_chunked(diag, ds_raw: xr.Dataset, chunk_days: int, verbose=True,
                    f2d_variant: str | None = None):
    """Compute all 21 diagnostics, optionally in overlapping time chunks.

    chunk_days = 0 means one pass over the whole file.

    Returns (diagnostics, failures). The output is time-identical to an
    unchunked run: each chunk is computed with OVERLAP extra steps on each
    side and then reindexed onto exactly its target timestamps, so a
    diagnostic that drops edge steps (F2D) gets its neighbours from the
    adjacent chunk rather than losing them.
    """
    tname = None
    for candidate in ("time", "valid_time"):
        if candidate in ds_raw.dims:
            tname = candidate
            break
    if tname is None:
        raise RuntimeError(f"no recognised time dimension in {list(ds_raw.dims)}")

    n_steps = ds_raw.sizes[tname]
    all_times = ds_raw[tname].values

    if chunk_days <= 0:
        bounds = [(0, n_steps)]
    else:
        step = chunk_days * STEPS_PER_DAY
        bounds = [(a, min(a + step, n_steps)) for a in range(0, n_steps, step)]

    if verbose:
        print(f"    {n_steps} timesteps -> {len(bounds)} chunk(s) "
              f"of up to {chunk_days * STEPS_PER_DAY if chunk_days > 0 else n_steps} "
              f"steps, overlap {OVERLAP}")

    per_chunk: dict[str, list[xr.DataArray]] = {}
    all_failures = []

    for i, (a, b) in enumerate(bounds):
        # Widen by the overlap buffer, clipped at the file's real ends.
        lo, hi = max(0, a - OVERLAP), min(n_steps, b + OVERLAP)
        target_times = all_times[a:b]

        sub = ds_raw.isel({tname: slice(lo, hi)})
        catdata = diag.prepare_for_rojak(sub)
        diagnostics, failures = diag.compute_all_21(
            catdata, target_level=200,
            f2d_variant=f2d_variant or diag.F2D_DEFAULT_VARIANT)
        all_failures.extend(failures)

        for key, da in diagnostics.items():
            dt = _time_dim(da)
            if dt is not None:
                # Reindex onto EXACTLY the target timestamps. Diagnostics that
                # kept every step get trimmed; F2D, which dropped the buffer
                # steps, lines up on the interior and takes NaN only where the
                # file genuinely has no neighbour (first and last step overall).
                da = da.reindex({dt: target_times})
            elif i > 0:
                # No time dimension at all -- one copy is enough.
                continue
            per_chunk.setdefault(key, []).append(da.astype(np.float32))

        if verbose:
            print(f"    chunk {i + 1}/{len(bounds)}  steps {a}:{b} "
                  f"(computed {lo}:{hi})   rss {peak_rss_gb():.1f} GB")

        del sub, catdata, diagnostics

    merged: dict[str, xr.DataArray] = {}
    for key, parts in per_chunk.items():
        dt = _time_dim(parts[0])
        merged[key] = xr.concat(parts, dim=dt) if (dt and len(parts) > 1) else parts[0]

    return merged, all_failures


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="ERA5 GRIB file")
    p.add_argument("output", help="output .zarr path")
    p.add_argument("--chunk-days", type=int, default=0,
                   help="process this many days at a time, with a 1-timestep "
                        "overlap buffer so the result is identical to an "
                        "unchunked run. 0 (default) = whole file at once. "
                        "Use 1 to keep peak memory near the North Atlantic "
                        "month's ~15 GB and schedule on any node.")
    p.add_argument("--f2d-variant", default=None,
                   help="which reading of Sharman A9 to use for #20 (see "
                        "2_diagnostics.F2D_VARIANTS and FORMULA_AUDIT.md 4). "
                        "Default is 2_diagnostics.F2D_DEFAULT_VARIANT. Only "
                        "change this once ada/check_f2d_variants.py has "
                        "reported -- and note the choice is recorded in the "
                        "output zarr's attributes either way.")
    p.add_argument("--skip-if-matching", action="store_true",
                   help="exit 0 immediately if a COMPLETE output already exists "
                        "that was built with this same configuration. An output "
                        "from an older configuration is recomputed, not skipped "
                        "— see existing_output_matches(). This is what makes the "
                        "production array safe to resubmit while the download is "
                        "still running.")
    args = p.parse_args()

    in_path, out_path = Path(args.input), Path(args.output)
    if not in_path.exists():
        print(f"!! input not found: {in_path}")
        return 1

    diag = _load("diagnostics", "2_diagnostics.py")

    if args.skip_if_matching:
        variant = args.f2d_variant or diag.F2D_DEFAULT_VARIANT
        matches, why = existing_output_matches(out_path, variant)
        print(f">>> existing output: {why}")
        if matches:
            print(f"    SKIP — {out_path.name} is already current")
            return 0

    t0 = time.time()
    print(f">>> loading {in_path.name}  ({in_path.stat().st_size / 1e9:.2f} GB)")
    ds_raw = diag.load_era5(in_path)
    print(f"    raw vars : {sorted(ds_raw.data_vars)}")
    print(f"    dims     : {dict(ds_raw.sizes)}")
    print(f"    [rss after load: {peak_rss_gb():.1f} GB]")

    f2d_variant = args.f2d_variant or diag.F2D_DEFAULT_VARIANT
    if f2d_variant not in diag.F2D_VARIANTS:
        print(f"!! unknown --f2d-variant {f2d_variant!r}; "
              f"expected one of {sorted(diag.F2D_VARIANTS)}")
        return 1

    print(f"\n>>> computing 21 diagnostics at 200 hPa "
          f"(chunk_days={args.chunk_days or 'whole file'})")
    print(f"    f2d variant {f2d_variant}: {diag.F2D_VARIANTS[f2d_variant]}")
    diagnostics, failures = compute_chunked(diag, ds_raw, args.chunk_days,
                                            f2d_variant=f2d_variant)
    print(f"    [rss after compute: {peak_rss_gb():.1f} GB]")

    if failures:
        # A failed diagnostic becomes an all-NaN placeholder, and
        # calibration.compute_thresholds REFUSES to build a threshold from one.
        # Better to see it here than to have Layer 2 die twelve months later.
        print(f"\n!! {len(failures)} diagnostic failure(s) on this file:")
        for f in failures:
            print(f"     #{f.wj_number:>2} {f.key:<23} {f.exception_type}: {f.message}")
        print("   These become all-NaN placeholders and CANNOT be calibrated.")

    print("\n>>> writing zarr (float32)")
    ds_out = xr.Dataset(diagnostics)

    # Dataset-level provenance. Both entries record a decision that is
    # invisible in the numbers themselves and unrecoverable afterwards: which
    # reading of A9 produced #20, and whether #8 holds DEF or DEF^2. A zarr
    # written before 2026-08-29 carries neither attribute, which is itself the
    # signal that it is the old convention.
    ds_out.attrs.update({
        "f2d_variant": f2d_variant,
        "f2d_variant_formula": diag.F2D_VARIANTS[f2d_variant],
        "deformation_convention": "DEF (un-squared, Sharman A17)",
        "source_file": in_path.name,
        "target_level_hPa": 200,
    })

    empty = [k for k, da in ds_out.data_vars.items()
             if not np.isfinite(da.values).any()]
    for key in empty:
        print(f"    !! {key}: NO finite values -- will break Layer 2")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds_out.to_zarr(out_path, mode="w", zarr_format=2)

    size_gb = sum(f.stat().st_size for f in out_path.rglob("*") if f.is_file()) / 1e9
    elapsed = time.time() - t0

    print("\n=== RESULT ===")
    print(f"  output    : {out_path}  ({size_gb:.2f} GB)")
    print(f"  variables : {len(ds_out.data_vars)} of 21")
    print(f"  wall time : {int(elapsed // 60)} min {int(elapsed % 60)} s")
    print(f"  PEAK RSS  : {peak_rss_gb():.1f} GB   <-- size the remaining months from this")
    return 0 if not (failures or empty) else 2


if __name__ == "__main__":
    raise SystemExit(main())
