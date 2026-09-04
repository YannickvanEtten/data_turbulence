"""
ada/tail_thresholds.py
======================
Percentile calibration that does not need the whole calibration year in memory.

WHY THIS EXISTS
---------------
`calibration.compute_thresholds` ravels the entire calibration field. On the
FULL global year 2000 that is 721 x 1440 x 2928 = 3.04e9 points per diagnostic:
24 GB of float64 values, 24 GB of broadcast weights, a 24 GB int64 argsort
index, and np.unique's own sort buffer on top -- over 100 GB peak, PER
DIAGNOSTIC (CALIBRATION_REFERENCE.md 8.4). That single fact is why the
year-2000 calibration has been a 48-day sub-sample since 2026-08-29, and it is
the only thing that has to change to run the full reference year.

WHAT THIS DOES INSTEAD -- AND WHY IT IS EXACT, NOT AN ESTIMATE
--------------------------------------------------------------
The Hazen plotting position `calib_weighted_percentile` uses is

    P_i = (C_i - 0.5 * w_i) / W_total ,    C_i = cumulative weight through i

Everything BELOW a retained upper tail enters that expression through exactly
two scalars: `W_total` and `W_below`. So one streaming pass that keeps

    * every (value, latitude row) with value >= cut,
    * W_total = sum of cos(phi) over FINITE cells,
    * W_below = sum of cos(phi) over finite cells with value < cut,

reproduces `weighted_percentile` for every percentile whose interpolation
bracket lies inside the retained tail. This is not a t-digest, not a histogram
estimator and not a sub-sample: it is the same arithmetic on the same numbers,
with the discarded bulk summarised LOSSLESSLY by the only two quantities it is
able to contribute.

The one caveat, stated honestly: agreement is to floating point, not bitwise.
Measured on a synthetic latitude-dependent heavy-tailed field with NaNs and a
mass of exact ties at zero (the `nva`/`ncsu1` clipping pattern), against the
project's own `weighted_percentile` on the whole array:

    W_total   identical to the last bit
    LOG       9.3e-12   LMOG 6.2e-12   MOG 9.7e-11
    MSOG      1.2e-09   SOG  2.0e-09

The residual grows toward the extreme tail because that is where the plotting
positions thin out and dv/dP is steepest, so a last-bit difference in a
cumulative weight is amplified by np.interp. Nine to ten significant figures
on a quantity whose real uncertainty is percent-level. `--tolerance` defaults
to 1e-6 for that reason -- tight enough to catch any actual bug by many orders,
loose enough not to fail on arithmetic that is already exact.

Three details that keep it exact rather than nearly-exact:

  1. Retain `v >= cut`, never `v > cut`. np.unique consolidates ties into ONE
     plotting position; a value sitting exactly at the cut must not have its
     weight split across the boundary.
  2. Cast retained float32 to float64 before np.unique, because
     weighted_percentile does. float32 -> float64 is lossless, so the tie
     consolidation is identical.
  3. Weights are recomputed from the stored latitude ROW INDEX, using the same
     np.cos(np.deg2rad(lat)) in float64. Nothing about the weighting is
     approximated -- int16 indexes 721 rows exactly.

W_total is accumulated PER DIAGNOSTIC, not once. The 21 do not share a NaN
mask: f2d loses the two timesteps at each file's true ends that the other 20
keep (STATUS 4d, 11.11).

THE GUARD IS THE POINT OF THE FILE
----------------------------------
None of the above is true if the cut lands above the percentile you asked for.
np.interp would then silently clamp to `left=v_unique[0]` and return the cut
itself -- which looks like a perfectly plausible threshold and is not one.
So the retained weighted fraction is asserted against the most extreme
percentile requested, with margin, and the run FAILS LOUDLY rather than
returning a clamped number. This is the same principle as
`aggregate.assert_thresholds_not_hardcoded_table2`: the failure mode that
produces a believable wrong answer gets an assertion, not a comment.

THE CUT
-------
Set from a cheap stride sub-sample of timesteps, using the project's own
`weighted_percentile` on exactly the workload that already runs today. The
default targets 10 % retention with the cut placed at p88 for margin, which
covers:

  * Prosser's ladder, p97.0 .. p99.9 (SEVERITIES below),
  * Lee et al. (2023) Table 1, which uses p95 for MOG -- the only ERA5
    threshold table there is, and not comparable without our own p95,
  * a RANGE of candidate thresholds for GPD fitting. That is not a luxury:
    mean-residual-life and parameter-stability plots over varying thresholds
    are how an EVT threshold gets chosen, and a single frozen 5 % cut
    forecloses it.

OUTPUT
------
  calibration/thresholds_<date>.json     identical schema to calibration.py,
                                         written by calibration.save_thresholds
  calibration/tails_<date>/<diag>.npz    values (float32, ascending),
                                         lat_index (int16), and the scalars
                                         w_below / w_total / cut / n_total

The .npz set IS the peaks-over-threshold archive phase 5 needs (STATUS 6,
phase 5 table) -- it is a by-product of calibrating, not extra work.

VALIDATE IT BEFORE TRUSTING IT
------------------------------
    pixi run python ada/tail_thresholds.py \
        --derived-subdir derived/global_sub48 \
        --validate-against calibration/thresholds_2026-09-03.json

That reproduces the EXISTING 48-day thresholds through the tail path, on data
already on disk. It costs minutes and needs no CDS queue. Run it before
spending a day downloading the full year: if it disagrees, the bug is here,
not in the new data.
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
sys.path.insert(0, str(REPO))

import calibration                                             # noqa: E402
from calib_weighted_percentile import weighted_percentile      # noqa: E402


def _load(name: str, filename: str):
    """Import a numbered module (2_diagnostics.py) by path.

    Registering in sys.modules is NOT optional: 2_diagnostics.py uses
    `from __future__ import annotations`, so @dataclass resolves its field
    annotations through sys.modules[cls.__module__], and leaving that None
    kills the import with a bare 'NoneType' object has no attribute
    '__dict__'. See STATUS.md 6, phase 0.
    """
    spec = importlib.util.spec_from_file_location(name, REPO / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Williams (2017) Table 1, adopted verbatim by Prosser (2023) s2. LOCKED.
SEVERITIES = {
    "light":              97.0,
    "light_to_moderate":  99.1,
    "moderate":           99.6,
    "moderate_to_severe": 99.8,
    "severe":             99.9,
}


# ---------------------------------------------------------------------------
# Streaming machinery
# ---------------------------------------------------------------------------
def peak_rss_gb() -> float:
    """Peak RSS of this process, in GB.

    Job-step accounting is DISABLED cluster-wide on ADA -- sacct and sstat both
    return empty MaxRSS -- so every driver has to measure and print its own
    (STATUS 11.9). Without this the only way to size the full-year run is to
    submit it and watch it die.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 ** 2)


def month_paths(base: Path, subdir: str) -> list[Path]:
    paths = sorted((base / subdir).glob("diagnostics_glob_*.zarr"))
    if not paths:
        raise SystemExit(f"no diagnostics_glob_*.zarr under {base / subdir}")
    return paths


def lat_weights_from(lat: np.ndarray) -> np.ndarray:
    """cos(phi) in float64 -- byte-for-byte what calibration_check builds."""
    return np.cos(np.deg2rad(np.asarray(lat, dtype=np.float64)))


def _blocks(da: xr.DataArray, time_chunk: int):
    """Yield (nlat, nlon, nt) float32 blocks, bounded in memory by time_chunk."""
    da = da.transpose("latitude", "longitude", "time")
    nt = da.sizes["time"]
    for a in range(0, nt, time_chunk):
        yield np.asarray(da.isel(time=slice(a, min(a + time_chunk, nt))).values,
                         dtype=np.float32)


def cut_for(name: str, paths: list[Path], w_lat: np.ndarray,
            stride: int, cut_pct: float) -> float:
    """Locate the retention cut on a stride sub-sample of timesteps.

    Uses the project's own weighted_percentile, on roughly the workload that
    already runs today, so this pass cannot be the thing that does not fit.
    The value only has to be BELOW the lowest percentile we will ask for --
    the guard in compute() checks that it is, so a sloppy cut fails loudly
    instead of biasing a threshold.
    """
    vals, wts = [], []
    for p in paths:
        da = xr.open_zarr(p)[name].isel(time=slice(None, None, stride))
        block = np.asarray(da.transpose("latitude", "longitude", "time").values,
                           dtype=np.float32)
        w = np.broadcast_to(w_lat[:, None, None], block.shape)
        finite = np.isfinite(block)
        vals.append(block[finite].astype(np.float64))
        wts.append(w[finite])
    v = np.concatenate(vals)
    w = np.concatenate(wts)
    if v.size == 0:
        raise ValueError(f"diagnostic {name!r} is entirely non-finite on the "
                         f"stride sub-sample -- refusing to guess a cut.")
    return float(np.atleast_1d(weighted_percentile(v, w, np.asarray([cut_pct])))[0])


def scan_tail(name: str, paths: list[Path], w_lat: np.ndarray,
              cut: float, time_chunk: int):
    """One pass over every timestep. Returns (values, lat_index, w_below,
    w_total, n_finite).

    Weight sums are accumulated as (per-latitude count) . (per-latitude
    weight) rather than by masking a broadcast weight array, which would
    allocate a float64 copy of every block.
    """
    keep_v: list[np.ndarray] = []
    keep_i: list[np.ndarray] = []
    w_total = 0.0
    w_below = 0.0
    n_finite = 0
    nlat = w_lat.size
    row = np.arange(nlat, dtype=np.int16)

    for p in paths:
        for block in _blocks(xr.open_zarr(p)[name], time_chunk):
            finite = np.isfinite(block)
            keep = finite & (block >= cut)          # >= : never split a tie

            cnt_fin = finite.sum(axis=(1, 2))
            cnt_keep = keep.sum(axis=(1, 2))
            w_total += float(np.dot(w_lat, cnt_fin))
            w_below += float(np.dot(w_lat, cnt_fin - cnt_keep))
            n_finite += int(cnt_fin.sum())

            if cnt_keep.sum():
                keep_v.append(block[keep])          # float32, one allocation
                keep_i.append(np.broadcast_to(
                    row[:, None, None], block.shape)[keep])

    if not keep_v:
        raise ValueError(f"diagnostic {name!r} retained nothing above cut "
                         f"{cut!r} -- the cut is wrong or the field is empty.")
    return (np.concatenate(keep_v), np.concatenate(keep_i),
            w_below, w_total, n_finite)


def thresholds_from_tail(values: np.ndarray, lat_index: np.ndarray,
                         w_lat: np.ndarray, w_below: float, w_total: float,
                         percentiles: np.ndarray, name: str, margin: float):
    """Reproduce weighted_percentile from the retained tail alone.

    Returns (thresholds, sorted_values_float32, sorted_lat_index, p_first)
    where p_first is the plotting position of the smallest retained value --
    the number the guard is about.
    """
    v = values.astype(np.float64)                   # lossless from float32
    w = w_lat[lat_index]

    order = np.argsort(v, kind="stable")
    v = v[order]
    w = w[order]

    v_u, inverse = np.unique(v, return_inverse=True)
    w_u = np.zeros(v_u.shape, dtype=np.float64)
    np.add.at(w_u, inverse, w)

    cumw = w_below + np.cumsum(w_u)                 # the ONLY use of the bulk
    plotting = (cumw - 0.5 * w_u) / w_total

    p = np.asarray(percentiles, dtype=np.float64) / 100.0
    p_first = float(plotting[0])

    # THE GUARD. Without it np.interp clamps to v_u[0] and hands back the cut
    # dressed up as a threshold.
    if p_first > p.min() - margin:
        raise SystemExit(
            f"\nGUARD FAILED for {name!r}: the retained tail does not reach the "
            f"requested percentile.\n"
            f"  smallest retained plotting position : {p_first:.6f}\n"
            f"  lowest percentile requested         : {p.min():.6f}\n"
            f"  required margin                     : {margin:.6f}\n"
            f"  retained weighted fraction          : {1.0 - w_below / w_total:.4%}\n"
            f"Lower --cut-percentile and re-run. Refusing to return a clamped "
            f"threshold, which would look entirely plausible and be wrong.")

    out = np.interp(p, plotting, v_u, left=v_u[0], right=v_u[-1])
    # v and lat_index[order] are now in ascending value order, which is the
    # form the .npz archive wants: an EVT fit re-thresholds by slicing a
    # sorted array, not by re-sorting 3e8 values.
    return out, v.astype(np.float32), lat_index[order], p_first


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="Exact percentile calibration from retained tails.")
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--derived-subdir", default="derived/global",
                    help="relative to --base; use derived/global_sub48 to "
                         "validate against the existing 48-day thresholds")
    ap.add_argument("--cut-percentile", type=float, default=88.0,
                    help="cut placed here; ~12%% retained, covering p97..p99.9 "
                         "plus Lee's p95 plus room to vary an EVT threshold")
    ap.add_argument("--cut-sample-stride", type=int, default=7,
                    help="every Nth timestep when locating the cut. MUST NOT be "
                         "a multiple of 8 -- see the check below")
    ap.add_argument("--time-chunk", type=int, default=64,
                    help="timesteps per block in the full pass; sets peak RSS")
    ap.add_argument("--margin", type=float, default=0.005,
                    help="required headroom below the lowest percentile")
    ap.add_argument("--no-archive", action="store_true",
                    help="compute thresholds but do not write the .npz tails")
    ap.add_argument("--compress-archive", action="store_true",
                    help="zlib the .npz tails. OFF by default: the payload is "
                         "float32 mantissas from a continuous distribution plus "
                         "a value-sorted (hence scrambled) int16 latitude index, "
                         "so it compresses poorly -- a few GB saved for tens of "
                         "minutes of CPU across the 21.")
    ap.add_argument("--validate-against", default=None,
                    help="path to an existing thresholds_*.json; compare and "
                         "exit non-zero on disagreement")
    ap.add_argument("--tolerance", type=float, default=1e-6,
                    help="max relative difference accepted in --validate-against; "
                         "see the module docstring for why this is not zero")
    args = ap.parse_args()

    if args.cut_sample_stride % 8 == 0:
        # The data is 3-hourly: 8 steps per day. A stride that is a multiple of
        # 8 samples ONE time of day and nothing else, so the cut is set from a
        # diurnally phase-locked slice of the distribution. Harmless for a
        # threshold (the cut only has to sit below the lowest percentile, and
        # the guard checks that) but it silently makes the retained fraction
        # wrong for anything with a diurnal cycle.
        #
        # Measured, run 1098461 at stride 8 on the 48-day stores: 20 of 21
        # diagnostics retained 11.8-12.2%, f2d retained 11.19%. On the 4-day
        # files that stride lands on steps 0, 8, 16, 24 -- the file start plus
        # the first step after each of the three 8-day gaps, i.e. exactly the
        # timesteps where f2d's time derivative is one-sided or straddles a
        # gap. 100% of its cut sample was drawn from the pathological steps.
        raise SystemExit(
            f"--cut-sample-stride {args.cut_sample_stride} is a multiple of 8, "
            f"the number of 3-hourly steps per day.\nThat samples a single time "
            f"of day. Use a stride coprime with 8 (7, 9, 11, ...) so the "
            f"sub-sample walks the diurnal cycle.")

    base = Path(args.base)
    paths = month_paths(base, args.derived_subdir)
    diag = _load("diagnostics_module", "2_diagnostics.py")
    reference_table = diag.REFERENCE_TABLE

    ds0 = xr.open_zarr(paths[0])
    names = sorted([str(v) for v in ds0.data_vars],
                   key=lambda k: reference_table[k]["num"])
    w_lat = lat_weights_from(ds0["latitude"].values)
    pcts = np.asarray([SEVERITIES[s] for s in SEVERITIES], dtype=float)

    print(f"=== tail_thresholds on {len(paths)} stores under "
          f"{args.derived_subdir} ===")
    print(f"    {len(names)} diagnostics, {w_lat.size} latitude rows, "
          f"cut at p{args.cut_percentile}")
    print(f"    total timesteps: "
          f"{sum(xr.open_zarr(p).sizes['time'] for p in paths):,}")

    stamp = dt.date.today().isoformat()
    archive = base / "calibration" / f"tails_{stamp}"
    if not args.no_archive:
        archive.mkdir(parents=True, exist_ok=True)

    thresholds: dict[str, dict[str, float]] = {}
    signs: dict[str, str] = {}
    sample_sizes: dict[str, int] = {}
    t0 = time.time()

    for k, name in enumerate(names, 1):
        t_d = time.time()
        cut = cut_for(name, paths, w_lat, args.cut_sample_stride,
                      args.cut_percentile)
        values, lat_index, w_below, w_total, n_finite = scan_tail(
            name, paths, w_lat, cut, args.time_chunk)
        out, v_sorted, i_sorted, p_first = thresholds_from_tail(
            values, lat_index, w_lat, w_below, w_total, pcts, name,
            args.margin)

        thresholds[name] = {s: float(v) for s, v in zip(SEVERITIES, out)}
        signs[name] = reference_table[name]["sign"]
        sample_sizes[name] = int(n_finite)

        if not args.no_archive:
            writer = np.savez_compressed if args.compress_archive else np.savez
            writer(
                archive / f"{name}.npz",
                values=v_sorted, lat_index=i_sorted,
                w_below=np.float64(w_below), w_total=np.float64(w_total),
                cut=np.float64(cut), n_finite=np.int64(n_finite))

        kept = 1.0 - w_below / w_total
        flag = ""
        if kept > 0.40:
            # Almost always a clipped diagnostic (nva, ncsu1 apply max(.,0)):
            # if more than (100 - cut_percentile)% of the field is EXACTLY
            # zero, the cut lands on that tie mass and `>= cut` retains the
            # whole field. Harmless here, fatal on the full year -- it is a
            # 3e9-element sort instead of a 3.6e8 one.
            flag = "  <-- RETENTION HIGH: cut probably landed on a tie mass"
        print(f"  [{k:2d}/{len(names)}] {name:24s} cut={cut: .6g}  "
              f"kept={kept:6.2%}  n={n_finite:,}  p1={p_first:.5f}  "
              f"{time.time() - t_d:5.1f}s  rss={peak_rss_gb():.1f}G{flag}")

    print(f"    all diagnostics in {(time.time() - t0) / 60:.1f} min, "
          f"peak RSS {peak_rss_gb():.1f} GB")
    print(f"    (a full year is 7.62x this data volume: scale time roughly "
          f"linearly, and the per-diagnostic tail with it)")

    if args.validate_against:
        ref_path = Path(args.validate_against)
        if not ref_path.is_absolute():
            ref_path = base / ref_path
        ref, _, prov = calibration.load_thresholds(ref_path)
        worst, worst_where = 0.0, ""
        for name in thresholds:
            if name not in ref:
                print(f"  !! {name} absent from {ref_path.name}")
                return 1
            for sev, got in thresholds[name].items():
                want = ref[name][sev]
                denom = max(abs(want), 1e-300)
                rel = abs(got - want) / denom
                if rel > worst:
                    worst, worst_where = rel, f"{name}/{sev}"
        print(f"\nVALIDATION vs {ref_path.name} "
              f"(domain={prov.get('calibration_domain')}, "
              f"period={prov.get('period')})")
        print(f"  worst relative difference: {worst:.3e}  at {worst_where}")
        print(f"  tolerance:                 {args.tolerance:.3e}")
        if worst <= args.tolerance:
            print("  PASS -- the tail path reproduces the whole-array "
                  "calibration to floating-point summation order.")
            return 0
        print("  FAIL -- do NOT run the full year until this is understood.")
        return 1

    path = calibration.save_thresholds(
        base / "calibration" / f"thresholds_{stamp}.json",
        thresholds, signs, SEVERITIES,
        domain="global",
        period=f"{args.derived_subdir} ({len(paths)} stores)",
        pressure_levels=["175", "200", "225"],
        sample_sizes=sample_sizes,
        notes=(f"Computed by ada/tail_thresholds.py from retained tails "
               f"(cut p{args.cut_percentile}). Exact reproduction of "
               f"calib_weighted_percentile on the full field; the discarded "
               f"bulk enters only through W_total and W_below. Evaluated at "
               f"200 hPa from the 175/200/225 stencil."))
    print(f"\n    saved -> {path}")
    if not args.no_archive:
        print(f"    tails -> {archive}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
