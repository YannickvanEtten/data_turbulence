#!/usr/bin/env python
"""
ada/check_composition.py
========================
Pointwise composition identities on a derived zarr. FORMULA_AUDIT.md §7.3.

    pixi run python ada/check_composition.py <diagnostics.zarr> [--max-points N]

WHAT THIS IS FOR
----------------
CALIBRATION_REFERENCE.md §5.1 tests TI1 = VWS x DEF by comparing MEDIAN RATIOS
and gets agreement to 0.15 %. That is a good smell test but it can only ever be
approximate, because the median of a product is not the product of the medians.
Every one of these diagnostics is an exact algebraic function of the others, so
the saved fields permit the exact version: cell by cell, timestep by timestep.

Any deviation beyond float32 rounding means a mis-wired argument, a level
mismatch, or a chunk-boundary error -- the class of fault that produces a
plausible-looking field and is invisible in every distributional check.

Minutes on one month, and it verifies nine of the 21 against each other using
no external reference at all.

THE FIVE IDENTITIES
-------------------
  I1  DEF three independent ways: ti1/vws, ngm1/wind_speed, and the saved
      `deformation` field itself. Also decides the DEF-vs-DEF^2 question
      (FORMULA_AUDIT.md §5) empirically rather than by reading the source.
  I2  brown2 / (brown1 * vws^2) = 1/24                        (Sharman A14)
  I3  |(ti1 - ti2) / vws| = horizontal_divergence             (A15, A16, A33)
  I4  (zeta_a - f)^2 = vorticity_squared, zeta_a recovered from brown1 and DEF
      via A13. Ties brown1, deformation, vorticity_squared and latitude
      together in one expression.
  I5  |dT/dz| = ngm2 / DEF lands in a physically sensible band  (A29, sanity)

I2 is the one that closes `brown2`: FORMULA_AUDIT.md §2 confirms A14 against
the paper and §6 explains its published magnitude, so an exact pointwise
composition here leaves nothing open about it.

WHY TI2 IS CHECKED BACKWARDS (I3)
---------------------------------
TI2 = Sv (DEF - delta) needs the SIGNED divergence, and what reaches disk is
`horizontal_divergence` = |delta| (rojak takes the absolute value, A33). So the
identity is inverted: recover delta from the pair, delta = (ti1 - ti2)/vws,
then check its magnitude against the saved field. Same content, and it has the
side benefit of confirming that ti1 and ti2 were computed from the same vws.

READING THE OUTPUT
------------------
Errors are relative and should sit at the float32 floor, ~1e-7 to 1e-5. A
median error near 1e-3 or above is a real disagreement. The p99 column matters
as much as the median: an identity that holds in the bulk and fails in the tail
is exactly the failure mode a percentile-calibrated pipeline would carry all
the way into the econometrics.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _checkutil import peak_rss_gb, rel_error, summarise_error  # noqa: E402

OMEGA = 7.292115e-05          # matches 2_diagnostics.OMEGA

# Relative-error thresholds. float32 carries ~7 decimal digits, and these
# quantities are products of three or four such numbers, so 1e-5 is a
# comfortable pass and 1e-3 is unambiguously a real disagreement.
TOL_PASS = 1e-5
TOL_FAIL = 1e-3

NEEDED = ["ti1", "ti2", "ngm1", "ngm2", "deformation", "vertical_wind_shear",
          "wind_speed", "brown1", "brown2", "vorticity_squared",
          "horizontal_divergence"]


def verdict(err_summary: dict) -> str:
    m = err_summary["median"]
    if not np.isfinite(m):
        return "NO DATA"
    if m <= TOL_PASS:
        return "PASS"
    return "CHECK" if m <= TOL_FAIL else "FAIL"


def flat(da: xr.DataArray, idx: np.ndarray | None) -> np.ndarray:
    a = np.asarray(da.values, dtype=np.float64).ravel()
    return a if idx is None else a[idx]


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("zarr", help="a derived diagnostics .zarr")
    ap.add_argument("--max-points", type=int, default=20_000_000,
                    help="randomly subsample to at most this many cells "
                         "(default 20M). The identities hold pointwise, so a "
                         "sample tests them exactly as well as the full field "
                         "and keeps peak memory near the file size.")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    path = Path(args.zarr)
    if not path.exists():
        print(f"!! not found: {path}")
        return 1

    print(f">>> {path}")
    ds = xr.open_zarr(path)

    missing = [k for k in NEEDED if k not in ds.data_vars]
    if missing:
        print(f"!! missing variables, cannot run: {missing}")
        return 1

    conv = ds.attrs.get("deformation_convention")
    print(f"    dims  : {dict(ds.sizes)}")
    print(f"    attrs : deformation_convention={conv!r} "
          f"f2d_variant={ds.attrs.get('f2d_variant')!r}")
    if conv is None:
        print("    (no provenance attributes -- written before 2026-08-29, so")
        print("     `deformation` is expected to hold DEF^2; I1 confirms below)")

    n_total = int(np.prod([ds.sizes[d] for d in ds["ti1"].dims]))
    idx = None
    if n_total > args.max_points:
        rng = np.random.default_rng(args.seed)
        idx = rng.choice(n_total, size=args.max_points, replace=False)
        print(f"    sampling {args.max_points:,} of {n_total:,} cells")
    else:
        print(f"    using all {n_total:,} cells")

    ti1 = flat(ds["ti1"], idx)
    ti2 = flat(ds["ti2"], idx)
    ngm1 = flat(ds["ngm1"], idx)
    ngm2 = flat(ds["ngm2"], idx)
    defo = flat(ds["deformation"], idx)
    vws = flat(ds["vertical_wind_shear"], idx)
    spd = flat(ds["wind_speed"], idx)
    br1 = flat(ds["brown1"], idx)
    br2 = flat(ds["brown2"], idx)
    z2 = flat(ds["vorticity_squared"], idx)
    hdiv = flat(ds["horizontal_divergence"], idx)

    lat2d = ds["latitude"].broadcast_like(ds["ti1"])
    lat = flat(lat2d, idx)
    f_cor = 2 * OMEGA * np.sin(np.deg2rad(lat))

    print(f"    [rss {peak_rss_gb():.1f} GB]\n")

    rows = []

    # ---------------------------------------------------------------- I1
    # DEF recovered two ways from products that never involved the DEF
    # diagnostic itself, then compared to each other and to the stored field.
    with np.errstate(divide="ignore", invalid="ignore"):
        def_from_ti1 = np.where(vws > 0, ti1 / vws, np.nan)
        def_from_ngm1 = np.where(spd > 0, ngm1 / spd, np.nan)

    e = summarise_error(rel_error(def_from_ti1, def_from_ngm1))
    rows.append(("I1a  ti1/vws  ==  ngm1/wind_speed", e, verdict(e)))

    e_unsq = summarise_error(rel_error(defo, def_from_ti1))
    e_sq = summarise_error(rel_error(defo, def_from_ti1 ** 2))
    stored_is_squared = (e_sq["median"] < e_unsq["median"])
    label = ("I1b  deformation == (ti1/vws)^2   [stored SQUARED]"
             if stored_is_squared else
             "I1b  deformation == ti1/vws       [stored un-squared]")
    rows.append((label, e_sq if stored_is_squared else e_unsq,
                 verdict(e_sq if stored_is_squared else e_unsq)))

    # The physical DEF, whichever convention the file uses. Everything below
    # that needs DEF uses this, so the script works on both old and new zarr.
    DEF = np.sqrt(np.abs(defo)) if stored_is_squared else defo

    # ---------------------------------------------------------------- I2
    with np.errstate(divide="ignore", invalid="ignore"):
        got = np.where((br1 > 0) & (vws > 0), br2 / (br1 * vws ** 2), np.nan)
    e = summarise_error(rel_error(got, np.full_like(got, 1.0 / 24.0)))
    rows.append(("I2   brown2/(brown1*vws^2) == 1/24        A14", e, verdict(e)))

    # ---------------------------------------------------------------- I3
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.where(vws > 0, (ti1 - ti2) / vws, np.nan)
    e = summarise_error(rel_error(np.abs(delta), hdiv))
    rows.append(("I3   |(ti1-ti2)/vws| == |divergence|  A15/A16", e, verdict(e)))

    # ---------------------------------------------------------------- I4
    # A13: brown1^2 = 0.3 zeta_a^2 + Dsh^2 + Dst^2, and Dsh^2 + Dst^2 = DEF^2.
    # zeta_a is recovered up to a sign, so take whichever root reproduces
    # vorticity_squared better -- the sign of absolute vorticity is not
    # knowable from these fields and is not what is being tested.
    inner = (br1 ** 2 - DEF ** 2) / 0.3
    zeta_a_mag = np.sqrt(np.maximum(inner, 0.0))
    z_pos = (zeta_a_mag - f_cor) ** 2
    z_neg = (-zeta_a_mag - f_cor) ** 2
    e_pos = rel_error(z_pos, z2)
    e_neg = rel_error(z_neg, z2)
    best = np.minimum(e_pos, e_neg) if e_pos.size == e_neg.size else e_pos
    e = summarise_error(best)
    rows.append(("I4   (zeta_a-f)^2 == vorticity_squared   A13/A21", e, verdict(e)))

    # ---------------------------------------------------------------- I5
    with np.errstate(divide="ignore", invalid="ignore"):
        dT_dz = np.where(DEF > 0, ngm2 / DEF, np.nan)
    good = dT_dz[np.isfinite(dT_dz)]

    # ------------------------------------------------------------- report
    print("=" * 78)
    print("POINTWISE COMPOSITION IDENTITIES")
    print("=" * 78)
    print(f"   {'identity':<52}{'median':>10}{'p99':>10}{'':>8}")
    for label, e, v in rows:
        print(f"   {label:<52}{e['median']:>10.2e}{e['p99']:>10.2e}  {v}")
    print()
    print(f"   pass <= {TOL_PASS:.0e}   check <= {TOL_FAIL:.0e}   fail above")
    print(f"   n compared: {rows[0][1]['n']:,}")

    print("\n" + "-" * 78)
    print("I5  |dT/dz| = ngm2/DEF  — physical sanity, not an identity")
    if good.size:
        q = np.quantile(good, [0.01, 0.5, 0.99])
        print(f"   p01 {q[0]:.3e}   median {q[1]:.3e}   p99 {q[2]:.3e}  K/m")
        print("   Expect a median of order 1e-3 K/m near the tropopause. Orders")
        print("   away means ngm2 and deformation disagree about which field")
        print("   they multiplied.")
    else:
        print("   no finite values")

    print("\n" + "-" * 78)
    if stored_is_squared:
        print("   `deformation` in this file holds DEF^2 (pre-2026-08-29).")
        print("   Exceedance fields are unaffected (STATUS.md §5.5); magnitudes")
        print("   are. Re-derive before any phase-5 magnitude work — see")
        print("   FORMULA_AUDIT.md §5.")
    else:
        print("   `deformation` in this file holds DEF, un-squared. Good.")

    worst = max((e["median"] for _, e, _ in rows if np.isfinite(e["median"])),
                default=float("nan"))
    print(f"\n   PEAK RSS {peak_rss_gb():.1f} GB")
    if not np.isfinite(worst):
        return 1
    return 0 if worst <= TOL_FAIL else 2


if __name__ == "__main__":
    raise SystemExit(main())
