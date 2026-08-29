#!/usr/bin/env python
"""
ada/check_ubf.py
================
Why does UBF behave unlike its 20 siblings? FORMULA_AUDIT.md §7.2.

    pixi run python ada/check_ubf.py <raw.grib> [--box prosser]

THE ANOMALY
-----------
STATUS.md §12.6 measured the MOG exceedance change for each of the 21 over
DJF 1979 -> 2020, and found a coherent picture with exactly one exception:

    ubf   0.249 % -> 0.276 %   +11 %      15th of 21 on level, LAST on trend

while every sibling ran +38 % to +226 %, and while the trend magnitude
otherwise anti-correlates neatly with the exceedance level.

Prosser's Figure S5 makes it worse rather than better. Its "magnitude of
residual of nonlinear balance equation" panel is strongly POSITIVE over the
North Atlantic and Europe -- UBF is supposed to be near the TOP of the 21 over
exactly the box we measure, not the bottom.

WHAT IS ALREADY RULED OUT
-------------------------
The equation. FORMULA_AUDIT.md §2 confirms our implementation against Sharman
A30 term by term:

    UBF = -grad^2(Phi) + 2 J(u,v) + f zeta - beta u                  (A30)

we compute the negative of that and take the absolute value, which is
identical. STATUS.md §4g separately fixed the two geometry errors and pinned
them with theorems. So the formula is not the fault, and this script looks at
the three things that remain.

WHAT THIS SCRIPT MEASURES
-------------------------
1. TERM MAGNITUDES AND CANCELLATION. UBF is a residual of four large terms.
   If the residual is orders smaller than the terms, then its value is set by
   how precisely they cancel, and anything that perturbs a term at the
   1e-3 level dominates the answer. Prints each term's weighted median and
   the ratio residual / largest-term.

2. FLOAT32. The pipeline computes and stores float32
   (ada/diagnostics_global.py). Seven significant digits is ample for a
   magnitude and can be nowhere near enough for a residual that cancels to
   three or four. Recomputes in float64 and compares by RANK -- which is all
   the exceedance counting uses -- and by the actual flagged set at p99.6.

3. ERA5 vorticity vs DERIVED vorticity. The f*zeta term uses ERA5's own `vo`
   product, while the Laplacian and the Jacobian are differentiated from u and
   v on our grid. Those are two different estimates of the same physical
   field. Any systematic disagreement between them does not cancel: it lands
   directly in the residual and is then indistinguishable from atmospheric
   imbalance. This is the most likely culprit and it is measured here against
   the size of the residual itself, which is the only comparison that matters.

READING IT
----------
A cancellation ratio near 1e-3 with a float32-vs-float64 rank correlation
below ~0.99 means precision is a real contributor and the production run
should compute UBF in float64 (storing float32 is fine -- it is the arithmetic
that matters, not the storage).

A vorticity discrepancy whose f*(zeta_era5 - zeta_derived) is comparable to
the residual means UBF is substantially measuring a product-versus-derivative
inconsistency rather than imbalance, and the honest fix is to compute the
f*zeta term from the SAME derived vorticity as the other terms -- consistency
mattering more here than which estimate is individually better.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _checkutil import (PROSSER_BOX, WILLIAMS_BOX, cos_phi_weights,  # noqa: E402
                        load_module, peak_rss_gb, spearman, subset_box)
from calib_weighted_percentile import weighted_percentile  # noqa: E402

MOG_PERCENTILE = 99.6


def ubf_terms(diag, ds, level: int, dtype) -> dict:
    """The four A30 terms, computed exactly as 2_diagnostics.ubf() does.

    Deliberately a COPY of that function's body rather than a call to it,
    because the whole point is to see the terms separately and at a chosen
    precision. Any edit to ubf() must be mirrored here; the residual printed
    below is checked against ubf()'s own output so the two cannot drift
    silently.
    """
    from rojak.core.derivatives import VelocityDerivative, vector_derivatives

    def cast(da):
        return da.astype(dtype)

    u = cast(diag._sel_level(ds, "eastward_wind", level))
    v = cast(diag._sel_level(ds, "northward_wind", level))
    phi = cast(diag._sel_level(ds, "geopotential", level))
    zeta = cast(diag._sel_level(ds, "vorticity", level))

    lat_rad = np.deg2rad(ds["latitude"])
    f = (2 * diag.OMEGA * np.sin(lat_rad)).broadcast_like(u.isel(time=0))
    beta = (2 * diag.OMEGA * np.cos(lat_rad) / diag.R_EARTH).broadcast_like(u.isel(time=0))

    dphi_dx, dphi_dy = diag._grad(phi)
    vd_phi = vector_derivatives(
        dphi_dx, dphi_dy, "deg",
        components=[VelocityDerivative.DU_DX, VelocityDerivative.DV_DY])
    d2phi = vd_phi[VelocityDerivative.DU_DX] + vd_phi[VelocityDerivative.DV_DY]

    vd = vector_derivatives(
        u, v, "deg",
        components=[VelocityDerivative.DU_DX, VelocityDerivative.DU_DY,
                    VelocityDerivative.DV_DX, VelocityDerivative.DV_DY])
    du_dx = vd[VelocityDerivative.DU_DX]
    du_dy = vd[VelocityDerivative.DU_DY]
    dv_dx = vd[VelocityDerivative.DV_DX]
    dv_dy = vd[VelocityDerivative.DV_DY]
    jac = du_dx * dv_dy - du_dy * dv_dx

    return {
        "grad2_phi": d2phi,
        "-2J": -2 * jac,
        "-f*zeta": -f * zeta,
        "+beta*u": beta * u,
        "_zeta_era5": zeta,
        "_zeta_derived": dv_dx - du_dy,     # the same vorticity, differentiated
        "_f": f,
    }


def wmed(da, w) -> float:
    v = np.asarray(da.values if hasattr(da, "values") else da,
                   dtype=np.float64).ravel()
    ok = np.isfinite(v)
    if ok.sum() < 100:
        return float("nan")
    return float(weighted_percentile(v[ok], w[ok], 50.0))


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("grib", help="raw ERA5 GRIB file")
    ap.add_argument("--box", default="prosser",
                    choices=["prosser", "williams", "none"])
    ap.add_argument("--level", type=int, default=200)
    args = ap.parse_args()

    box = {"williams": WILLIAMS_BOX, "prosser": PROSSER_BOX, "none": None}[args.box]
    path = Path(args.grib)
    if not path.exists():
        print(f"!! not found: {path}")
        return 1

    diag = load_module("diagnostics", "2_diagnostics.py")
    print(f">>> {path.name}  ({path.stat().st_size / 1e9:.2f} GB), box {args.box}")
    ds = diag.prepare_for_rojak(diag.load_era5(path))._dataset

    t32 = ubf_terms(diag, ds, args.level, np.float32)
    t64 = ubf_terms(diag, ds, args.level, np.float64)
    print(f"    [rss {peak_rss_gb():.1f} GB]")

    term_keys = ["grad2_phi", "-2J", "-f*zeta", "+beta*u"]

    def clip(da):
        if box is None:
            return da
        return subset_box(da.to_dataset(name="x"), **box)["x"]

    t32 = {k: clip(v) for k, v in t32.items()}
    t64 = {k: clip(v) for k, v in t64.items()}

    w = cos_phi_weights(t64["grad2_phi"])

    res32 = sum(t32[k] for k in term_keys)
    res64 = sum(t64[k] for k in term_keys)

    # Cross-check against the production function, so this script cannot
    # quietly diverge from the thing it is diagnosing.
    ref = clip(diag.ubf(ds, target_level=args.level))
    ref_v = np.asarray(ref.values, dtype=np.float64).ravel()
    mine_v = np.abs(np.asarray(res32.values, dtype=np.float64).ravel())
    agree = spearman(ref_v, mine_v)
    print(f"    sanity: rho(this decomposition, 2_diagnostics.ubf) = {agree:.6f}")
    if not (agree > 0.999):
        print("    !! decomposition does NOT match ubf() — the copy in "
              "ubf_terms() has drifted from 2_diagnostics.ubf(). Fix that "
              "before reading anything below.")
        return 2

    # ------------------------------------------------------------ 1. terms
    print("\n" + "=" * 78)
    print("1. TERM MAGNITUDES AND CANCELLATION")
    print("=" * 78)
    print(f"   {'term':<14}{'weighted median |term|':>26}")
    mags = {}
    for k in term_keys:
        m = wmed(np.abs(t64[k]), w)
        mags[k] = m
        print(f"   {k:<14}{m:>26.4e}")
    res_med = wmed(np.abs(res64), w)
    biggest = max(v for v in mags.values() if np.isfinite(v))
    print(f"   {'|residual|':<14}{res_med:>26.4e}")
    print(f"\n   cancellation ratio  |residual| / largest term = "
          f"{res_med / biggest:.3e}")
    print("   Below ~1e-3 means the answer is set by how precisely four large")
    print("   terms cancel, and every source of small error is amplified by")
    print("   the reciprocal of this number.")

    # -------------------------------------------------------- 2. precision
    print("\n" + "=" * 78)
    print("2. FLOAT32 vs FLOAT64")
    print("=" * 78)
    a32 = np.abs(np.asarray(res32.values, dtype=np.float64).ravel())
    a64 = np.abs(np.asarray(res64.values, dtype=np.float64).ravel())
    rho = spearman(a32, a64)
    ok = np.isfinite(a32) & np.isfinite(a64)
    thr32 = weighted_percentile(a32[ok], w[ok], MOG_PERCENTILE)
    thr64 = weighted_percentile(a64[ok], w[ok], MOG_PERCENTILE)
    m32, m64 = a32 >= thr32, a64 >= thr64
    inter = np.count_nonzero(m32 & m64 & ok)
    union = np.count_nonzero((m32 | m64) & ok)
    print(f"   Spearman rho(float32, float64)          {rho:.6f}")
    print(f"   MOG threshold  float32 {float(thr32):.4e}   float64 {float(thr64):.4e}")
    print(f"   flagged-set overlap at p{MOG_PERCENTILE}            "
          f"{inter / union if union else float('nan'):.4f}")
    print("\n   Exceedance counting sees RANKS, so rho is the number that")
    print("   matters. Below ~0.99, or an overlap below ~0.95, means float32")
    print("   arithmetic is materially reordering the tail and the production")
    print("   run should compute this diagnostic in float64. Storage can stay")
    print("   float32 — it is the arithmetic, not the file.")

    # ------------------------------------------------------- 3. vorticity
    print("\n" + "=" * 78)
    print("3. ERA5 vorticity vs DERIVED vorticity")
    print("=" * 78)
    z_era = t64["_zeta_era5"]
    z_der = t64["_zeta_derived"]
    f = t64["_f"]
    diff_term = np.abs(f * (z_era - z_der))
    rho_z = spearman(np.asarray(z_era.values).ravel(),
                     np.asarray(z_der.values).ravel())
    print(f"   median |zeta_era5|                      {wmed(np.abs(z_era), w):.4e}")
    print(f"   median |zeta_derived|                   {wmed(np.abs(z_der), w):.4e}")
    print(f"   median |zeta_era5 - zeta_derived|       {wmed(np.abs(z_era - z_der), w):.4e}")
    print(f"   Spearman rho                            {rho_z:.6f}")
    print(f"\n   median |f * (zeta_era5 - zeta_derived)| {wmed(diff_term, w):.4e}")
    print(f"   median |UBF residual|                   {res_med:.4e}")
    ratio = wmed(diff_term, w) / res_med if res_med else float("nan")
    print(f"   ratio                                   {ratio:.3f}")
    print("\n   This is the key number. UBF's f*zeta term uses ERA5's own `vo`")
    print("   product while the Laplacian and Jacobian are differentiated from")
    print("   u and v on our grid — two estimates of the same field. Their")
    print("   disagreement does not cancel; it lands in the residual. A ratio")
    print("   approaching or exceeding 1 means UBF is substantially measuring")
    print("   that inconsistency rather than atmospheric imbalance, which")
    print("   would explain both the flat trend in STATUS.md §12.6 and the")
    print("   disagreement with Prosser Figure S5.")
    print("\n   If it does: compute f*zeta from the DERIVED vorticity instead,")
    print("   so all four terms share one geometry and one estimator. That is")
    print("   the same argument STATUS.md §4g used to fix the Laplacian —")
    print("   consistency across the terms of a residual matters more than the")
    print("   individual accuracy of any one of them.")

    print(f"\n   PEAK RSS {peak_rss_gb():.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
