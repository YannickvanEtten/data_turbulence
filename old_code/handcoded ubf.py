"""
handcoded_ubf.py
=================
Handcoded reference implementation of UBF (Unbalanced Flow diagnostic)
from Sharman et al. (2006), *Wea. Forecasting*, Eq. A30 — the
magnitude of the residual of the nonlinear balance equation (NBE).

Charney's nonlinear balance equation (Koch & Caracena 2002, following
Charney 1955):

    ∇²Φ = 2·J(u,v) + f·ζ − β·u

so the NBE residual is:

    R = ∇²Φ − 2·J(u,v) − f·ζ + β·u
      = −(−∇²Φ + 2·J(u,v) + f·ζ − β·u)      (Sharman A30 form)

Both forms give the same |R|.  We compute the magnitude:

    UBF = |R|

with:
    ∇²Φ   Laplacian of geopotential          m/s²/m = s⁻²
    J(u,v) Jacobian = ∂u/∂x·∂v/∂y − ∂u/∂y·∂v/∂x    s⁻²
    f     Coriolis parameter = 2Ω sin(φ)     s⁻¹
    ζ     Relative vorticity                 s⁻¹
    β     df/dy = 2Ω cos(φ) / R_earth        s⁻¹ m⁻¹
    u     zonal wind                         m s⁻¹

Every term is s⁻².  Typical mid-latitude synoptic-scale values:
    ∇²Φ, f·ζ ~ 1e-9 s⁻²      (large)
    2·J     ~ 1e-10 s⁻²      (small)
    β·u     ~ 1e-10 s⁻²      (small)
    residual ~ 1e-11 – 1e-10 s⁻²    (the NBE-balance residual is a
                                     partial cancellation of the two
                                     large terms).

Williams & Joshi (2013) Table 1 gives DJF median = 161 × 10⁻¹² s⁻² ≈
1.6 × 10⁻¹⁰ s⁻², perfectly consistent with that scale.

Rojak's `UBF._compute` has a sign-flip bug: it returns
    |∇²Φ + 2·J + f·ζ − β·u|
i.e. the sum of the two large terms rather than their residual, which
gives values ~1e-3 s⁻² — 7 orders of magnitude too large.  This script
demonstrates the fix.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# rojak helpers we still use, so the comparison is apples-to-apples
from rojak.core.derivatives import (
    grid_spacing, first_derivative,
)
# rojak's UBF implementation, to compare against
from rojak.core.data import CATData
from rojak.turbulence.diagnostic import DiagnosticFactory
from rojak.orchestrator.configuration import TurbulenceDiagnostics


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
OMEGA = 7.2921159e-5  # Earth angular velocity   s⁻¹
R_EARTH = 6.371e6     # Earth radius             m


# ---------------------------------------------------------------------------
# Helpers — geospatial derivatives via rojak's grid_spacing
# ---------------------------------------------------------------------------
def _spacings_1d(ds: xr.Dataset):
    """Return 1-D dx (nlon-1) and dy (nlat-1) arrays in metres."""
    gs = grid_spacing(ds["latitude"], ds["longitude"], units="degrees")
    dx = np.asarray(gs.dx)
    dy = np.asarray(gs.dy)
    dx_1d = dx.mean(axis=0) if dx.ndim == 2 else dx
    dy_1d = dy.mean(axis=1) if dy.ndim == 2 else dy
    return dx_1d, dy_1d


def _ddx(field: xr.DataArray, dx_1d) -> xr.DataArray:
    return first_derivative(field, dx_1d, axis=field.dims.index("longitude"))


def _ddy(field: xr.DataArray, dy_1d) -> xr.DataArray:
    return first_derivative(field, dy_1d, axis=field.dims.index("latitude"))


# ---------------------------------------------------------------------------
# UBF terms — one function per term, so the physics is auditable
# ---------------------------------------------------------------------------
def compute_ubf_terms(ds: xr.Dataset, target_level: int = 200) -> dict[str, xr.DataArray]:
    """Compute all four NBE terms at `target_level` (hPa). Returns a dict."""
    u   = ds["eastward_wind"    ].sel(pressure_level=target_level)
    v   = ds["northward_wind"   ].sel(pressure_level=target_level)
    phi = ds["geopotential"     ].sel(pressure_level=target_level)
    zeta= ds["vorticity"        ].sel(pressure_level=target_level)

    # Coriolis f and its meridional derivative β
    lat_rad = np.deg2rad(ds["latitude"])
    f = 2 * OMEGA * np.sin(lat_rad)                       # (nlat,)
    beta = 2 * OMEGA * np.cos(lat_rad) / R_EARTH          # (nlat,)
    # Broadcast to (lat, lon)
    f_2d    = f.broadcast_like(u.isel(time=0))
    beta_2d = beta.broadcast_like(u.isel(time=0))

    # --- ∇²Φ  — TRUE second-derivative Laplacian, not rojak's broken helper ---
    # (rojak's `spatial_laplacian` = `divergence(∂Φ/∂x, ∂Φ/∂y)` = ∂Φ/∂x + ∂Φ/∂y,
    # NOT ∂²Φ/∂x² + ∂²Φ/∂y² — it just sums the first derivatives.)
    dx_1d, dy_1d = _spacings_1d(ds)
    dphi_dx  = _ddx(phi, dx_1d)
    dphi_dy  = _ddy(phi, dy_1d)
    d2phi_dx2 = _ddx(dphi_dx, dx_1d)
    d2phi_dy2 = _ddy(dphi_dy, dy_1d)
    laplacian_phi = d2phi_dx2 + d2phi_dy2

    # --- Jacobian J(u,v) = ∂u/∂x · ∂v/∂y − ∂u/∂y · ∂v/∂x ---
    du_dx = _ddx(u, dx_1d);  du_dy = _ddy(u, dy_1d)
    dv_dx = _ddx(v, dx_1d);  dv_dy = _ddy(v, dy_1d)
    jacobian = du_dx * dv_dy - du_dy * dv_dx

    # --- Coriolis-vorticity term f·ζ ---
    f_zeta = f_2d * zeta

    # --- β·u ---
    beta_u = beta_2d * u

    return {
        "laplacian_phi": laplacian_phi,          # ∇²Φ
        "jacobian":       jacobian,               # J(u,v)
        "two_jacobian":   2 * jacobian,           # 2·J(u,v)
        "f_zeta":         f_zeta,                 # f·ζ
        "beta_u":         beta_u,                 # β·u
    }


def compute_ubf_handcoded(ds: xr.Dataset, target_level: int = 200) -> xr.DataArray:
    """
    Correct UBF per Sharman (2006) A30 / Koch & Caracena (2002):
        UBF = | ∇²Φ − 2·J(u,v) − f·ζ + β·u |
    """
    terms = compute_ubf_terms(ds, target_level)
    residual = (terms["laplacian_phi"]
                - terms["two_jacobian"]
                - terms["f_zeta"]
                + terms["beta_u"])
    ubf = np.abs(residual).rename("ubf_handcoded")
    ubf.attrs["long_name"] = "UBF (handcoded, Sharman A30, Charney NBE residual)"
    ubf.attrs["units"] = "s^-2"
    return ubf


# ---------------------------------------------------------------------------
# Wrapper: run both rojak and handcoded on one GRIB file, produce comparison
# ---------------------------------------------------------------------------
def compare(input_path: str | Path, out_dir: str | Path = "/home/claude/work/ubf_out"):
    from old_code.cat_pipeline import load_era5, prepare_for_rojak  # reuse pipeline plumbing
    input_path = Path(input_path); out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    ds_raw = load_era5(input_path)
    catdata = prepare_for_rojak(ds_raw)
    ds_rojak = catdata._dataset

    # --- rojak's UBF ---
    factory = DiagnosticFactory(catdata)
    ubf_rojak = factory.create(TurbulenceDiagnostics.UBF).computed_value.sel(pressure_level=200)

    # --- handcoded UBF ---
    ubf_hand = compute_ubf_handcoded(ds_rojak, target_level=200)

    # --- individual terms (for a magnitude audit) ---
    terms = compute_ubf_terms(ds_rojak, target_level=200)

    # 1) Print term-by-term magnitudes (day-mean absolute median)
    print("\n=== UBF term magnitudes (median |·| over full day) ===")
    print("Sharman (2006) A30 / Koch & Caracena (2002):")
    print("    NBE:      ∇²Φ = 2J(u,v) + fζ − βu")
    print("    Residual: R = ∇²Φ − 2J − fζ + βu")
    print("    UBF = |R|")
    print()
    for name, da in terms.items():
        arr = np.abs(np.asarray(da).ravel())
        arr = arr[np.isfinite(arr)]
        print(f"    {name:20s} median|·| = {np.median(arr):.3e}  s⁻²")
    print()
    print("Expected residual magnitude ≈ 1e-11 to 1e-10 s⁻²  (W&J: 1.6e-10)")

    # 2) Compare rojak vs handcoded
    rojak_arr = np.abs(np.asarray(ubf_rojak)).ravel()
    hand_arr  = np.abs(np.asarray(ubf_hand)).ravel()
    rojak_arr = rojak_arr[np.isfinite(rojak_arr)]
    hand_arr  = hand_arr[np.isfinite(hand_arr)]

    print("\n=== Comparison: rojak's UBF vs handcoded UBF ===")
    print(f"    rojak     median = {np.median(rojak_arr):.3e}  p99 = {np.quantile(rojak_arr, 0.99):.3e}")
    print(f"    handcoded median = {np.median(hand_arr):.3e}  p99 = {np.quantile(hand_arr, 0.99):.3e}")
    ratio = np.median(rojak_arr) / max(np.median(hand_arr), 1e-30)
    print(f"    rojak / handcoded ratio (median) = {ratio:.3e}")
    print(f"    W&J Table 1 median = 1.61e-10  s⁻²  (161 × 10⁻¹² s⁻²)")

    # 3) Rank correlation between rojak and handcoded (day mean, over grid)
    r_flat = np.abs(ubf_rojak.mean("time")).values.ravel()
    h_flat = np.abs(ubf_hand.mean("time")).values.ravel()
    ok = np.isfinite(r_flat) & np.isfinite(h_flat)
    rho, _ = spearmanr(r_flat[ok], h_flat[ok])
    print(f"\n    Spearman rank correlation ρ(rojak, handcoded) = {rho:+.3f}")

    # 4) Compare handcoded against the working shear cluster
    # Reload the previously computed diagnostics (from cat_pipeline)
    print("\n=== Rank correlation of handcoded UBF with working shear diagnostics ===")
    diag_nc = out_dir.parent / "diagnostics.nc"
    if diag_nc.exists():
        old = xr.open_dataset(diag_nc)
        for shear_key in ["vertical_wind_shear", "ti1", "ti2", "ngm1"]:
            s = np.abs(old[shear_key].mean("time")).values.ravel()
            h = h_flat
            m = np.isfinite(s) & np.isfinite(h)
            r_rojak_sk = spearmanr(r_flat[m], s[m])[0]
            r_hand_sk  = spearmanr(h[m], s[m])[0]
            print(f"    ρ(handcoded, {shear_key:22s}) = {r_hand_sk:+.3f}    (rojak's was {r_rojak_sk:+.3f})")

    # 5) Plots
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    day_mean_hand = np.abs(ubf_hand.mean("time"))
    day_mean_rojak = np.abs(ubf_rojak.mean("time"))
    p2, p98 = np.quantile(day_mean_hand, [0.02, 0.98])
    day_mean_hand.plot.pcolormesh(ax=axes[0], x="longitude", y="latitude",
                                  cmap="magma", vmin=p2, vmax=p98)
    axes[0].set_title(f"Handcoded UBF (day mean)\nmedian={np.median(hand_arr):.2e} s⁻²")
    p2, p98 = np.quantile(day_mean_rojak, [0.02, 0.98])
    day_mean_rojak.plot.pcolormesh(ax=axes[1], x="longitude", y="latitude",
                                   cmap="magma", vmin=p2, vmax=p98)
    axes[1].set_title(f"Rojak UBF (day mean)\nmedian={np.median(rojak_arr):.2e} s⁻²")
    # scatter
    idx = np.random.default_rng(0).choice(ok.sum(), size=min(5000, ok.sum()), replace=False)
    axes[2].loglog(h_flat[ok][idx], r_flat[ok][idx], ".", markersize=2, alpha=0.4)
    axes[2].set_xlabel("Handcoded UBF (s⁻²)")
    axes[2].set_ylabel("Rojak UBF (s⁻²)")
    axes[2].set_title(f"Scatter (Spearman ρ = {rho:+.3f})")
    axes[2].grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "ubf_comparison.png", dpi=140, bbox_inches="tight")
    plt.close()
    print(f"\nSaved plots: {out_dir / 'ubf_comparison.png'}")

    # 6) Save handcoded UBF
    ubf_hand.to_netcdf(out_dir / "ubf_handcoded.nc")
    print(f"Saved: {out_dir / 'ubf_handcoded.nc'}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python handcoded_ubf.py <path/to/era5.grib>")
        sys.exit(1)
    compare(sys.argv[1])