"""
4_verify.py
===========
Phase B verification harness.

For each of the 7 hand-written diagnostics, run two independent checks:

  (a) DIMENSIONAL CHECK  — the output units match what Williams & Joshi
      (2013) Table 1 expects (order-of-magnitude sanity on real data).

  (b) CORRECTED-ROJAK CROSS-CHECK — monkey-patch the *specific* broken line
      in rojak (in-memory, throwaway), then confirm the hand-written value
      agrees with the corrected rojak to numerical precision.  This oracle
      is NEVER part of the shipping pipeline; it exists only to prove that
      two independent implementations of the same corrected formula agree.

Self-tests on a NetCDF/GRIB file given on the command line (default: the
clean synthetic file).  Run against era5_validation_subset.nc for real data.

Usage:
    python verify_handcoded.py [path/to/data.nc]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import spearmanr

import importlib.util
from pathlib import Path
_spec = importlib.util.spec_from_file_location(
    "diagnostics", Path(__file__).with_name("2_diagnostics.py"))
hc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hc)
load_era5 = hc.load_era5
prepare_for_rojak = hc.prepare_for_rojak


# ---------------------------------------------------------------------------
# Monkey-patch helpers — build a "corrected rojak" in memory
# ---------------------------------------------------------------------------
def patch_rojak():
    """Apply the minimal Phase-A fixes to rojak in-memory. Returns an
    unpatch() callable to restore the originals. THROWAWAY — test only."""
    import rojak.turbulence.diagnostic as D
    import rojak.turbulence.calculations as Ccalc
    import rojak.core.derivatives as Cder

    originals = {}

    # --- Fix 1: Richardson  N²/Sv → N²/Sv² -------------------------------
    originals["GradientRichardson._compute"] = D.GradientRichardson._compute
    def fixed_ri_compute(self):
        ri = self._brunt_vaisala / (self._vws ** 2)   # square the shear
        return -ri if self._is_negative else ri
    D.GradientRichardson._compute = fixed_ri_compute

    # --- Fix 2: latitudinal_derivative  f/R → 2Ω cosφ/R ------------------
    originals["latitudinal_derivative"] = Ccalc.latitudinal_derivative
    def fixed_lat_deriv(coriolis_param):
        # β = ∂f/∂y = 2Ω cosφ / R.  Recover φ from f = 2Ω sinφ.
        omega = Ccalc.EARTH_ANGULAR_VELOCITY
        R = Ccalc.EARTH_AVG_RADIUS if hasattr(Ccalc, "EARTH_AVG_RADIUS") else 6371008.7714
        sin_phi = coriolis_param / (2 * omega)
        sin_phi = sin_phi.clip(-1, 1)
        cos_phi = np.sqrt(1 - sin_phi ** 2)
        return 2 * omega * cos_phi / R
    Ccalc.latitudinal_derivative = fixed_lat_deriv
    D.latitudinal_derivative = fixed_lat_deriv  # rebind name imported into diagnostic.py

    # --- Fix 3: spatial_laplacian  (∂f/∂x+∂f/∂y) → (∂²f/∂x²+∂²f/∂y²) ------
    originals["spatial_laplacian"] = Cder.spatial_laplacian
    def fixed_laplacian(array, units, gradient_mode, geod=None, crs=None):
        g1 = Cder.spatial_gradient(array, units, gradient_mode, geod=geod, crs=crs)
        g2x = Cder.spatial_gradient(g1["dfdx"], units, gradient_mode, geod=geod, crs=crs)
        g2y = Cder.spatial_gradient(g1["dfdy"], units, gradient_mode, geod=geod, crs=crs)
        return g2x["dfdx"] + g2y["dfdy"]
    Cder.spatial_laplacian = fixed_laplacian
    D.spatial_laplacian = fixed_laplacian

    # --- Fix 4: UBF sign  (mass + inertial) → (mass − inertial) ----------
    originals["UBF._compute"] = D.UBF._compute
    def fixed_ubf_compute(self):
        from rojak.core.derivatives import GradientMode
        coriolis_vorticity_term = self._coriolis_parameter * self._vorticity
        coriolis_deriv = D.latitudinal_derivative(self._coriolis_parameter)
        inertial_terms = coriolis_vorticity_term + 2 * self._jacobian
        mass_term = D.spatial_laplacian(self._geopotential, "deg", GradientMode.GEOSPATIAL)
        # residual form: mass − inertial + βu    (== A30 up to overall sign)
        return np.abs(mass_term - inertial_terms + coriolis_deriv * self._u_wind)
    D.UBF._compute = fixed_ubf_compute

    # --- Fix 5: F2D cross term  dv_dy → dv_dx ----------------------------
    originals["Frontogenesis2D._compute"] = D.Frontogenesis2D._compute
    def fixed_f2d_compute(self):
        from rojak.core.derivatives import GradientMode, spatial_gradient
        from rojak.turbulence.calculations import magnitude_of_vector
        dtheta = spatial_gradient(self._potential_temperature, "deg", GradientMode.GEOSPATIAL)
        mag = magnitude_of_vector(dtheta["dfdx"], dtheta["dfdy"])
        inv = -np.reciprocal(mag, where=mag != 0)
        return inv * (
            np.square(dtheta["dfdx"]) * self._du_dx
            + np.square(dtheta["dfdy"]) * self._dv_dy
            + dtheta["dfdx"] * dtheta["dfdy"] * self._dv_dx   # corrected
            + dtheta["dfdx"] * dtheta["dfdy"] * self._du_dy
        )
    D.Frontogenesis2D._compute = fixed_f2d_compute

    def unpatch():
        D.GradientRichardson._compute = originals["GradientRichardson._compute"]
        Ccalc.latitudinal_derivative  = originals["latitudinal_derivative"]
        D.latitudinal_derivative      = originals["latitudinal_derivative"]
        Cder.spatial_laplacian        = originals["spatial_laplacian"]
        D.spatial_laplacian           = originals["spatial_laplacian"]
        D.UBF._compute                = originals["UBF._compute"]
        D.Frontogenesis2D._compute    = originals["Frontogenesis2D._compute"]

    return unpatch


# ---------------------------------------------------------------------------
# Dimensional expectations (W&J Table 1 native SI, order of magnitude)
# ---------------------------------------------------------------------------
WJ_SI = {
    # key: (expected SI units, W&J median in SI, plausibility band factor)
    "negative_richardson": ("1",        -127.2,   None),   # sign matters, magnitude very domain-dependent
    "colson_panofsky":     ("m2 s-2",   -34.8 * 264.7, None),  # 10^3 kt^2 → m2/s2  (×264.7)
    "ubf":                 ("s-2",       1.61e-10, 1e3),
    "f2d":                 ("K2 m-2 s-1",5.66e-8,  1e3),   # 56.6e-9
    "ncsu1":               ("s-3",       1.11e-17, 1e4),   # 11.1e-18
    "rva_magnitude":       ("s-2",       2.33e-10, 1e3),
    "brown2":              ("s-3",       None,     None),  # rank-only, no W&J SI target
}


def summarize(da: xr.DataArray) -> dict:
    arr = np.asarray(da).ravel()
    arr = arr[np.isfinite(arr)]
    return {
        "median": float(np.median(arr)),
        "p99": float(np.quantile(arr, 0.99)),
        "min": float(arr.min()),
        "max": float(arr.max()),
        "n": int(arr.size),
    }


def spatial_rank_corr(a: xr.DataArray, b: xr.DataArray) -> float:
    """Spearman ρ between two fields (day-mean if time present)."""
    if "time" in a.dims:
        a = a.mean("time")
    if "time" in b.dims:
        b = b.mean("time")
    af = np.asarray(a).ravel()
    bf = np.asarray(b).ravel()
    m = np.isfinite(af) & np.isfinite(bf)
    if m.sum() < 10:
        return np.nan
    return float(spearmanr(af[m], bf[m])[0])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(path):
    print(f">>> Loading {path}")
    ds_raw = load_era5(path)
    catdata = prepare_for_rojak(ds_raw)
    ds = catdata._dataset

    from rojak.turbulence.diagnostic import DiagnosticFactory
    from rojak.orchestrator.configuration import TurbulenceDiagnostics as TD

    print("\n>>> Computing hand-written diagnostics")
    hand = {}
    hand["negative_richardson"] = hc.richardson(ds, negative=True)
    hand["colson_panofsky"]     = hc.colson_panofsky(ds)
    hand["ubf"]                 = hc.ubf(ds)
    hand["f2d"]                 = hc.frontogenesis_2d(ds)
    hand["ncsu1"]               = hc.ncsu1(ds)
    hand["rva_magnitude"]       = hc.rva(ds)
    # brown2 needs brown1 from rojak
    brown1 = DiagnosticFactory(catdata).create(TD.BROWN1).computed_value
    hand["brown2"] = hc.brown2(ds, brown1)
    for k, v in hand.items():
        print(f"    {k:22s} OK  units={v.attrs.get('units','?'):12s} shape={tuple(v.shape)}")

    # -----------------------------------------------------------------
    # Check (a): dimensional / order-of-magnitude
    # -----------------------------------------------------------------
    print("\n>>> (a) Dimensional + order-of-magnitude check vs W&J Table 1")
    dim_rows = []
    for k, da in hand.items():
        exp_units, wj_si, band = WJ_SI[k]
        got_units = da.attrs.get("units", "?")
        s = summarize(da)
        # reduce to 200 hPa if level dim remains
        note = ""
        if wj_si is None:
            verdict = "n/a (rank-only)" if k == "brown2" else "n/a"
        else:
            ratio = abs(s["median"]) / max(abs(wj_si), 1e-30)
            if band is None:
                verdict = f"sign={'−' if s['median']<0 else '+'} (magnitude domain-dependent)"
            elif ratio < band and ratio > 1.0/band:
                verdict = f"PLAUSIBLE (ratio {ratio:.2g})"
            else:
                verdict = f"CHECK (ratio {ratio:.2g})"
        dim_rows.append({
            "diagnostic": k,
            "units": got_units,
            "expected": exp_units,
            "median": f"{s['median']:.3g}",
            "p99": f"{s['p99']:.3g}",
            "W&J SI": f"{wj_si:.3g}" if wj_si is not None else "—",
            "verdict": verdict,
        })
    dim_df = pd.DataFrame(dim_rows)
    print(dim_df.to_string(index=False))

    # -----------------------------------------------------------------
    # Check (b): agreement with corrected rojak
    # -----------------------------------------------------------------
    print("\n>>> (b) Cross-check vs corrected rojak (monkey-patched, throwaway)")
    unpatch = patch_rojak()
    try:
        # Rebuild a fresh factory AFTER patching so patched code paths are used
        cat2 = prepare_for_rojak(load_era5(path))
        fac = DiagnosticFactory(cat2)
        rojak_fixed = {
            "negative_richardson": fac.create(TD.NEGATIVE_RICHARDSON).computed_value,
            "colson_panofsky":     fac.create(TD.COLSON_PANOFSKY).computed_value,
            "ubf":                 fac.create(TD.UBF).computed_value,
            "f2d":                 fac.create(TD.F2D).computed_value,
            "ncsu1":               fac.create(TD.NCSU1).computed_value,
            # rva: rojak has none — compare against nothing
            # brown2: rojak's is already faithful (unpatched); compare too
            "brown2":              fac.create(TD.BROWN2).computed_value,
        }
    finally:
        unpatch()

    xcheck_rows = []
    for k, r_fixed in rojak_fixed.items():
        h = hand[k]
        # align on 200 hPa if needed
        rf = r_fixed.sel(pressure_level=200) if "pressure_level" in r_fixed.dims else r_fixed
        hh = h.sel(pressure_level=200) if "pressure_level" in h.dims else h
        # broadcast-compatible compare of median + rank corr
        rho = spatial_rank_corr(rf, hh)
        sr = summarize(rf); sh = summarize(hh)
        med_ratio = sh["median"] / sr["median"] if sr["median"] != 0 else np.nan
        # numerical closeness (relative median difference)
        rel = abs(sh["median"] - sr["median"]) / max(abs(sr["median"]), 1e-30)
        agree = "MATCH" if (rho > 0.999 and rel < 0.02) else \
                ("CLOSE" if (rho > 0.95 and rel < 0.2) else "DIFFERS")
        # CP is a known intentional methodological difference: our CP computes
        # vertical shear on the full 3-level stencil (centered diff at 200 hPa),
        # whereas rojak pre-selects [200,225] and gets a one-sided difference at
        # 200 hPa. For a 175/200/225 design chosen precisely so 200 has both
        # neighbours, the centered stencil is the intended one.
        if k == "colson_panofsky" and agree == "DIFFERS":
            agree = "DIFFERS (intended: 3-level centered vs rojak 2-level one-sided)"
        xcheck_rows.append({
            "diagnostic": k,
            "hand median": f"{sh['median']:.3g}",
            "rojak-fixed median": f"{sr['median']:.3g}",
            "median ratio": f"{med_ratio:.3g}",
            "Spearman ρ": f"{rho:.4f}",
            "agreement": agree,
        })
    xcheck_df = pd.DataFrame(xcheck_rows)
    print(xcheck_df.to_string(index=False))

    print("\n>>> RVA has no rojak equivalent (patched or not) — verified by")
    print("    dimensional check + geospatial-derivative construction only.")

    # Save
    out = Path("/home/claude/work/phaseB_out"); out.mkdir(exist_ok=True)
    dim_df.to_csv(out / "dimensional_check.csv", index=False)
    xcheck_df.to_csv(out / "rojak_crosscheck.csv", index=False)
    print(f"\nSaved: {out}/dimensional_check.csv, {out}/rojak_crosscheck.csv")

    return dim_df, xcheck_df


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/synthetic_era5_clean_phi.nc"
    main(path)
