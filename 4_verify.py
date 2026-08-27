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

  (c) NEGATIVE-RI SUBSET CHECK (Q-DATA-1) — on a winter day, confirm the
      negative-Ri branches of #2/#11/#21 (never exercised on the summer
      trial day) actually agree with the corrected-rojak oracle where they
      fire. Located automatically wherever the negative-Ri cells actually
      are in the pressure-level stack -- do NOT assume 200 hPa.

Self-tests on a NetCDF/GRIB file given on the command line. The default is
era5_validation_subset.nc sitting beside this script -- real data, so the
default run is the meaningful one.

Results are written to verification/<YYYY-MM-DD>/ beside this script: the two
CSV tables PLUS a full console log, because checks (b) and (c) print their
Spearman rho values to stdout and nowhere else. Each run gets its own dated
folder, so a new run never overwrites the evidence from an old one.

Usage:
    python 4_verify.py                          # era5_validation_subset.nc
    python 4_verify.py path/to/data.nc          # a different input
    python 4_verify.py path/to/data.nc out_dir  # a different output folder
"""
from __future__ import annotations

import datetime as _dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import spearmanr


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "era5_validation_subset.nc"


def default_out_dir():
    """verification/<YYYY-MM-DD>/ beside this script. Dated so that re-running
    never overwrites an earlier run's evidence."""
    return HERE / "verification" / _dt.date.today().isoformat()


class _Tee:
    """Mirror stdout into a log file. Checks (b) and (c) print their Spearman
    rho values and nothing else captures them -- without this the strongest
    evidence in the project stays trapped in a terminal scrollback."""

    def __init__(self, stream, handle):
        self._stream = stream
        self._handle = handle

    def write(self, text):
        self._stream.write(text)
        self._handle.write(text)
        return len(text)

    def flush(self):
        self._stream.flush()
        self._handle.flush()

import importlib.util
from pathlib import Path
_spec = importlib.util.spec_from_file_location(
    "diagnostics", Path(__file__).with_name("2_diagnostics.py"))
hc = importlib.util.module_from_spec(_spec)
# MUST register in sys.modules BEFORE exec_module. 2_diagnostics.py uses
# `from __future__ import annotations`, so @dataclass resolves its field
# annotations through sys.modules[cls.__module__] -- which is None if we
# never registered it, and the import dies with a bare
# "'NoneType' object has no attribute '__dict__'" on Python 3.12.
# Same fix is needed in 3_pipeline.py (both its loaders).
sys.modules[_spec.name] = hc
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

    # --- Fix 3: spatial_laplacian -- NO LONGER NEEDED, fixed upstream -----
    # At the old pin (1a65326) rojak's spatial_laplacian returned
    # (df/dx + df/dy), which is not a Laplacian at all, so this harness used
    # to monkey-patch it.
    #
    # Upstream commit 25b8685 ("FIX: Geospatial laplacian was missing an extra
    # derivative", #234) repaired it, and repaired it CORRECTLY: it now routes
    # through divergence(), which uses vector_derivatives() and therefore
    # carries the spherical curvature term. That is the same correction this
    # project derived independently from the divergence theorem on 2026-08-26
    # (STATUS.md 4g) -- two independent routes to the same answer.
    #
    # The patch is deliberately NOT reinstated as a no-op equivalent. Every
    # monkey-patch here puts more of OUR code inside the "independent" oracle
    # and makes the cross-check weaker. Letting rojak use its own, now-correct
    # implementation makes the ubf comparison a genuinely independent check of
    # the Laplacian rather than a comparison of two copies of our own formula.
    #
    # If the pin is ever moved BACK before 25b8685, this patch must return.

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

    # --- Fix 6: ColsonPanofsky reduced form (Q-CP-3) ----------------------
    # Root cause (Q-CP-1): ColsonPanofsky.__init__ pre-slices u/v/geopotential
    # to 2 levels before recomputing its OWN Sv^2, while the Ri it receives
    # was built from GradientRichardson/VWS on the FULL unsliced 3-level
    # data. These two Sv^2's are not the same number (verified up to ~3000x
    # apart at low-shear cells), so Sv^2*(1-Ri/Ri_crit) never algebraically
    # reduces to Sv^2-N^2/Ri_crit -- it blows up. Fix: bypass the class's
    # own re-sliced Sv^2 entirely; use the SAME full-stencil VWS/BruntVaisala
    # diagnostics that feed Ri, then combine directly with no Ri round-trip.
    def patched_colson_panofsky(catdata, ri_crit=0.5):
        import rojak.turbulence.diagnostic as _D
        from rojak.orchestrator.configuration import TurbulenceDiagnostics as _TD
        fac = _D.DiagnosticFactory(catdata)
        vws = fac.create(_TD.VWS).computed_value            # Sv, full 3-level stencil
        n2 = fac.create(_TD.BRUNT_VAISALA).computed_value    # N^2, full 3-level stencil
        alt = catdata._dataset["altitude"]
        length_scale = alt.diff("pressure_level", label="upper")
        lvl = length_scale["pressure_level"]
        sv2 = (vws ** 2).sel(pressure_level=lvl)
        n2_sel = n2.sel(pressure_level=lvl)
        return np.square(length_scale) * (sv2 - n2_sel / ri_crit)

    def unpatch():
        D.GradientRichardson._compute = originals["GradientRichardson._compute"]
        Ccalc.latitudinal_derivative  = originals["latitudinal_derivative"]
        D.latitudinal_derivative      = originals["latitudinal_derivative"]
        D.UBF._compute                = originals["UBF._compute"]
        D.Frontogenesis2D._compute    = originals["Frontogenesis2D._compute"]

    return unpatch, patched_colson_panofsky


# ---------------------------------------------------------------------------
# Dimensional expectations (W&J Table 1 native SI, order of magnitude)
# ---------------------------------------------------------------------------
WJ_SI = {
    # key: (expected SI units, W&J median in SI, plausibility band factor)
    "negative_richardson": ("1",        -127.2,   None),   # sign matters, magnitude very domain-dependent
    "colson_panofsky":     ("m2 s-2",   -34.8 * 264.7, None),  # 10^3 kt^2 → m2/s2  (×264.7)
    "ubf":                 ("s-2",       1.61e-10, 1e3),
    "f2d":                 ("m2 s-3 K-2",5.66e-8,  1e3),   # 56.6e-9, Q-F2D-5: A9, was K2 m-2 s-1 Miller form
    "ncsu1":               ("s-3",       1.11e-17, 1e4),   # 11.1e-18
    "rva_magnitude":       ("s-2",       2.33e-10, 1e3),
    "brown2":              ("s-3",       None,     None),  # rank-only, no W&J SI target
}


def summarize(da: xr.DataArray) -> dict:
    arr = np.asarray(da).ravel()
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        # Report "no data" instead of crashing in np.quantile. A diagnostic
        # can legitimately be empty on a short input -- f2d drops the first
        # and last timestep for its centred d/dt, so a 2-timestep file leaves
        # it with nothing. That is a fact about the INPUT, and the harness
        # should say so rather than die halfway through the run.
        return {"median": np.nan, "p99": np.nan, "min": np.nan,
                "max": np.nan, "n": 0}
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


def rho_and_rel(h: xr.DataArray, r: xr.DataArray, mask: xr.DataArray):
    """Spearman rho and median relative difference, restricted to `mask`."""
    h_sub, r_sub = h.where(mask), r.where(mask)
    rho = spatial_rank_corr(r_sub, h_sub)
    hf = np.asarray(h_sub).ravel()
    rf = np.asarray(r_sub).ravel()
    m = np.isfinite(hf) & np.isfinite(rf)
    if m.sum() == 0:
        return rho, float("nan")
    rel = float(np.nanmedian(np.abs(hf[m] - rf[m]) / np.maximum(np.abs(rf[m]), 1e-30)))
    return rho, rel


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(path, out_dir=None):
    out = Path(out_dir) if out_dir is not None else default_out_dir()
    out.mkdir(parents=True, exist_ok=True)
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
    # Q-F2D-5: interior timesteps only (drops one-sided 00Z/21Z edges) -> true centered d/dt
    hand["f2d"]                 = hc.frontogenesis_isentropic(ds).isel(time=slice(1, -1))
    hand["ncsu1"]               = hc.ncsu1(ds)
    hand["rva_magnitude"]       = hc.rva(ds)
    # brown2 needs brown1 from rojak
    brown1 = DiagnosticFactory(catdata).create(TD.BROWN1).computed_value
    hand["brown2"] = hc.brown2(ds, brown1)
    for k, v in hand.items():
        print(f"    {k:22s} OK  units={v.attrs.get('units','?'):12s} shape={tuple(v.shape)}")

    neg_ri_all = hand["negative_richardson"]
    frac_negative_Ri = float((neg_ri_all > 0).mean())
    max_val = float(neg_ri_all.max())
    print(f"Fraction of cells with Ri<0 (pooled across levels): {frac_negative_Ri:.4%}")
    print(f"Max of -Ri (positive => genuine Ri<0 present): {max_val:.4g}")

    # -----------------------------------------------------------------
    # Check (a): dimensional / order-of-magnitude
    # -----------------------------------------------------------------
    print("\n>>> (a) Dimensional + order-of-magnitude check vs W&J Table 1")
    dim_rows = []
    for k, da in hand.items():
        exp_units, wj_si, band = WJ_SI[k]
        got_units = da.attrs.get("units", "?")
        s = summarize(da)
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
    unpatch, patched_colson_panofsky = patch_rojak()
    try:
        # Rebuild a fresh factory AFTER patching so patched code paths are used
        cat2 = prepare_for_rojak(load_era5(path))
        fac = DiagnosticFactory(cat2)
        rojak_fixed = {
            "negative_richardson": fac.create(TD.NEGATIVE_RICHARDSON).computed_value,
            # Q-CP-3: reduced form, bypasses ColsonPanofsky class entirely
            # (Fix 6 -- see patch_rojak docstring for root cause)
            "colson_panofsky":     patched_colson_panofsky(cat2),
            "ubf":                 fac.create(TD.UBF).computed_value,
            # f2d: Q-F2D-1 confirmed neither Frontogenesis2D (Miller) nor
            # Frontogenesis3D (physical-space) is isentropic A9 -- no rojak
            # equivalent to cross-check against. Dimensional check only,
            # same treatment as RVA (Q-RVA-1).
            "ncsu1":               fac.create(TD.NCSU1).computed_value,
            # rva: rojak has none — compare against nothing
            # brown2: rojak's is already faithful (unpatched); compare too
            "brown2":              fac.create(TD.BROWN2).computed_value,
        }

        xcheck_rows = []
        for k, r_fixed in rojak_fixed.items():
            h = hand[k]
            rf = r_fixed.sel(pressure_level=200) if "pressure_level" in r_fixed.dims else r_fixed
            hh = h.sel(pressure_level=200) if "pressure_level" in h.dims else h
            rho = spatial_rank_corr(rf, hh)
            sr = summarize(rf); sh = summarize(hh)
            med_ratio = sh["median"] / sr["median"] if sr["median"] != 0 else np.nan
            rel = abs(sh["median"] - sr["median"]) / max(abs(sr["median"]), 1e-30)
            agree = "MATCH" if (rho > 0.999 and rel < 0.02) else \
                    ("CLOSE" if (rho > 0.95 and rel < 0.2) else "DIFFERS")
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

        # ---------------------------------------------------------------
        # Check (c): Q-DATA-1 -- negative-Ri subset, wherever it actually is
        # ---------------------------------------------------------------
        print("\n>>> (c) Negative-Ri-only subset check (Q-DATA-1)")

        # negative_richardson keeps its native pressure_level coord (175/200/225)
        # -- find which level(s) actually have negative-Ri cells, don't assume 200.
        neg_ri_by_level = {}
        for lvl in np.atleast_1d(hand["negative_richardson"]["pressure_level"].values):
            n = int((hand["negative_richardson"].sel(pressure_level=lvl) > 0).sum())
            neg_ri_by_level[float(lvl)] = n
            print(f"    pressure_level={lvl:>6}: negative-Ri cells = {n}")

        target_lvl = max(neg_ri_by_level, key=neg_ri_by_level.get)
        n_neg = neg_ri_by_level[target_lvl]
        print(f"    -> using pressure_level={target_lvl} (most negative-Ri cells: {n_neg})")

        if n_neg == 0:
            print("    No negative-Ri cells at any level -- Q-DATA-1 not exercised on this "
                  "day/domain. Try a different date, don't force a check on empty data.")
        else:
            neg_mask = hand["negative_richardson"].sel(pressure_level=target_lvl) > 0

            # -- negative_richardson: native level, straightforward .sel()
            rho, rel = rho_and_rel(
                hand["negative_richardson"].sel(pressure_level=target_lvl),
                rojak_fixed["negative_richardson"].sel(pressure_level=target_lvl),
                neg_mask,
            )
            print(f"    negative_richardson    n={n_neg:5d}  rho={rho:.4f}  median_rel_diff={rel:.3g}")

            # -- colson_panofsky: pressure_level is DIFF-DERIVED (from alt.diff()),
            # its labels are not guaranteed to be 175/200/225. Print them, and use
            # method="nearest" rather than assuming a literal label match.
            cp_hand_levels = hand["colson_panofsky"]["pressure_level"].values
            cp_rojak_levels = rojak_fixed["colson_panofsky"]["pressure_level"].values
            print(f"    CP hand pressure_level coord:  {cp_hand_levels}")
            print(f"    CP rojak pressure_level coord: {cp_rojak_levels}")
            cp_h = hand["colson_panofsky"].sel(pressure_level=target_lvl, method="nearest")
            cp_r = rojak_fixed["colson_panofsky"].sel(pressure_level=target_lvl, method="nearest")
            print(f"    CP nearest-match levels used: hand={float(cp_h['pressure_level']):.1f}  "
                  f"rojak={float(cp_r['pressure_level']):.1f}  (target was {target_lvl})")
            rho, rel = rho_and_rel(cp_h, cp_r, neg_mask)
            print(f"    colson_panofsky        n={n_neg:5d}  rho={rho:.4f}  median_rel_diff={rel:.3g}")

            # -- ncsu1: hc.ncsu1(ds) already collapses to a single target_level
            # internally (default 200), so hand["ncsu1"] has NO pressure_level
            # dim to re-select. Must call it FRESH at target_lvl.
            ncsu1_hand_target = hc.ncsu1(ds, target_level=int(target_lvl))
            ncsu1_rojak_raw = rojak_fixed["ncsu1"]
            print(f"    NCSU1 rojak dims: {ncsu1_rojak_raw.dims}")
            ncsu1_rojak_target = (
                ncsu1_rojak_raw.sel(pressure_level=target_lvl, method="nearest")
                if "pressure_level" in ncsu1_rojak_raw.dims else ncsu1_rojak_raw
            )
            rho, rel = rho_and_rel(ncsu1_hand_target, ncsu1_rojak_target, neg_mask)
            print(f"    ncsu1                  n={n_neg:5d}  rho={rho:.4f}  median_rel_diff={rel:.3g}")

    finally:
        unpatch()

    print("\n>>> RVA and F2D (Q-F2D-5, isentropic A9) have no rojak equivalent")
    print("    (patched or not) — verified by dimensional check only.")

    # Save
    dim_df.to_csv(out / "dimensional_check.csv", index=False)
    xcheck_df.to_csv(out / "rojak_crosscheck.csv", index=False)
    print(f"\nSaved: {out}/dimensional_check.csv, {out}/rojak_crosscheck.csv")

    return dim_df, xcheck_df


def _run(path, out_dir):
    """Create the output folder, then run main() with stdout mirrored into
    verify_log.txt so checks (b) and (c) leave a record on disk."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "verify_log.txt"

    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write(f"# 4_verify.py\n# input : {path}\n")
        handle.write(f"# run at: {_dt.datetime.now().isoformat(timespec='seconds')}\n\n")
        original_stdout = sys.stdout
        sys.stdout = _Tee(original_stdout, handle)
        try:
            result = main(path, out_dir)
        finally:
            sys.stdout = original_stdout

    print(f"Console log saved: {log_path}")
    return result


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT
    out_dir = sys.argv[2] if len(sys.argv) > 2 else default_out_dir()
    if not Path(path).exists():
        raise SystemExit(
            f"input file not found: {path}\n"
            f"Pass one explicitly: python 4_verify.py path/to/data.nc"
        )
    _run(path, out_dir)