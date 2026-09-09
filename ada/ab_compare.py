"""
ada/ab_compare.py
=================
Stage 1 of the post-audit plan: run the 21 diagnostics TWICE on the SAME input
-- baseline versus a named set of audit fixes -- and answer one question per
diagnostic:

    does this fix change which cells exceed the threshold?

WHY THAT QUESTION, AND NOT "HOW MUCH DID THE VALUES CHANGE"
-----------------------------------------------------------
`4_verify.py` reports Spearman rho and median relative difference. Those are the
right numbers for "are these two implementations the same formula". They are the
wrong numbers for "does this change the science", because severity thresholds
are PERCENTILES OF THIS PROJECT'S OWN DATA. A constant scale error -- units,
resolution, an unknown length^2, a uniform bias -- cancels EXACTLY out of the
exceedance field (AUDIT_diagnostics_vs_literature.md 5.2, STATUS 17.5). A fix
can move every value by 8000x and change nothing that is ever published.

So the decision metric here is the SYMMETRIC DIFFERENCE OF THE EXCEEDANCE SET:

    flip = W(A xor B) / W(A or B)

with A the baseline exceedance mask, B the variant's, each against ITS OWN
cos(phi)-weighted percentile, and W a cos(phi)-weighted count. A fix with
rho = 0.9999 can still flip several per cent of the p97 tail; a fix with
flip = 0 cannot change any published number no matter what it did to the values.

Note what is deliberately NOT reported as the headline: the overall exceedance
FREQUENCY. Against a self-calibrated percentile the frequency is fixed at
(100 - p)% by construction, for both versions. What a fix can change is WHERE
those cells are -- so the frequency is reported PER LATITUDE BAND, which is the
quantity that actually differs between a globally-calibrated threshold and a
regionally-measured exceedance. F1 is predicted to show a monotone tilt across
those bands and almost nothing anywhere else; if it shows a uniform shift
instead, the audit's 3.1 derivation is wrong.

PRE-REGISTERED DECISION RULE (fix it before you see the numbers)
-----------------------------------------------------------------
  flip < 0.5 % at every severity   -> INERT. Merge it because the citation says
                                      so, without re-deriving 42 years.
  flip >= 0.5 % at any severity    -> MATERIAL. It goes into the single batched
                                      re-run, and its predicted direction (audit
                                      6) must be checked afterwards.

0.5 % is a judgement, not a theorem, and it is written down here so that it is
chosen before the numbers rather than after them.

USAGE
-----
    python ada/ab_compare.py <in.grib|in.nc> --fix meridional_metric
    python ada/ab_compare.py <in.grib> --fix ubf_computed_vorticity --fix pv_sharman_a18
    python ada/ab_compare.py <in.grib> --fix all --out data_prosser/ab/
    python ada/ab_compare.py <in.grib> --f2d-variant-baseline C \
                                       --f2d-variant-variant A

One input file. This is a SMALL-SAMPLE instrument -- one month is ~248 timesteps
and is plenty to measure a flip rate, but the thresholds it derives are not the
production thresholds and must not be used as such.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from calib_weighted_percentile import weighted_percentile  # noqa: E402


def _load(module_name: str, filename: str):
    """Load a numbered module by path. Registering in sys.modules BEFORE
    exec_module is mandatory -- see tests/conftest.py for why."""
    spec = importlib.util.spec_from_file_location(module_name, REPO / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


# Williams (2017) Table 1, p. 580 -- the percentile cut for each severity onset.
SEVERITIES: dict[str, float] = {
    "light": 97.0,
    "light_moderate": 99.1,
    "moderate": 99.6,
    "moderate_severe": 99.8,
    "severe": 99.9,
}

# AUDIT 4.0: derived from Williams (2017) Table 1 + Table 2 and corroborated by
# Sharman Table B1 p. 285 -- for ALL 21 the turbulent tail is the UPPER one,
# including the two with negative thresholds (-Ri and Colson-Panofsky). So the
# exceedance test is `>= threshold` for every diagnostic, with no per-diagnostic
# sign table. If that ever stops being true, it stops being true here first.
FLIP_INERT = 0.005          # the pre-registered 0.5 % rule
MIN_EXCEED_CELLS = 2000     # below this a flip rate is quantised noise, not a rate
LAT_BAND_WIDTH = 10.0       # degrees


def _flat(da: xr.DataArray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flatten to (values, cos-phi weights, latitudes), matching the weighting
    calibration.py applies (lat_weights = cos(phi), globally normalised)."""
    lat = da["latitude"]
    weights = np.cos(np.deg2rad(lat)).clip(min=0.0)
    w = weights.broadcast_like(da)
    latb = lat.broadcast_like(da)
    return (np.asarray(da.values, dtype=np.float64).ravel(),
            np.asarray(w.values, dtype=np.float64).ravel(),
            np.asarray(latb.values, dtype=np.float64).ravel())


def _spearman(a: np.ndarray, b: np.ndarray, cap: int = 400_000) -> float:
    """Spearman rho on the common finite support. Subsampled above `cap`
    because argsort on 10^8 elements is the whole runtime otherwise; the
    subsample is deterministic so two runs of this script agree."""
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float("nan")
    av, bv = a[m], b[m]
    if av.size > cap:
        idx = np.linspace(0, av.size - 1, cap).astype(np.int64)
        av, bv = av[idx], bv[idx]
    return float(spearmanr(av, bv).statistic)


def compare_one(base: xr.DataArray, var: xr.DataArray) -> dict:
    """All metrics for one diagnostic."""
    a, wa, lat = _flat(base)
    b, wb, _ = _flat(var)
    assert a.shape == b.shape, "baseline and variant must share a grid"

    finite = np.isfinite(a) & np.isfinite(b)
    rho = _spearman(a, b)

    denom = np.abs(a[finite])
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = np.abs(b[finite] - a[finite]) / np.where(denom > 0, denom, np.nan)
    median_rel = float(np.nanmedian(rel)) if rel.size else float("nan")
    identical = bool(np.array_equal(a[finite], b[finite]))

    out = {
        "n_finite": int(finite.sum()),
        "identical": identical,
        "spearman_rho": rho,
        "median_rel_diff": median_rel,
        "severities": {},
        "lat_bands": {},
    }

    edges = np.arange(-90.0, 90.0 + LAT_BAND_WIDTH, LAT_BAND_WIDTH)

    # ONE call with all five percentiles, not one call per severity.
    # `weighted_percentile` argsorts (and np.uniques) its whole input on every
    # call, so asking for the severities one at a time sorts the same 30-million
    # element array ten times instead of twice. STATUS.md 11.5 #2 records this
    # exact mistake taking out job 1092585 on walltime at 1:01:55; do not
    # reintroduce it by moving these back inside the loop.
    _pcts = np.asarray(list(SEVERITIES.values()), dtype=float)
    _thr_a = np.atleast_1d(weighted_percentile(a, wa, _pcts))
    _thr_b = np.atleast_1d(weighted_percentile(b, wb, _pcts))

    for _k, (name, pct) in enumerate(SEVERITIES.items()):
        thr_a = float(_thr_a[_k])
        thr_b = float(_thr_b[_k])
        # Upper tail for all 21 -- audit 4.0.
        A = finite & (a >= thr_a)
        B = finite & (b >= thr_b)
        w = wa
        w_or = float(w[A | B].sum())
        w_xor = float(w[A ^ B].sum())
        flip = (w_xor / w_or) if w_or > 0 else float("nan")
        n_exceed = int(A.sum())
        out["severities"][name] = {
            "percentile": pct,
            "n_exceed_baseline": n_exceed,
            "threshold_baseline": thr_a,
            "threshold_variant": thr_b,
            "threshold_rel_shift": (thr_b - thr_a) / abs(thr_a) if thr_a != 0 else float("nan"),
            "exceedance_freq_baseline": float(w[A].sum() / w[finite].sum()),
            "exceedance_freq_variant": float(w[B].sum() / w[finite].sum()),
            "flip_rate": flip,
            "material": bool(np.isfinite(flip) and flip >= FLIP_INERT),
        }
        if name == "light":
            bands = {}
            for lo, hi in zip(edges[:-1], edges[1:]):
                band = finite & (lat >= lo) & (lat < hi)
                wsum = float(wa[band].sum())
                if wsum <= 0:
                    continue
                fa = float(wa[band & A].sum() / wsum)
                fb = float(wa[band & B].sum() / wsum)
                bands[f"{lo:+.0f}..{hi:+.0f}"] = {
                    "freq_baseline_pct": 100 * fa,
                    "freq_variant_pct": 100 * fb,
                    "delta_pp": 100 * (fb - fa),
                }
            out["lat_bands"] = bands
    return out


def _cap_timesteps(ds, max_timesteps: int | None):
    """Take the first N timesteps, if asked.

    WHY THIS EXISTS. A flip rate needs enough cells in the tail to be a rate
    rather than a quantised count (see MIN_EXCEED_CELLS), and it needs the full
    LATITUDE range for the F1 tilt table -- but it does not need a long record.
    A global month at 0.25 deg is 721 x 1440 x 248, and this tool holds TWO
    complete result sets plus two persisted input datasets, so an uncapped
    global month is ~43 GB of diagnostic arrays before any intermediate. 32
    timesteps is 33 million cells per diagnostic, whose p99.9 set alone is
    ~33,000 cells -- an order of magnitude above the guard -- and it fits the
    48 GB request that schedules immediately on this cluster (STATUS 11.9:
    120 GB waited 22 hours where 24 GB started in 13 seconds).
    """
    if max_timesteps is None:
        return ds
    if max_timesteps < 3:
        raise SystemExit(
            f"--max-timesteps {max_timesteps} is unusable: f2d takes a material "
            f"derivative and xarray.differentiate needs at least 3 timesteps for "
            f"a centred difference (2 gives the same one-sided slope at both "
            f"points; 1 fails outright). Use 32 or more.")
    for name in ("time", "valid_time"):
        if name in ds.dims:
            n = ds.sizes[name]
            if n > max_timesteps:
                print(f"  capping {name}: {n} -> {max_timesteps} timesteps "
                      f"(--max-timesteps)")
                return ds.isel({name: slice(0, max_timesteps)})
            return ds
    return ds


def run(path: Path, fixes_on: list[str], target_level: int, out_dir: Path | None,
        f2d_variant: str | None, f2d_pair: tuple[str, str] | None = None,
        max_timesteps: int | None = None) -> dict:
    diag = _load("diagnostics", "2_diagnostics.py")

    baseline_fixes = diag.FixSet()
    variant_fixes = diag.FixSet(**{k: True for k in fixes_on})
    if variant_fixes == baseline_fixes and f2d_pair is None:
        raise SystemExit("ab_compare: no fixes requested -- nothing to compare. "
                         "Use --fix <name>, --fix all, or an "
                         "--f2d-variant-baseline/--f2d-variant-variant pair.")

    variant_label = variant_fixes.label()
    if f2d_pair is not None:
        variant_label += f"_f2d{f2d_pair[0]}to{f2d_pair[1]}"
    print(f"input        {path}")
    print(f"baseline     {baseline_fixes.label()}")
    print(f"variant      {variant_label}")
    print(f"target level {target_level} hPa")
    print(f"decision rule: flip >= {100 * FLIP_INERT:.1f} % at any severity = MATERIAL\n")

    ds_raw = _cap_timesteps(diag.load_era5(path), max_timesteps)
    kw_base = {"target_level": target_level}
    kw_var = {"target_level": target_level}
    if f2d_variant:
        kw_base["f2d_variant"] = kw_var["f2d_variant"] = f2d_variant
    if f2d_pair is not None:
        kw_base["f2d_variant"], kw_var["f2d_variant"] = f2d_pair
        print(f"f2d variant  baseline {f2d_pair[0]} -> variant {f2d_pair[1]}")

    print("computing baseline ...")
    base, fail_a = diag.compute_all_21(diag.prepare_for_rojak(ds_raw),
                                       fixes=baseline_fixes, **kw_base)
    print("computing variant ...")
    var, fail_b = diag.compute_all_21(diag.prepare_for_rojak(ds_raw),
                                      fixes=variant_fixes, **kw_var)

    if fail_a or fail_b:
        print(f"!! diagnostic failures -- baseline {[f.key for f in fail_a]}, "
              f"variant {[f.key for f in fail_b]}")

    keys = sorted(set(base) & set(var))
    results = {k: compare_one(base[k], var[k]) for k in keys}

    # ---------------- report ----------------
    print(f"\n{'diagnostic':<24}{'rho':>9}{'med|rel|':>11}"
          + "".join(f"{s[:9]:>10}" for s in SEVERITIES) + "   verdict")
    print("-" * (24 + 9 + 11 + 10 * len(SEVERITIES) + 11))
    material, inert, untouched = [], [], []
    for k in keys:
        r = results[k]
        flips = [r["severities"][s]["flip_rate"] for s in SEVERITIES]
        if r["identical"]:
            verdict = "UNTOUCHED"
            untouched.append(k)
        elif any(np.isfinite(f) and f >= FLIP_INERT for f in flips):
            verdict = "MATERIAL"
            material.append(k)
        else:
            verdict = "inert"
            inert.append(k)
        print(f"{k:<24}{r['spearman_rho']:>9.5f}{r['median_rel_diff']:>11.2e}"
              + "".join(f"{100 * f:>9.2f}%" for f in flips)
              + f"   {verdict}")

    # ---- small-sample guard ------------------------------------------------
    # A flip rate is a ratio of counts. On a small input the severe set can be a
    # handful of cells, and then "6.57 %" is one cell and means nothing. Say so
    # rather than letting the table imply a precision it does not have.
    n_light = min(results[k]["severities"]["light"]["n_exceed_baseline"] for k in keys)
    n_severe = min(results[k]["severities"]["severe"]["n_exceed_baseline"] for k in keys)
    if n_severe < MIN_EXCEED_CELLS:
        print(f"\n!! SMALL SAMPLE -- the severe (p99.9) exceedance set is only "
              f"{n_severe} cells (light: {n_light}). One cell is "
              f"{100 / max(n_severe, 1):.1f} % of it, so every flip rate in the "
              f"right-hand columns is quantised at that step and MATERIAL/inert "
              f"is not meaningful yet. Re-run on a full month "
              f"(~248 timesteps) before believing any verdict; use this run only "
              f"to confirm which diagnostics the fix touches at all.")

    print(f"\nMATERIAL  ({len(material):>2}): {', '.join(material) or '-'}")
    print(f"inert     ({len(inert):>2}): {', '.join(inert) or '-'}")
    print(f"untouched ({len(untouched):>2}): {', '.join(untouched) or '-'}")

    # Latitude tilt, at light severity, for the diagnostics that moved at all.
    moved = material + inert
    if moved:
        print(f"\nLATITUDE TILT at light (p97) -- change in exceedance frequency, "
              f"percentage points\nA monotone ramp here is the signature of a "
              f"latitude-structured error; a flat row is a scale change.")
        band_names = list(results[moved[0]]["lat_bands"])
        print(f"{'diagnostic':<24}" + "".join(f"{b:>12}" for b in band_names))
        for k in moved:
            row = results[k]["lat_bands"]
            print(f"{k:<24}" + "".join(f"{row[b]['delta_pp']:>+12.3f}" for b in band_names))

    payload = {
        "input": str(path),
        "target_level": target_level,
        "baseline": baseline_fixes.as_attrs(),
        "variant": variant_fixes.as_attrs(),
        "flip_inert_threshold": FLIP_INERT,
        "severity_percentiles": SEVERITIES,
        "material": material, "inert": inert, "untouched": untouched,
        "per_diagnostic": results,
    }

    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        stem = f"ab_{path.stem}_{variant_label}"
        if max_timesteps:
            stem += f"_n{max_timesteps}"
        (out_dir / f"{stem}.json").write_text(json.dumps(payload, indent=2))
        rows = ["diagnostic,identical,spearman_rho,median_rel_diff," +
                ",".join(f"flip_{s}" for s in SEVERITIES) + ",verdict"]
        for k in keys:
            r = results[k]
            v = ("UNTOUCHED" if k in untouched else
                 "MATERIAL" if k in material else "inert")
            rows.append(",".join([
                k, str(r["identical"]), f"{r['spearman_rho']:.6f}",
                f"{r['median_rel_diff']:.6e}",
                *[f"{r['severities'][s]['flip_rate']:.6f}" for s in SEVERITIES],
                v]))
        (out_dir / f"{stem}.csv").write_text("\n".join(rows) + "\n")
        print(f"\nwritten: {out_dir / (stem + '.json')}\n         "
              f"{out_dir / (stem + '.csv')}")

    return payload


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path)
    p.add_argument("--fix", action="append", default=[],
                   choices=["meridional_metric", "ubf_computed_vorticity",
                            "pv_sharman_a18", "endlich_component_shear",
                            "ncsu1_computed_vorticity", "all"],
                   help="repeatable; 'all' enables every audit fix")
    p.add_argument("--target-level", type=int, default=200)
    p.add_argument("--f2d-variant", default=None,
                   help="same variant on BOTH sides; default is "
                        "2_diagnostics.F2D_DEFAULT_VARIANT")
    p.add_argument("--f2d-variant-baseline", default=None,
                   help="f2d variant for the baseline side. Use with "
                        "--f2d-variant-variant to A/B two readings of Sharman "
                        "A9 (audit 6 F6), e.g. C against A.")
    p.add_argument("--f2d-variant-variant", default=None,
                   help="f2d variant for the variant side")
    p.add_argument("--max-timesteps", type=int, default=None,
                   help="use only the first N timesteps. REQUIRED in practice "
                        "for a global month: this tool holds two full result "
                        "sets, and 248 global timesteps is ~43 GB of arrays "
                        "before intermediates. 32 is ample for a flip rate.")
    p.add_argument("--out", type=Path, default=None,
                   help="directory for the JSON and CSV; omit to print only")
    a = p.parse_args()

    fixes = a.fix
    if "all" in fixes:
        fixes = ["meridional_metric", "ubf_computed_vorticity",
                 "pv_sharman_a18", "endlich_component_shear",
                 "ncsu1_computed_vorticity"]
    f2d_pair = (a.f2d_variant_baseline, a.f2d_variant_variant)
    if any(f2d_pair) and not all(f2d_pair):
        p.error("--f2d-variant-baseline and --f2d-variant-variant must be "
                "given together")
    if not fixes and not all(f2d_pair):
        p.error("nothing to compare -- pass at least one --fix, or an "
                "--f2d-variant-baseline/--f2d-variant-variant pair")

    run(a.input, sorted(set(fixes)), a.target_level, a.out, a.f2d_variant,
        f2d_pair if all(f2d_pair) else None, a.max_timesteps)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
