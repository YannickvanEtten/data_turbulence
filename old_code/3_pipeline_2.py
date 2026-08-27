"""
3_pipeline.py
=============
The pipeline: read an ERA5 GRIB/NetCDF file, compute all 21 diagnostics
(via 2_diagnostics.py), print summary statistics, build the Williams & Joshi
(2013) comparison table, and save one diagnostics.nc + diagnostics.zarr with
all 21 variables.

This file is deliberately THIN — all diagnostic logic lives in 2_diagnostics.py.
Here we only orchestrate: load -> compute_all_21 -> stats -> compare -> save.

Usage:
    python 3_pipeline.py <path/to/era5.grib> [output_dir]
"""
from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# import the diagnostics module (filename starts with a digit, so load by path)
_spec = importlib.util.spec_from_file_location(
    "diagnostics", Path(__file__).with_name("2_diagnostics.py"))
diag = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(diag)

REFERENCE_TABLE = diag.REFERENCE_TABLE

# --- Q-DISPATCH-2 (this task): wire in the 4 standalone Layer 2-6 modules ---
def _load_by_path(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

calib = _load_by_path("calib_weighted_percentile", "calib_weighted_percentile.py")   # Layer 2
agg = _load_by_path("aggregate", "aggregate.py")                                      # Layers 3-4
annual_agg = _load_by_path("annual_aggregate", "annual_aggregate.py")                 # Layer 5
trend_mod = _load_by_path("trend", "trend.py")                                        # Layer 6

@dataclass
class DiagStats:
    key: str
    units_native: str
    n_finite: int
    min: float
    max: float
    mean: float
    median: float
    p95: float
    p99: float
    p99_9: float | None


def summarise(name: str, da: xr.DataArray) -> DiagStats:
    arr = np.asarray(da).ravel()
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return DiagStats(name, str(da.attrs.get("units", "?")), 0, *([np.nan] * 7))
    p99_9 = float(np.quantile(arr, 0.999)) if arr.size >= 1000 else None
    return DiagStats(
        key=name,
        units_native=str(da.attrs.get("units", "?")),
        n_finite=int(arr.size),
        min=float(np.min(arr)),
        max=float(np.max(arr)),
        mean=float(np.mean(arr)),
        median=float(np.median(arr)),
        p95=float(np.quantile(arr, 0.95)),
        p99=float(np.quantile(arr, 0.99)),
        p99_9=p99_9,
    )


# Scale factors that convert from native (SI) rojak output into the
# units used in Williams & Joshi (2013) Table 1.  Multiply native value
# by this to get the table-units value.
#
# A small number of diagnostics need a pre-transform first:
#   - `deformation`: rojak returns DEF² (s⁻²), so we take √ first then × 1e6
#     to compare against W&J's DEF in 10⁻⁶ s⁻¹.
SCALE_TO_TABLE: dict[str, float] = {
    "magnitude_pv":          1.0e6,   # K m^2/(kg s) → PVU
    "colson_panofsky":       1.0,     # rojak returns m²/s²; W&J table is 10³ kt². Conversion factor is non-trivial — see note below
    "brown1":                1.0e6,   # s^-1 → 10^-6 s^-1
    "temperature_gradient":  1.0e6,   # K/m → 10^-6 K/m
    "horizontal_divergence": 1.0e6,   # s^-1 → 10^-6 s^-1
    "vertical_wind_shear":   1.0e3,   # s^-1 → 10^-3 s^-1
    "endlich":               1.0e3,   # rad/s → 10^-3 rad/s
    "deformation":           1.0e6,   # special: √(rojak) × 1e6  (see PRETRANSFORM below)
    "wind_speed":            1.0,
    "ngm2":                  1.0e9,
    "negative_richardson":   1.0,
    "rva_magnitude":         1.0e10,
    "ubf":                   1.0e12,
    "nva":                   1.0e10,
    "brown2":                1.0e6,
    "vorticity_squared":     1.0e9,
    "ti1":                   1.0e9,
    "ngm1":                  1.0e3,
    "ti2":                   1.0e9,
    "f2d":                   1.0e9,
    "ncsu1":                 1.0e18,
}

# Optional pre-transform applied to the native rojak value *before* scaling
# (used where rojak returns the squared form of a quantity but W&J Table 1
#  lists the un-squared form).
PRETRANSFORM: dict[str, callable] = {
    "deformation": lambda x: np.sqrt(np.abs(x)),
}


def build_comparison_table(stats: dict[str, DiagStats], diagnostics: dict[str, xr.DataArray]) -> pd.DataFrame:
    rows = []
    for key, ref in REFERENCE_TABLE.items():
        s = stats.get(key)
        if s is None:
            rows.append({
                "#": ref["num"],
                "Diagnostic": ref["name"],
                "W&J units": ref["units"],
                "W&J median": ref["wj_median"],
                "Trial median (table units)": "—",
                "Trial p99 (table units)": "—",
                "Status": "NOT COMPUTED",
            })
            continue
        scale = SCALE_TO_TABLE[key]
        # Apply any pre-transform first (e.g. √ for rojak's DEF²)
        if key in PRETRANSFORM and key in diagnostics:
            arr = np.asarray(diagnostics[key]).ravel()
            arr = arr[np.isfinite(arr)]
            arr = PRETRANSFORM[key](arr)
            med_native = float(np.median(arr))
            p99_native = float(np.quantile(arr, 0.99))
        else:
            med_native = s.median
            p99_native = s.p99
        med_t = med_native * scale
        p99_t = p99_native * scale
        wj = ref["wj_median"]
        ratio = abs(med_t) / max(abs(wj), 1e-30)
        if 0.1 <= ratio <= 10:
            status = "PLAUSIBLE"
        elif 0.01 <= ratio <= 100:
            status = "1-2 ORDERS OFF"
        else:
            status = "FLAG: magnitude mismatch"
        # Q-AGG-3 fix: this used to check `ref["sign"] == "+" and med_native < 0`,
        # which assumed sign="+" implies a positive median. That held for the
        # other 19 diagnostics but is false for colson_panofsky/negative_richardson
        # (#2, #11) -- both are correctly sign="+" (single-tailed, value>=threshold)
        # AND have negative W&J medians (Table 2's ladder is entirely negative-
        # valued). Comparing against wj's own sign generalizes correctly: it still
        # flags the 19 "should be positive" diagnostics if they go negative, and
        # now also correctly flags #2/#11 only if THEIR median flips positive
        # (genuinely wrong), not merely for being negative (expected).
        if wj != 0 and np.sign(med_native) != np.sign(wj):
            status += " | SIGN MISMATCH vs W&J"
        rows.append({
            "#": ref["num"],
            "Diagnostic": ref["name"],
            "W&J units": ref["units"],
            "W&J median": wj,
            "Trial median (table units)": f"{med_t:.3g}",
            "Trial p99 (table units)": f"{p99_t:.3g}",
            "Status": status,
        })
    return pd.DataFrame(rows).set_index("#")


# ---------------------------------------------------------------------------
# 5. Save outputs
# ---------------------------------------------------------------------------

def save_outputs(diagnostics: dict[str, xr.DataArray], out_dir: Path) -> tuple[Path, Path]:
    """Save all diagnostics to one NetCDF file and one Zarr store.

    Zarr writing needs a clean target directory. On Windows, a leftover
    store from a previous run (or a file handle OneDrive/Explorer still
    has open) can make shutil.rmtree raise PermissionError. We retry a
    few times with a short delay, and if it still fails, write to a
    fresh timestamped zarr path instead of crashing the whole run --
    the NetCDF file (the primary output) is saved either way.
    """
    import shutil, time

    out_dir.mkdir(parents=True, exist_ok=True)
    ds_out = xr.Dataset(diagnostics)
    nc_path = out_dir / "diagnostics.nc"
    zarr_path = out_dir / "diagnostics.zarr"

    # NetCDF first -- this is the primary output and rarely has locking issues
    ds_out.to_netcdf(nc_path)

    # Zarr: remove any existing store, tolerating Windows file-lock delays
    if zarr_path.exists():
        for attempt in range(5):
            try:
                shutil.rmtree(zarr_path)
                break
            except PermissionError:
                if attempt == 4:
                    # give up cleaning the old store; write to a fresh path instead
                    import datetime
                    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                    zarr_path = out_dir / f"diagnostics_{stamp}.zarr"
                    print(f"    (could not remove old zarr store -- Windows file lock; "
                          f"writing to {zarr_path.name} instead. You can delete the old "
                          f"'diagnostics.zarr' folder manually once no program has it open.)")
                    break
                time.sleep(1.0)

    ds_out.to_zarr(zarr_path, mode="w", zarr_format=2)
    return nc_path, zarr_path


# ---------------------------------------------------------------------------
# 6. Top-level entry point
# ---------------------------------------------------------------------------



# ===========================================================================
# Top-level entry point
# ===========================================================================
def run(input_path: str | Path, output_dir: str | Path = "./cat_outputs") -> dict:
    """End-to-end: load -> compute 21 -> save -> stats -> comparison table."""
    input_path = Path(input_path); output_dir = Path(output_dir)
    print(f"\n>>> Loading {input_path}")
    ds_raw = diag.load_era5(input_path)
    print(f"    raw vars   : {list(ds_raw.data_vars)}")

    catdata = diag.prepare_for_rojak(ds_raw)
    print(f"    rojak dims : {dict(catdata._dataset.sizes)}")

    print("\n>>> Computing all 21 diagnostics (14 rojak + 7 hand-written)")
    diagnostics = diag.compute_all_21(catdata, target_level=200)
    for k in diagnostics:
        tag = "hand" if k in diag.HANDCODED_KEYS else "rojak"
        print(f"    [{tag:5s}] {k:22s} shape={tuple(diagnostics[k].shape)}")

    print("\n>>> Saving outputs")
    nc_path, zarr_path = save_outputs(diagnostics, output_dir)
    print(f"    NetCDF: {nc_path}")
    print(f"    Zarr  : {zarr_path}")

    print("\n>>> Summary statistics")
    stats = {k: summarise(k, v) for k, v in diagnostics.items()}
    stats_df = pd.DataFrame([s.__dict__ for s in stats.values()])
    stats_df.to_csv(output_dir / "summary_stats.csv", index=False)

    print(">>> Williams & Joshi (2013) comparison table")
    comp_df = build_comparison_table(stats, diagnostics)
    comp_df.to_csv(output_dir / "comparison_table.csv")

    return {"diagnostics": diagnostics, "stats_df": stats_df,
            "comparison_df": comp_df, "nc_path": nc_path,
            "zarr_path": zarr_path, "output_dir": output_dir}


# ===========================================================================
# Layers 2-6: calibration -> exceedance -> annual aggregation -> trend.
# This wiring didn't exist before this task. Two seams made EXPLICIT below,
# not implicit, per the ask:
#   (a) box / pressure_levels / calibration_domain are PipelineConfig fields,
#       not hardcoded constants (Q-FIDELITY-1 interim decision).
#   (b) the annual_aggregate 0-100 (%) -> trend 0-1 (probability) conversion
#       is its own visibly-logged step, not folded silently into either
#       neighboring layer.
# ===========================================================================
@dataclass
class PipelineConfig:
    """Q-FIDELITY-1: box/levels/calibration-domain as config, not constants.

    box:                the TREND/EXCEEDANCE/COMPARISON box (Layers 3-6).
                         Q-FIDELITY-1 point 1: Prosser's exact box for a real
                         comparison is 36-60N/55-10W -- NOT necessarily the
                         same as the download superset (which can stay wider
                         "for computational convenience" per STATUS_24, with
                         subsetting happening at comparison time only).
    pressure_levels:     Q-FIDELITY-1 point 2 -- settled at [175, 200, 225]
                         for the full 42-year pull, kept as a config field
                         (not a literal in this module) so calibration and
                         application are forced to visibly share one value,
                         per the "same stencil" self-consistency argument.
    calibration_domain:  Q-FIDELITY-1 point 3 -- "global" (Prosser's real
                         year-2000 domain, pending Q-GLOBAL-1's feasibility
                         sizing) or "regional-plumbing-test" (current NA-box
                         run, explicitly flagged non-scientific -- Opus's
                         near-tautological finding in Q-CALIB-2). This field
                         exists so a plumbing-test run can NEVER silently be
                         mistaken for a real calibration -- see run below.
    severities:          named severity levels and their percentile cutoffs,
                         e.g. {"light": 97.0, "light_moderate": 99.1, ...}.
    """
    box: dict[str, float]
    pressure_levels: list[int]
    calibration_domain: str
    severities: dict[str, float]
    start_year: int
    end_year: int


def run_layers_2_to_6(
    diagnostic_fields: dict[str, xr.DataArray],
    calibration_fields: dict[str, xr.DataArray],
    calibration_lat_weights: xr.DataArray,
    config: PipelineConfig,
) -> dict:
    """Layer 2 (weighted-percentile calibration) through Layer 6 (trend).

    diagnostic_fields:      {diagnostic_name: DataArray} -- the 21 raw
                             diagnostic values over config.box, full
                             multi-year time range (3-hourly). This is
                             Layer 1's output (compute_all_21, run over
                             the years config.start_year..end_year).
    calibration_fields:     {diagnostic_name: DataArray} -- the SAME 21
                             diagnostics computed over the SEPARATE
                             calibration domain (config.calibration_domain;
                             Q-GLOBAL-1's global year-2000 pull, or the
                             regional plumbing-test data). Deliberately a
                             DIFFERENT dataset from diagnostic_fields --
                             this is seam (a) made concrete: calibration
                             domain and trend/exceedance box are NOT
                             required to be (and per Q-FIDELITY-1, are
                             NOT) the same spatial domain.
    calibration_lat_weights: cos(phi) weights matching calibration_fields'
                             grid, for Layer 2's weighted_percentile.

    Returns per-severity: thresholds, exceedance_mean (per timestep),
    annual_probability_pct (0-100), annual_probability (0-1), trend.
    """
    if config.calibration_domain not in ("global", "regional-plumbing-test"):
        raise ValueError(
            f"calibration_domain must be 'global' or 'regional-plumbing-test', "
            f"got {config.calibration_domain!r} -- Q-FIDELITY-1 point 3 requires "
            f"this to be an explicit, deliberate choice, not a default."
        )
    if config.calibration_domain == "regional-plumbing-test":
        print(">>> WARNING: calibration_domain='regional-plumbing-test' -- this run "
              "is NOT a real calibration (Q-CALIB-2, near-tautological on the same "
              "box it's tested against). Results here validate PLUMBING ONLY.")

    # ---- Layer 2: weighted-percentile calibration -------------------------
    print(f"\n>>> Layer 2: weighted-percentile calibration (domain={config.calibration_domain})")
    thresholds: dict[str, dict[str, float]] = {}
    signs: dict[str, str] = {}
    for name in diagnostic_fields:
        ref = REFERENCE_TABLE[name]
        # Q-AGG-3: sign is read from REFERENCE_TABLE, not hardcoded here.
        # colson_panofsky/negative_richardson are correctly "+" as of the
        # Q-AGG-3 fix -- sign="either" capability still exists in aggregate.py
        # for any diagnostic that genuinely needs it, but nothing here
        # special-cases those two into it.
        signs[name] = ref["sign"]
        values = calibration_fields[name]
        thresholds[name] = {
            sev_name: calib.weighted_percentile(
                np.asarray(values).ravel(),
                np.asarray(calibration_lat_weights.broadcast_like(values)).ravel(),
                pct,
            )
            for sev_name, pct in config.severities.items()
        }
        print(f"    [{ref['sign']}] {name:22s} thresholds: "
              f"{ {k: f'{v:.4g}' for k, v in thresholds[name].items()} }")

    # ---- Layers 3-4: exceedance-first, then average across 21 diagnostics --
    print("\n>>> Layers 3-4: per-timestep exceedance-mean (exceedance-first ordering, Q-AGG-1)")
    severity_names = list(config.severities.keys())
    exceedance_mean = agg.exceedance_mean_all_severities(
        diagnostic_fields, thresholds, severity_names, signs
    )  # Q-CALIB-6 guard runs automatically inside this call

    # ---- Layer 5: annual aggregation, leap-aware normalization (Q-AGG-5) ---
    print("\n>>> Layer 5: annual aggregation (leap-aware 3-hour-period normalization)")
    annual_pct: dict[str, xr.DataArray] = {}
    for sev_name, field in exceedance_mean.items():
        annual_pct[sev_name] = annual_agg.annual_exceedance_probability(
            field, time_dim="time", year_dim="year"
        )
        print(f"    {sev_name}: annual % field, years "
              f"{int(annual_pct[sev_name]['year'].min())}-{int(annual_pct[sev_name]['year'].max())}")

    # ---- Seam (b): EXPLICIT 0-100 (%) -> 0-1 (probability) conversion ------
    # annual_aggregate.py returns PERCENTAGE (0-100, matching Prosser's
    # "percentage probabilities of exceedance" wording literally).
    # trend.py's fit_annual_trend() was built and synthetic-tested (Q-AGG-2)
    # on a 0-1 probability scale. This conversion is deliberately its own
    # visible step -- not folded silently into either annual_aggregate.py
    # or trend.py -- per the Q-AGG-5 flag.
    print("\n>>> Converting annual % (0-100) -> probability (0-1) before trend fitting")
    annual_prob: dict[str, xr.DataArray] = {
        sev_name: field / 100.0 for sev_name, field in annual_pct.items()
    }

    # ---- Layer 6: linear trend regression, negative-clipped ---------------
    print("\n>>> Layer 6: per-gridpoint linear trend regression (Q-AGG-2)")
    trend_results: dict[str, xr.Dataset] = {}
    for sev_name, prob_field in annual_prob.items():
        trend_results[sev_name] = trend_mod.fit_annual_trend(
            prob_field, year_dim="year",
            start_year=config.start_year, end_year=config.end_year,
        )

    return {
        "thresholds": thresholds,
        "exceedance_mean": exceedance_mean,
        "annual_probability_pct": annual_pct,
        "annual_probability": annual_prob,
        "trend": trend_results,
    }


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python 3_pipeline.py <path/to/era5.grib> [output_dir]")
        sys.exit(1)
    out = sys.argv[2] if len(sys.argv) >= 3 else "./cat_outputs"
    result = run(sys.argv[1], out)
    print("\n=== Comparison vs Williams & Joshi (2013) Table 1 ===")
    pd.set_option("display.max_colwidth", None); pd.set_option("display.width", 200)
    print(result["comparison_df"].to_string())