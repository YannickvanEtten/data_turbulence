"""
Check Stage 3 (thresholds) and Stage 4 (exceedance frequencies, trends) against
the ORIGINAL pipeline on small synthetic stores. Needs the original repo (the
parent folder of essentials/); no rojak, no ERA5.

    pixi run python essentials/tests/compare_stages_3_4_with_original.py [workdir]

Stage 3: 12 synthetic "global" stores -> original ada/tail_thresholds.py and
         essentials thresholds.calibrate -> the two JSON files must agree exactly.
         Also: tail method == whole-array weighted_percentile (to 1e-6, the original tolerance;
         the residue is np.interp amplifying last-bit differences in cumulative sums).
Stage 4: 3 synthetic years of "North Atlantic" stores -> per-year frequencies
         from the original aggregate.py + full_trend_check.weighted_rate, and from
         essentials trends.year_frequencies -> must agree to 1e-12.
"""
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve()
ESSENTIALS, ORIGINAL = HERE.parents[1], HERE.parents[2]
sys.path.insert(0, str(ESSENTIALS))
sys.path.insert(0, str(ORIGINAL))

import config as C            # noqa: E402
import thresholds as T        # noqa: E402
import trends as TR           # noqa: E402

rng = np.random.default_rng(42)
ORIGINAL_AUDIT_ATTRS = {"audit_fixes": "meridional_metric+four_variable_provenance",
                        "f2d_variant": "A", "deformation_convention": "DEF (un-squared, Sharman A17)"}


def synthetic_field(shape, key, lat):
    """Heavy-tailed, latitude-dependent, with NaNs and (for some) ties at zero."""
    x = rng.lognormal(0, 1 + 0.3 * np.cos(np.deg2rad(lat))[:, None, None], size=shape)
    if key in ("nva", "ncsu1"):
        x = np.where(rng.random(shape) < 0.5, 0.0, x)          # clipped diagnostics
    if key in ("negative_richardson", "colson_panofsky", "f2d"):
        x = x - 3.0                                            # negative thresholds
    x[rng.random(shape) < 0.001] = np.nan
    return x.astype(np.float32)


def write_store(path, lat, lon, times, attrs):
    shape = (lat.size, lon.size, times.size)
    ds = xr.Dataset({k: (("latitude", "longitude", "time"), synthetic_field(shape, k, lat))
                     for k in C.DIAGNOSTIC_KEYS},
                    coords={"latitude": lat, "longitude": lon, "time": times})
    ds.attrs.update(attrs)
    ds.to_zarr(path, mode="w", zarr_format=2, consolidated=True)


def month_times(year, month):
    start = np.datetime64(f"{year}-{month:02d}-01")
    end = np.datetime64(f"{year + (month == 12)}-{month % 12 + 1:02d}-01")
    return np.arange(start, end, np.timedelta64(3, "h")).astype("datetime64[ns]")


def stage3(work: Path) -> bool:
    base = work / "base3"
    sub = base / "derived" / "global_audit"
    lat, lon = np.arange(90, -90.1, -5.0), np.arange(-180, 180, 5.0)
    for m in range(1, 13):
        times = month_times(2000, m)[:40]
        write_store(sub / f"diagnostics_glob_2000-{m:02d}.zarr", lat, lon, times, ORIGINAL_AUDIT_ATTRS)

    subprocess.run([sys.executable, str(ORIGINAL / "ada" / "tail_thresholds.py"), "--base", str(base),
                    "--derived-subdir", "derived/global_audit", "--no-archive"], check=True,
                   stdout=subprocess.DEVNULL)
    old = sorted((base / "calibration").glob("thresholds_*.json"))[-1]
    new = work / "new_thresholds.json"
    T.calibrate(sub, new, "audit")

    a = json.loads(old.read_text())["thresholds"]
    b = json.loads(new.read_text())["thresholds"]
    exact = a == b
    print(f"\nSTAGE 3: original vs essentials thresholds identical: {exact}")

    # tail method against a whole-array sort
    worst = 0.0
    for k in C.DIAGNOSTIC_KEYS:
        vals = np.concatenate([xr.open_zarr(s)[k].transpose("latitude", "longitude", "time").values.ravel()
                               for s in sorted(sub.glob("*.zarr"))])
        wts = np.tile(np.repeat(T.latitude_weights(lat), lon.size * 40), 12)
        full = T.weighted_percentile(vals, wts, list(C.SEVERITIES.values()))
        tail = np.array([b[k][s] for s in C.SEVERITIES])
        worst = max(worst, float(np.max(np.abs(full - tail) / np.abs(full))))
    print(f"STAGE 3: tail method vs whole-array sort, largest relative difference {worst:.2e}")
    return exact and worst < 1e-6, a


def stage4(work: Path, thresholds: dict) -> bool:
    spec = importlib.util.spec_from_file_location("ftc", ORIGINAL / "ada" / "full_trend_check.py")
    ftc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ftc)
    import aggregate

    derived = work / "derived4"
    lat, lon = np.arange(60, 29.9, -1.0), np.arange(-75, 0.1, 1.0)
    years = [1979, 1980, 1981]
    for y in years:
        for m in range(1, 13):
            write_store(derived / f"diagnostics_na_{y}-{m:02d}.zarr", lat, lon, month_times(y, m),
                        ORIGINAL_AUDIT_ATTRS)

    names = C.DIAGNOSTIC_KEYS
    signs = {k: "+" for k in names}
    worst = 0.0
    for y in years:
        ds = ftc.subset_box(ftc.load_year(derived.parent, y, derived.name), **ftc.PROSSER_BOX)
        w = ftc.lat_weights_for(ds)
        new = TR.year_frequencies(TR.open_year(derived, y, "audit"), thresholds, list(C.SEASON_MONTHS), names)
        for season, months in C.SEASON_MONTHS.items():
            part = ds.isel(time=ds["time"].dt.month.isin(months).values)
            fields = {k: part[k] for k in names}
            per_sev = aggregate.exceedance_mean_all_severities(fields, thresholds, list(C.SEVERITIES), signs)
            for sev in C.SEVERITIES:
                old_ens = ftc.weighted_rate(per_sev[sev]["exceedance_mean"], w)
                worst = max(worst, abs(old_ens - new[season][sev]["ensemble"]) / old_ens)
                for k in names:
                    old_k = ftc.weighted_rate(aggregate.exceedance_field(fields[k], thresholds[k][sev], "+"), w)
                    if old_k > 0:
                        worst = max(worst, abs(old_k - new[season][sev][k]) / old_k)
    print(f"STAGE 4: frequencies, largest relative difference to the original {worst:.2e}")

    # OLS: same numbers as the original ols()
    x = np.arange(1979, 2021, dtype=float)
    y = 0.01 + 1e-4 * (x - 1979) + rng.normal(0, 5e-4, x.size)
    a, b = ftc.ols(x, y), TR.ols(x, y)
    ols_same = all(np.isclose(a[k], b[k], rtol=0, atol=0) for k in ("slope", "intercept", "se", "t", "half", "r2"))
    print(f"STAGE 4: OLS identical to the original: {ols_same}")
    return worst < 1e-12 and ols_same


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp())
    ok3, thr = stage3(work)
    ok4 = stage4(work, thr)
    print("\nALL CHECKS PASSED" if ok3 and ok4 else "\nCHECKS FAILED")
    return 0 if ok3 and ok4 else 1


if __name__ == "__main__":
    raise SystemExit(main())
