"""
trends.py -- Stage 4: exceedance, seasonal turbulence frequency, 1979-2020 trends.

What is computed (Prosser et al. 2023, following Williams 2017)
---------------------------------------------------------------
1. Exceedance. For diagnostic D, severity s, grid cell c and time step t:
       E_{D,s}(c,t) = 1  if D(c,t) >= threshold_{D,s}   else 0
   (undefined where D is NaN). Thresholds come from Stage 3.

2. Ensemble. The 21 exceedance fields are averaged per cell and time step,
   ignoring diagnostics that are NaN there:
       Ebar_s(c,t) = sum_D E_{D,s}(c,t) / #(defined D at c,t)
   Exceedance first, average second -- averaging raw values first would be a
   different (and wrong) method.

3. Frequency. Over Prosser's box (36-60N, 55-10W) and all time steps of a
   season in year Y, weighted by cell area w = cos(latitude):
       P_s(Y) = sum_{c,t} w_c Ebar_s(c,t) / sum_{c,t: defined} w_c
   Multiply by 24 x days in the season for "hours at an average point",
   which is what Prosser's Table 1 reports. The same is done for every
   diagnostic on its own (E_{D,s} instead of Ebar_s).
   Seasons: DJF(Y) = Jan, Feb and Dec of the same calendar year Y.

4. Trend. Ordinary least squares P_s(Y) = a + b Y over Y = 1979..2020 (n = 42).
       fitted change   = (P(2020) - P(1979)) / P(1979)       (fitted values)
       95 % interval   = change +/- t_{0.975, n-2} * se(b) * 41 / P(1979)
   and it is compared with Prosser's published relative change.

Memory: one year of one diagnostic is held at a time; the exceedance counts
per cell are kept as small integers. A year of all 21 diagnostics is never in
memory at once.
"""
from __future__ import annotations

import calendar
import csv
import json
from pathlib import Path

import numpy as np
import xarray as xr

import config as C
from thresholds import load_thresholds, original_store_convention

# Two-sided 95 % Student-t critical values (same lookup as the original).
_TCRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 10: 2.228, 15: 2.131,
          20: 2.086, 25: 2.060, 30: 2.042, 40: 2.021, 60: 2.000, 120: 1.980}


def t_critical(df: int) -> float:
    if df <= 0:
        return float("nan")
    if df in _TCRIT:
        return _TCRIT[df]
    return 1.960 if df > 120 else _TCRIT[min(k for k in _TCRIT if k >= df)]


def season_days(year: int, season: str) -> int:
    feb = 29 if calendar.isleap(year) else 28
    return {"djf": 31 + feb + 31, "mam": 92, "jja": 92, "son": 91,
            "annual": 337 + feb}[season]


# ---------------------------------------------------------------------------
# Reading one year of the North Atlantic series
# ---------------------------------------------------------------------------
def store_convention(store: Path) -> str:
    attrs = json.loads((Path(store) / ".zattrs").read_text())
    return attrs.get("convention") or original_store_convention(attrs)


def open_year(derived: Path, year: int, convention_name: str) -> xr.Dataset:
    """The 12 monthly stores of one year, lazily, cut to Prosser's box."""
    stores = [Path(derived) / C.month_store_name("north_atlantic", year, m) for m in range(1, 13)]
    missing = [s.name for s in stores if not s.exists()]
    if missing:
        raise FileNotFoundError(f"{year}: missing {missing}")
    wrong = {s.name: store_convention(s) for s in stores if store_convention(s) != convention_name}
    if wrong:
        raise SystemExit(f"stores not built under '{convention_name}': {wrong}")
    ds = xr.concat([xr.open_zarr(s) for s in stores], dim="time")
    return subset_box(ds, **C.PROSSER_BOX)


def subset_box(ds: xr.Dataset, lat, lon) -> xr.Dataset:
    """Cut a lat/lon box, whichever way the coordinates run (ERA5 latitudes descend)."""
    la, lo = ds["latitude"].values, ds["longitude"].values
    lat_slice = slice(lat[1], lat[0]) if la[0] > la[-1] else slice(lat[0], lat[1])
    lon_slice = slice(lon[1], lon[0]) if lo[0] > lo[-1] else slice(lon[0], lon[1])
    return ds.sel(latitude=lat_slice, longitude=lon_slice)


# ---------------------------------------------------------------------------
# Exceedance frequencies for one year
# ---------------------------------------------------------------------------
def year_frequencies(box: xr.Dataset, thresholds: dict, seasons, names) -> dict:
    """{season: {severity: {"ensemble": P, name: P_D, ...}}} for one year."""
    w = np.cos(np.deg2rad(box["latitude"].values))                  # (lat,)
    months = box["time"].dt.month.values
    masks = {s: np.isin(months, C.SEASON_MONTHS[s]) for s in seasons}
    shape = (box.sizes["latitude"], box.sizes["longitude"])
    n_exceed = {(s, sev): np.zeros(shape + (int(masks[s].sum()),), np.uint8)
                for s in seasons for sev in C.SEVERITIES}
    n_valid = {s: np.zeros(shape + (int(masks[s].sum()),), np.uint8) for s in seasons}
    out = {s: {sev: {} for sev in C.SEVERITIES} for s in seasons}

    for name in names:
        x_all = box[name].transpose("latitude", "longitude", "time").values
        for s in seasons:
            x = x_all[:, :, masks[s]]
            valid = np.isfinite(x)
            n_valid[s] += valid
            w_valid = float(np.sum(w * valid.sum(axis=(1, 2))))
            for sev in C.SEVERITIES:
                exceed = valid & (x >= thresholds[name][sev])        # sign "+" for all 21
                n_exceed[(s, sev)] += exceed
                out[s][sev][name] = float(np.sum(w * exceed.sum(axis=(1, 2)))) / w_valid
        del x_all

    for s in seasons:
        defined = n_valid[s] > 0
        w_defined = float(np.sum(w * defined.sum(axis=(1, 2))))
        for sev in C.SEVERITIES:
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = np.where(defined, n_exceed[(s, sev)] / n_valid[s], 0.0)
            out[s][sev]["ensemble"] = float(np.sum(mean * w[:, None, None])) / w_defined
    return out


# ---------------------------------------------------------------------------
# Trend
# ---------------------------------------------------------------------------
def ols(x, y) -> dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    df = x.size - 2
    sxx = float(((x - x.mean()) ** 2).sum())
    ss_res, ss_tot = float((resid ** 2).sum()), float(((y - y.mean()) ** 2).sum())
    se = float(np.sqrt(ss_res / df / sxx)) if df > 0 and sxx > 0 else float("nan")
    return {"slope": float(slope), "intercept": float(intercept), "se": se, "df": df,
            "t": float(slope / se) if se > 0 else float("inf"),
            "half": t_critical(df) * se, "r2": 1 - ss_res / ss_tot if ss_tot else float("nan")}


def trend_row(years, series) -> dict:
    fit = ols(years, series)
    y0 = fit["slope"] * years[0] + fit["intercept"]
    y1 = fit["slope"] * years[-1] + fit["intercept"]
    span = float(years[-1] - years[0])
    change = (y1 - y0) / y0 if y0 > 0 else float("nan")
    half = span * fit["half"] / y0 if y0 > 0 else float("nan")
    return {"level_first": y0, "level_last": y1, "change": change,
            "ci_lo": change - half, "ci_hi": change + half,
            "t": fit["t"], "r2": fit["r2"],
            "significant": abs(fit["t"]) >= t_critical(fit["df"])}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def threshold_convention(provenance: dict, stated) -> str:
    """The convention a threshold file belongs to (older files do not say)."""
    if stated:
        return stated
    period = str(provenance.get("period", ""))
    if "global_audit" in period:
        return "audit"
    if "derived/global" in period:
        return "baseline"
    raise SystemExit(f"cannot tell which convention these thresholds belong to (period={period!r})")


def run(derived: Path, thresholds_path: Path, convention_name: str, out_dir: Path,
        seasons=tuple(C.SEASON_MONTHS), years=None) -> dict:
    years = list(years or range(C.TREND_YEARS[0], C.TREND_YEARS[1] + 1))
    thresholds, stated, prov = load_thresholds(thresholds_path)
    tconv = threshold_convention(prov, stated)
    if tconv != convention_name:
        raise SystemExit(f"thresholds are '{tconv}', series is '{convention_name}': refusing to mix")
    names = [k for k in C.DIAGNOSTIC_KEYS if k in thresholds]
    print(f">>> series {derived}\n>>> thresholds {thresholds_path} ({tconv})\n"
          f">>> {len(years)} years, seasons {list(seasons)}, {len(names)} diagnostics")

    rates = {s: {sev: {} for sev in C.SEVERITIES} for s in seasons}   # -> series -> [per year]
    for year in years:
        box = open_year(derived, year, convention_name)
        for s in seasons:
            expected = season_days(year, s) * C.STEPS_PER_DAY
            got = int(np.isin(box["time"].dt.month.values, C.SEASON_MONTHS[s]).sum())
            if got != expected:
                print(f"   !! {year} {s}: {got} time steps, expected {expected}")
        freq = year_frequencies(box, thresholds, seasons, names)
        for s in seasons:
            for sev in C.SEVERITIES:
                for series, p in freq[s][sev].items():
                    rates[s][sev].setdefault(series, []).append(p)
        print(f"   {year} done", flush=True)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{convention_name}_{years[0]}-{years[-1]}"
    with open(out_dir / f"frequencies_{tag}.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["season", "severity", "series", *years])
        for s in seasons:
            for sev in C.SEVERITIES:
                for series, vals in rates[s][sev].items():
                    wr.writerow([s, sev, series, *[f"{v:.10g}" for v in vals]])

    trends = {}
    with open(out_dir / f"trends_{tag}.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["season", "severity", "series", "level_first", "level_last", "change",
                     "ci_lo", "ci_hi", "t", "r2", "significant", "prosser_change", "prosser_inside_ci"])
        for s in seasons:
            for sev in C.SEVERITIES:
                for series, vals in rates[s][sev].items():
                    row = trend_row(np.array(years, float), vals)
                    p = C.PROSSER_TABLE1[s][sev][2] if series == "ensemble" else ""
                    inside = (row["ci_lo"] <= p <= row["ci_hi"]) if p != "" else ""
                    trends[(s, sev, series)] = {**row, "prosser": p, "inside": inside}
                    wr.writerow([s, sev, series, *[f"{row[k]:.8g}" for k in
                                 ("level_first", "level_last", "change", "ci_lo", "ci_hi", "t", "r2")],
                                 row["significant"], p, inside])
    print_summary(trends, seasons, years)
    print(f"\n    wrote {out_dir / f'frequencies_{tag}.csv'}\n    wrote {out_dir / f'trends_{tag}.csv'}")
    return trends


def print_summary(trends, seasons, years):
    print(f"\n{'season':<8}{'sev':<6}{'1979':>9}{'2020':>9}{'ours':>7}{'95% CI':>15}"
          f"{'Prosser':>9}{'ratio':>7}{'t':>7}{'sig':>5}{'in CI':>7}")
    n_sig = n_in = 0
    for s in seasons:
        for sev in C.SEVERITIES:
            r = trends[(s, sev, "ensemble")]
            n_sig += r["significant"]
            n_in += bool(r["inside"])
            print(f"{s:<8}{C.SEVERITY_LABEL[sev]:<6}{r['level_first']:>9.3%}{r['level_last']:>9.3%}"
                  f"{r['change']:>7.0%}   [{r['ci_lo']:>4.0%},{r['ci_hi']:>4.0%}]"
                  f"{r['prosser']:>9.0%}{r['change'] / r['prosser']:>7.2f}{r['t']:>7.2f}"
                  f"{'yes' if r['significant'] else 'no':>5}{'yes' if r['inside'] else 'no':>7}")
            hours0 = r["level_first"] * season_days(years[0], s) * 24
            hours1 = r["level_last"] * season_days(years[-1], s) * 24
            p0, p1, _ = C.PROSSER_TABLE1[s][sev]
            print(f"{'':<14}hours {hours0:7.1f} -> {hours1:7.1f}   Prosser {p0:7.1f} -> {p1:7.1f}")
    total = len(seasons) * len(C.SEVERITIES)
    print(f"\nsignificant at 5 %: {n_sig}/{total}    Prosser inside 95 % CI: {n_in}/{total}")
