#!/usr/bin/env python
"""
ada/ubf_stencil_test.py
=======================
(UBF-STENCIL-2026-11-20)
Is `ubf`'s disagreement with Prosser an artefact of HOW the residual is
discretised? Three versions of the same formula, side by side, on data already
on the share. Nothing is written outside --out; no production store is touched.

    pixi run python ada/ubf_stencil_test.py            # paths default to the share

THE QUESTION
    `ubf` is the one class-C diagnostic whose disagreement with Prosser was put
    down to our own code (FINAL_DATASET.md L10). Its 1979 North Atlantic level is
    0.39 / 0.19 / 0.08 of his (LOG / MOG / SOG), and its reference-year tail is
    POLAR (56 % of the SOG mass poleward of 60 deg), where Williams & Storer
    (2022, p. 1428-9) describe UBF as stronger in the Northern Hemisphere.

    The formula (Sharman A30) and the spherical geometry are both verified, the
    geometry now also on an exact balanced flow over the poles (2026-11-20, see
    stage 0 below). What is NOT balance-preserving is the discretisation:

        grad^2 Phi = D(D Phi)    two first differences  -> response (sin kh/kh)^2
        f zeta     = f D(v)      one first difference   -> response  sin kh/kh

    so even an exactly geostrophic wave leaves a NUMERICAL residual, 36 % of
    f zeta at 4 grid steps, 10 % at 8, 2.6 % at 16. At 0.25 deg ERA5 carries
    real energy at 4-8 grid steps, and FORMULA_AUDIT.md section 10.5 measured
    the residual at 73 % of its largest term on real data -- for a balance
    residual at 200 hPa that is large. If the mismatch is concentrated in
    particular regions it takes over a GLOBAL percentile tail, which would
    explain both the polar tail and a NA deficit that deepens with severity.

THE THREE VERSIONS (only ubf changes; the other 20 diagnostics are untouched)
    V0  current      grad^2 = div(grad) with 2nd-order centred differences.
                     Must reproduce the production stores (control, stage 1).
    V1  compact      grad^2 with the compact 3-point second difference
                     (MetPy's `laplacian` convention); J and zeta unchanged.
    V2  4th-order    every horizontal derivative in ubf with the 4th-order
                     centred difference (-f2 + 8f1 - 8f-1 + f-2)/12h; same
                     construction as V0. Mismatch falls to 15 / 1.2 / 0.1 %
                     at 4 / 8 / 16 grid steps.
    All three share the audit conventions (map dy, computed zeta), the same
    map-scale factors and curvature terms, and the same edge handling.

STAGES
    0  analytic: each version on an exact balanced flow (Williamson et al. 1992
       test case 2, tilted over the poles) and on balanced plane waves.
    1  CONTROL: V0 on the first global file against
       derived/global_audit/diagnostics_glob_2000-01.zarr at the same times.
       Exit 2 if it does not reproduce production -- nothing below is then
       interpretable.
    2  global reference year, the 48-day sub-sample (raw/global_sub48, days
       1/9/17/25 of 2000): cos-weighted thresholds at the 5 severities, zonal
       tail profile, NH/SH ratio, Prosser-box density, flip rate vs V0.
    3  North Atlantic 1979 (12 months): annual exceedance level in Prosser's
       box under each version's own thresholds, against his 1979 level.

PRE-REGISTERED, BEFORE THE NUMBERS
    P1  V0 reproduces production: Spearman >= 0.99999 against the zarr.
    P2  V2's numerical residual on balanced 4-8 grid-step waves is at least 2x
        smaller than V0's; V1 lies between at 6-8 steps.
    P3  IF the operator mismatch drives the polar tail, then under V2, against V0:
          - Prosser-box density rises at LOG, MOG and SOG, most at SOG;
          - the polar-cap densities fall and NH/SH at SOG rises above 1;
          - the NA 1979 level ratio vs Prosser moves toward 1 at all three.
    P4  FALSIFIED if V2's box density is within +-20 % of V0's at both MOG and
        SOG. Then the mismatch is not the cause; ubf stays class C with its
        cause recorded as unknown, and this line of enquiry stops.
    No class changes on this test alone. If P3 holds, the next step is a
    ubf-only 42-year re-derive with the chosen version (504 NA months + the
    full reference year), scored against Prosser S4 at three severities, and
    only then the class rule of FINAL_DATASET.md section 4.3.

    This script prints statistics, not verdicts (STATUS 12.8, 15.6).
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "essentials"))
sys.path.insert(0, str(HERE))
import config as C                       # noqa: E402  essentials/config.py
import diagnostics as E                  # noqa: E402  essentials/diagnostics.py

VERSION = "UBF-STENCIL-2026-11-20"
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
PROSSER_BOX = dict(lat=(36.0, 60.0), lon=(-55.0, -10.0))      # Prosser Fig. 2
SEVERITIES = {"light": 97.0, "light_moderate": 99.1, "moderate": 99.6,
              "moderate_severe": 99.8, "severe": 99.9}
REPORT = ("light", "moderate", "severe")
# Prosser (2023) Fig. S4-LOG, Fig. 4, Fig. S4-SOG: ubf 1979 level, % (read off his
# figures; ada/prosser_published*.py), and our production fitted level ratio
# (FINAL_DATASET.md section 5, audit series).
PROSSER_1979 = {"light": 6.3, "moderate": 0.72, "severe": 0.135}
PRODUCTION_RATIO = {"light": 0.39, "moderate": 0.19, "severe": 0.08}
BANDS = [("90S-60S", -90.1, -60), ("60S-30S", -60, -30), ("30S-0", -30, 0),
         ("0-30N", 0, 30), ("30N-60N", 30, 60), ("60N-90N", 60, 90.1)]
VERSIONS = ("V0", "V1", "V2")


# ===========================================================================
# Operators. V0 uses essentials' own grid_difference, which reproduces rojak.
# ===========================================================================
def diff4(f: xr.DataArray, spacing_m: np.ndarray, dim: str) -> xr.DataArray:
    """4th-order centred first difference along `dim` on a UNIFORM grid.

    Interior: (-f[i+2] + 8 f[i+1] - 8 f[i-1] + f[i-2]) / (12 h).
    The two rows/columns next to each edge fall back to np.gradient (2nd-order
    centred one in, 1st-order one-sided at the edge), i.e. exactly V0's edges.
    """
    h = np.asarray(spacing_m, dtype=float)
    if not np.allclose(h, h[0], rtol=1e-9, atol=0.0):
        raise ValueError(f"diff4 needs uniform spacing along {dim}")
    out = E.grid_difference(f, spacing_m, dim)          # edges + dtype as V0
    ax = f.get_axis_num(dim)
    a = np.moveaxis(f.values, ax, 0)
    o = np.moveaxis(out.values.copy(), ax, 0)
    o[2:-2] = (-a[4:] + 8 * a[3:-1] - 8 * a[1:-3] + a[:-4]) / (12.0 * h[0])
    return out.copy(data=np.moveaxis(o, 0, ax))


def second_diff_compact(f: xr.DataArray, spacing_m: np.ndarray, dim: str) -> xr.DataArray:
    """Compact 3-point second difference (f[i+1] - 2 f[i] + f[i-1]) / h^2.
    Edge rows/columns: V0's D(D f), so only the interior differs."""
    h = np.asarray(spacing_m, dtype=float)
    if not np.allclose(h, h[0], rtol=1e-9, atol=0.0):
        raise ValueError(f"compact second difference needs uniform spacing along {dim}")
    dd = E.grid_difference(E.grid_difference(f, spacing_m, dim), spacing_m, dim)
    ax = f.get_axis_num(dim)
    a = np.moveaxis(f.values, ax, 0)
    o = np.moveaxis(dd.values.copy(), ax, 0)
    o[1:-1] = (a[2:] - 2 * a[1:-1] + a[:-2]) / (h[0] ** 2)
    return dd.copy(data=np.moveaxis(o, 0, ax))


def wind_derivs(u, v, geo, diff):
    """essentials.wind_derivatives with the first difference as a parameter."""
    k_x, k_y = geo["k_x"], geo["k_y"]
    x_curv = (k_y / k_x) * geo["dkx_dY"]
    y_curv = (k_x / k_y) * geo["dky_dX"]
    return {
        "du_dx": k_x * diff(u, geo["dx"], "longitude") - v * x_curv,
        "du_dy": k_y * diff(u, geo["dy"], "latitude") + v * y_curv,
        "dv_dx": k_x * diff(v, geo["dx"], "longitude") + u * x_curv,
        "dv_dy": k_y * diff(v, geo["dy"], "latitude") - u * y_curv,
    }


def ubf_versions(u, v, phi, geo, f=None, f_x=None, f_y=None, return_terms=False):
    """|grad^2 Phi - 2 J(u,v) - f zeta + beta u| under V0, V1, V2.

    u, v, phi at ONE level, dims (latitude, longitude[, time]). zeta is computed
    from u, v with the same operator as J (audit convention, F12).
    f defaults to 2 Omega sin(phi) and the beta u term to 2 Omega cos(phi)/R * u,
    exactly as production; stage 0 passes a rotated f with its gradient instead.
    """
    lat_rad = np.deg2rad(u["latitude"])
    production_f = f is None
    if production_f:   # exactly essentials.ubf / 2_diagnostics.ubf, same order of operations
        f = (2 * C.OMEGA * np.sin(lat_rad)).broadcast_like(u.isel(time=0) if "time" in u.dims else u)
        beta = (2 * C.OMEGA * np.cos(lat_rad) / C.R_EARTH).broadcast_like(u.isel(time=0) if "time" in u.dims else u)

    out, terms = {}, {}
    D2 = E.grid_difference
    for name in VERSIONS:
        diff = diff4 if name == "V2" else D2
        wd = wind_derivs(u, v, geo, diff)
        jac = wd["du_dx"] * wd["dv_dy"] - wd["du_dy"] * wd["dv_dx"]
        zeta = wd["dv_dx"] - wd["du_dy"]
        if name == "V1":
            k_x, k_y = geo["k_x"], geo["k_y"]
            phi_X = D2(phi, geo["dx"], "longitude")
            phi_Y = D2(phi, geo["dy"], "latitude")
            dky_dY = D2(k_y, geo["dy"], "latitude")
            x_curv = (k_y / k_x) * geo["dkx_dY"]
            y_curv = (k_x / k_y) * geo["dky_dX"]
            lap = (k_x ** 2 * second_diff_compact(phi, geo["dx"], "longitude")
                   - (k_y * phi_Y) * x_curv
                   + k_y * (k_y * second_diff_compact(phi, geo["dy"], "latitude") + dky_dY * phi_Y)
                   - (k_x * phi_X) * y_curv)
        else:
            gp = wind_derivs(diff(phi, geo["dx"], "longitude") * geo["k_x"],
                             diff(phi, geo["dy"], "latitude") * geo["k_y"], geo, diff)
            lap = gp["du_dx"] + gp["dv_dy"]
        if production_f:
            out[name] = np.abs(lap - 2 * jac - f * zeta + beta * u)
            cor = f * zeta - beta * u
        else:   # div(f grad psi) = f zeta + f_x v - f_y u  (rotated f, stage 0)
            cor = f * zeta + f_x * v - f_y * u
            out[name] = np.abs(lap - 2 * jac - cor)
        if return_terms:
            terms[name] = dict(lap=lap, J2=2 * jac, cor=cor)
    return (out, terms) if return_terms else out


# ===========================================================================
# Weighted statistics, without a full-size weights array
# ===========================================================================
def tail_thresholds(x: np.ndarray, w_lat: np.ndarray, pcts) -> np.ndarray:
    """cos-weighted Hazen percentiles for x of shape (lat, ...), exact for any
    percentile above the retained cut. Same arithmetic as
    calib_weighted_percentile.weighted_percentile (and ada/tail_thresholds.py):
    P_i = (W_below + C_i - 0.5 w_i) / W_total over the sorted tail."""
    pcts = np.asarray(pcts, dtype=float)
    n_per_row = np.prod(x.shape[1:])
    finite_rows = np.isfinite(x).reshape(x.shape[0], -1).sum(axis=1)
    w_total = float((w_lat * finite_rows).sum())
    sample = x.reshape(-1)[:: max(1, x.size // 2_000_000)]
    for q in (90.0, 85.0, 75.0, 50.0, 0.0):
        cut = np.nanpercentile(sample, q) if q > 0 else np.nanmin(x)
        above = x >= cut
        cnt = above.reshape(x.shape[0], -1).sum(axis=1)
        w_above = float((w_lat * cnt).sum())
        if w_above / w_total > (100.0 - pcts.min()) / 100.0 * 1.2:
            break
    vals = x[above]
    lat_idx = np.broadcast_to(np.arange(x.shape[0]).reshape((-1,) + (1,) * (x.ndim - 1)),
                              x.shape)[above]
    order = np.argsort(vals, kind="stable")
    v = vals[order].astype(np.float64)
    w = w_lat[lat_idx[order]].astype(np.float64)
    pos = ((w_total - w_above) + np.cumsum(w) - 0.5 * w) / w_total
    lo = pos[0]
    if lo > pcts.min() / 100.0:
        raise RuntimeError("tail cut above the lowest requested percentile")
    del n_per_row
    return np.interp(pcts / 100.0, pos, v)


def lat_counts(mask: np.ndarray) -> np.ndarray:
    return mask.reshape(mask.shape[0], -1).sum(axis=1)


def spearman_sub(a: np.ndarray, b: np.ndarray, n: int = 2_000_000, seed: int = 0) -> float:
    rng = np.random.default_rng(seed)
    idx = rng.choice(a.size, size=min(n, a.size), replace=False)
    a, b = a.reshape(-1)[idx].astype(np.float64), b.reshape(-1)[idx].astype(np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    ra = np.argsort(np.argsort(a[ok])); rb = np.argsort(np.argsort(b[ok]))
    return float(np.corrcoef(ra, rb)[0, 1])


def peak_rss_gb() -> float:
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    except ImportError:
        return float("nan")


# ===========================================================================
# Stage 0 -- analytic
# ===========================================================================
def stage0() -> dict:
    print("=" * 78 + "\nSTAGE 0  analytic checks\n" + "=" * 78)
    res = {}
    # (a) balanced plane waves, f-plane at 45N, along longitude on a 0.25 deg grid
    print("  (a) geostrophic plane wave, numerical residual as a share of |f zeta|")
    print("      wavelength (grid steps)   " + "   ".join(f"{n:>5}" for n in VERSIONS))
    lat0 = 45.0
    for nsteps in (4, 6, 8, 12, 16, 32):
        lon = np.arange(0, 0.25 * 32 * nsteps, 0.25) - 180.0
        lat = np.array([lat0 + 0.5, lat0 + 0.25, lat0, lat0 - 0.25, lat0 - 0.5])
        geo = E.grid_geometry(lat, lon, "map")
        hx = float(geo["dx"][0] / geo["k_x"].values[2, 0])     # physical spacing = map spacing / k_x
        k = 2 * np.pi / (nsteps * hx)
        f0 = 2 * C.OMEGA * np.sin(np.deg2rad(lat0))
        X = np.arange(lon.size) * hx
        co = {"latitude": lat, "longitude": lon}
        da = lambda a: xr.DataArray(np.tile(a, (lat.size, 1)), dims=("latitude", "longitude"), coords=co)
        phi = da(1000.0 * np.cos(k * X))
        v = da(-1000.0 * k * np.sin(k * X) / f0)
        u = xr.zeros_like(v)
        f = xr.full_like(v, f0)
        z0 = xr.zeros_like(v)
        out, terms = ubf_versions(u, v, phi, geo, f=f, f_x=z0, f_y=z0, return_terms=True)
        sl = dict(latitude=2, longitude=slice(3 * nsteps, -3 * nsteps))
        row = {}
        for name in VERSIONS:
            r = out[name].isel(**sl).values
            s = np.abs(terms[name]["cor"].isel(**sl).values)
            row[name] = float(r.max() / s.max())
        res[f"wave_{nsteps}"] = row
        print(f"      {nsteps:>5}                     " + "   ".join(f"{row[n]:5.3f}" for n in VERSIONS))

    # (b) Williamson TC2 tilted over the poles, global 0.25 deg; exact value of the
    #     Cartesian-form residual on the sphere is -u0^2 cos^2(phi')/a^2.
    a, Om, u0, gh0, al = C.A_WGS84, C.OMEGA, 40.0, 2.94e4, np.pi / 2
    lat = np.arange(90, -90.001, -0.25); lon = np.arange(-180, 180, 0.25)
    P, L = np.meshgrid(np.deg2rad(lat), np.deg2rad(lon), indexing="ij")
    co = {"latitude": lat, "longitude": lon}
    da = lambda x: xr.DataArray(x, dims=("latitude", "longitude"), coords=co)
    s_rot = -np.cos(L) * np.cos(P) * np.sin(al) + np.sin(P) * np.cos(al)
    u = da(u0 * (np.cos(P) * np.cos(al) + np.sin(P) * np.cos(L) * np.sin(al)))
    v = da(-u0 * np.sin(L) * np.sin(al))
    phi = da(gh0 - (a * Om * u0 + u0 ** 2 / 2) * s_rot ** 2)
    f = da(2 * Om * s_rot)
    geo = E.grid_geometry(lat, lon, "map")
    f_x, f_y = E.d_dx(f, geo), E.d_dy(f, geo)
    out = ubf_versions(u, v, phi, geo, f=f, f_x=f_x, f_y=f_y)
    exact = np.abs(-(u0 ** 2 / a ** 2) * (1 - s_rot ** 2))
    inner = dict(longitude=slice(2, -2))
    print("\n  (b) exact balanced flow over the poles (Williamson TC2, alpha=90):")
    print("      max |error| / 3.6e-9 (median real UBF), seam columns excluded")
    print("      |lat| band      " + "   ".join(f"{n:>7}" for n in VERSIONS))
    for lo_, hi_ in [(0, 60), (60, 85), (85, 89.5), (89.5, 90)]:
        m = (np.abs(lat) >= lo_) & (np.abs(lat) < hi_)
        row = {n: float(np.abs(out[n].values - exact)[m][:, 2:-2].max() / 3.6e-9) for n in VERSIONS}
        res[f"tc2_{lo_}_{hi_}"] = row
        print(f"      {lo_:>4}-{hi_:<6}     " + "   ".join(f"{row[n]:7.4f}" for n in VERSIONS))
    del inner
    return res


# ===========================================================================
# Loading
# ===========================================================================
def load_200(path) -> xr.Dataset:
    ds = E.prepare(E.open_era5(path))
    ds = ds[["eastward_wind", "northward_wind", "geopotential"]].sel(
        pressure_level=C.TARGET_LEVEL_HPA)
    return ds


def per_day(ds: xr.Dataset, steps: int = 8):
    for i in range(0, ds.sizes["time"], steps):
        yield ds.isel(time=slice(i, i + steps)).load()


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]


# ===========================================================================
# Stage 1 -- control
# ===========================================================================
def stage1(global_file: Path, control_zarr: Path) -> dict:
    print("\n" + "=" * 78 + "\nSTAGE 1  CONTROL: V0 against the production store\n" + "=" * 78)
    ds = load_200(global_file)
    z = xr.open_zarr(control_zarr)["ubf"]
    tdim = next(d for d in ("time", "valid_time") if d in z.dims)
    if tdim != "time":
        z = z.rename({tdim: "time"})
    if float(z["longitude"].max()) > 180.0:      # same normalisation as essentials.prepare
        z = z.assign_coords(longitude=(((z["longitude"] + 180) % 360) - 180)).sortby("longitude")
    times = np.intersect1d(ds["time"].values, z["time"].values)[:8]
    if times.size == 0:
        print("!! no common timesteps between the global file and the control store")
        return {"ok": False}
    day = ds.sel(time=times).load()
    geo = E.grid_geometry(day["latitude"], day["longitude"], "map")
    mine = ubf_versions(day["eastward_wind"], day["northward_wind"], day["geopotential"], geo)["V0"]
    ref = z.sel(time=times).load()
    mine, ref = xr.align(mine, ref.transpose(*mine.dims), join="inner")
    # The longitude seam (one-sided differences at the array edge) may sit at a
    # different meridian in the production store; exclude 2 columns either side
    # of 0 and +-180 deg so the control tests the operators, not the seam.
    lonv = mine["longitude"].values
    far = (np.abs(lonv) > 0.6) & (np.abs(np.abs(lonv) - 180.0) > 0.6)
    mine, ref = mine.isel(longitude=far), ref.isel(longitude=far)
    a = mine.values.astype(np.float64); b = ref.values.astype(np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    rel = np.abs(a[ok] - b[ok]) / np.maximum(np.abs(b[ok]), 1e-30)
    rho = spearman_sub(a, b)
    out = dict(ok=bool(rho >= 0.99999), n_steps=int(times.size), spearman=rho,
               rel_median=float(np.median(rel)), rel_p99=float(np.percentile(rel, 99)),
               rel_max=float(rel.max()))
    print(f"  {times.size} timesteps, {ok.sum():,} cells")
    print(f"  Spearman {rho:.7f}   relative difference median {out['rel_median']:.2e}  "
          f"p99 {out['rel_p99']:.2e}  max {out['rel_max']:.2e}")
    print("  P1 requires Spearman >= 0.99999" + ("" if out["ok"] else "   -- NOT MET, stopping"))
    return out


# ===========================================================================
# Stage 2 -- global reference-year sub-sample
# ===========================================================================
def stage2(files) -> tuple[dict, dict]:
    print("\n" + "=" * 78 + "\nSTAGE 2  global reference year (48-day sub-sample)\n" + "=" * 78)
    chunks = {n: [] for n in VERSIONS}
    lat = None
    t0 = time.time()
    for k, path in enumerate(files):
        ds = load_200(path)
        for day in per_day(ds):
            if lat is None:
                lat = day["latitude"].values.astype(np.float64)
                lon = day["longitude"].values.astype(np.float64)
                geo = E.grid_geometry(day["latitude"], day["longitude"], "map")
            out = ubf_versions(day["eastward_wind"], day["northward_wind"], day["geopotential"], geo)
            for n in VERSIONS:
                chunks[n].append(out[n].transpose("latitude", "longitude", "time").values.astype(np.float32))
        print(f"  [{k + 1:>2}/{len(files)}] {Path(path).name}  "
              f"{(time.time() - t0) / 60:5.1f} min  rss {peak_rss_gb():.1f} GB")
    data = {n: np.concatenate(chunks[n], axis=2) for n in VERSIONS}
    del chunks
    w_lat = np.cos(np.deg2rad(lat)); w_lat[np.abs(lat) >= 90.0] = 0.0
    pcts = list(SEVERITIES.values())
    thr = {n: dict(zip(SEVERITIES, map(float, tail_thresholds(data[n], w_lat, pcts)))) for n in VERSIONS}

    lon_in = (lon >= PROSSER_BOX["lon"][0]) & (lon <= PROSSER_BOX["lon"][1])
    lat_in = (lat >= PROSSER_BOX["lat"][0]) & (lat <= PROSSER_BOX["lat"][1])
    area_band = np.array([w_lat[(lat >= lo) & (lat < hi)].sum() for _, lo, hi in BANDS]) / w_lat.sum()
    nt = data["V0"].shape[2]
    w_all = w_lat.sum() * lon.size * nt
    w_box = w_lat[lat_in].sum() * lon_in.sum() * nt

    stats = {}
    exceed = {}
    for n in VERSIONS:
        for sev in REPORT:
            m = data[n] >= thr[n][sev]
            exceed[(n, sev)] = m
            mass_lat = lat_counts(m) * w_lat
            gfrac = mass_lat.sum() / w_all
            band = np.array([mass_lat[(lat >= lo) & (lat < hi)].sum() for _, lo, hi in BANDS]) / mass_lat.sum()
            box = float((m[lat_in][:, lon_in].reshape(lat_in.sum(), -1).sum(axis=1) * w_lat[lat_in]).sum())
            box_frac = box / w_box
            stats[(n, sev)] = dict(
                version=n, severity=sev, threshold=thr[n][sev], global_pct=100 * gfrac,
                **{f"density_{b}": float(d) for (b, _, _), d in zip(BANDS, band / area_band)},
                nh_sh=float(mass_lat[lat > 0].sum() / max(mass_lat[lat < 0].sum(), 1e-30)),
                box_pct=100 * box_frac, box_density=float(box_frac / gfrac))
    for n in ("V1", "V2"):
        rho = spearman_sub(data[n], data["V0"])
        for sev in REPORT:
            a, b = exceed[(n, sev)], exceed[("V0", sev)]
            flip = (lat_counts(a ^ b) * w_lat).sum() / max((lat_counts(a | b) * w_lat).sum(), 1e-30)
            stats[(n, sev)].update(spearman_vs_V0=rho, flip_vs_V0=float(flip))

    print(f"\n  {nt} timesteps x {lat.size} x {lon.size};  rss {peak_rss_gb():.1f} GB")
    print("  tail DENSITY by band (1.00 = even), NH/SH mass ratio, and Prosser-box density")
    print(f"  {'ver':<4}{'sev':<10}{'global%':>8}" + "".join(f"{b:>9}" for b, _, _ in BANDS)
          + f"{'NH/SH':>7}{'box':>7}{'rho':>9}{'flip':>7}")
    for sev in REPORT:
        for n in VERSIONS:
            s = stats[(n, sev)]
            print(f"  {n:<4}{sev:<10}{s['global_pct']:>8.3f}"
                  + "".join(f"{s['density_' + b]:>9.2f}" for b, _, _ in BANDS)
                  + f"{s['nh_sh']:>7.2f}{s['box_density']:>7.2f}"
                  + (f"{s['spearman_vs_V0']:>9.5f}{100 * s['flip_vs_V0']:>6.1f}%" if n != "V0" else ""))
        print()
    def rel(n, sev):
        b0 = stats[("V0", sev)]["box_density"]
        return stats[(n, sev)]["box_density"] / b0 if b0 > 0 else float("nan")
    for n in ("V1", "V2"):
        for sev in REPORT:
            stats[(n, sev)]["box_density_vs_V0"] = rel(n, sev)
    print("  box density, V / V0:   " + "   ".join(
        f"{sev}: V1 {rel('V1', sev):.2f} V2 {rel('V2', sev):.2f}" for sev in REPORT))
    print("  P3: V2 box density up at all three, most at SOG; polar densities down; NH/SH(SOG) > 1")
    print("  P4: falsified if V2/V0 box density within 0.80-1.20 at BOTH moderate and severe")
    return thr, stats


# ===========================================================================
# Stage 3 -- North Atlantic 1979
# ===========================================================================
def stage3(files, thr) -> dict:
    print("\n" + "=" * 78 + "\nSTAGE 3  North Atlantic 1979, Prosser's box, each version's own thresholds\n" + "=" * 78)
    num = {(n, s): 0.0 for n in VERSIONS for s in REPORT}
    den = 0.0
    t0 = time.time()
    geo = None
    for k, path in enumerate(files):
        ds = load_200(path)
        for day in per_day(ds):
            if geo is None:
                geo = E.grid_geometry(day["latitude"], day["longitude"], "map")
                lat = day["latitude"].values
                lon = day["longitude"].values
                lat_in = (lat >= PROSSER_BOX["lat"][0]) & (lat <= PROSSER_BOX["lat"][1])
                lon_in = (lon >= PROSSER_BOX["lon"][0]) & (lon <= PROSSER_BOX["lon"][1])
                w = np.cos(np.deg2rad(lat[lat_in]))
            out = ubf_versions(day["eastward_wind"], day["northward_wind"], day["geopotential"], geo)
            nt = day.sizes["time"]
            den += w.sum() * lon_in.sum() * nt
            for n in VERSIONS:
                x = out[n].transpose("latitude", "longitude", "time").values[lat_in][:, lon_in]
                for s in REPORT:
                    num[(n, s)] += float(((x >= thr[n][s]).reshape(x.shape[0], -1).sum(axis=1) * w).sum())
        print(f"  [{k + 1:>2}/{len(files)}] {Path(path).name}  {(time.time() - t0) / 60:5.1f} min")
    res = {}
    print(f"\n  {'ver':<4}" + "".join(f"{s + ' %':>16}" for s in REPORT)
          + "".join(f"{'ratio ' + s[:3]:>12}" for s in REPORT))
    for n in VERSIONS:
        lv = {s: 100 * num[(n, s)] / den for s in REPORT}
        rt = {s: lv[s] / PROSSER_1979[s] for s in REPORT}
        res[n] = dict(level=lv, ratio_vs_prosser=rt)
        print(f"  {n:<4}" + "".join(f"{lv[s]:>16.4f}" for s in REPORT) + "".join(f"{rt[s]:>12.3f}" for s in REPORT))
    print(f"  Prosser 1979  " + "  ".join(f"{s} {PROSSER_1979[s]} %" for s in REPORT))
    print("\n  implied FITTED level ratio = production fitted ratio x (V / V0, raw 1979):")
    for n in ("V1", "V2"):
        imp = {s: (PRODUCTION_RATIO[s] * res[n]["level"][s] / res["V0"]["level"][s]
                   if res["V0"]["level"][s] > 0 else float("nan")) for s in REPORT}
        res[n]["implied_fitted_ratio"] = imp
        print(f"  {n:<4}" + "   ".join(f"{s} {PRODUCTION_RATIO[s]:.2f} -> {imp[s]:.2f}" for s in REPORT))
    print("  (class rule band for a level ratio: 0.75-1.33)")
    return res


# ===========================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--global-glob", default="raw/global_sub48/*.grib")
    ap.add_argument("--na-glob", default="raw/north_atlantic/era5_na_1979-*.grib")
    ap.add_argument("--control-zarr", default="derived/global_audit/diagnostics_glob_2000-01.zarr")
    ap.add_argument("--out", default=str(REPO / "cat_outputs" / "ubf_stencil"))
    ap.add_argument("--stages", default="0123")
    args = ap.parse_args()

    base = Path(args.base)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    print(f">>> {VERSION}   script sha256 {file_sha(Path(__file__))}")
    gfiles = sorted(glob.glob(str(base / args.global_glob)))
    nfiles = sorted(glob.glob(str(base / args.na_glob)))
    print(f">>> global files {len(gfiles)}  ({args.global_glob})")
    print(f">>> NA files     {len(nfiles)}  ({args.na_glob})\n")
    summary = {"version": VERSION}

    if "0" in args.stages:
        summary["stage0"] = stage0()
    if "1" in args.stages:
        if not gfiles:
            print("!! no global files"); return 2
        summary["stage1"] = stage1(Path(gfiles[0]), base / args.control_zarr)
        if not summary["stage1"].get("ok"):
            (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
            return 2
    if "2" in args.stages:
        if len(gfiles) != 12:
            print(f"!! expected 12 global files, found {len(gfiles)}"); return 2
        thr, stats = stage2(gfiles)
        summary["thresholds"] = thr
        with open(out / "global_tail.csv", "w", newline="") as fh:
            rows = list(stats.values())
            keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[-1]).index(k) if k in rows[-1] else 0)
            wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader()
            for r in rows:
                wr.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
        if "3" in args.stages:
            if len(nfiles) != 12:
                print(f"!! expected 12 NA files, found {len(nfiles)}"); return 2
            summary["stage3"] = stage3(nfiles, thr)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\n>>> wrote {out}/summary.json and global_tail.csv   peak rss {peak_rss_gb():.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
