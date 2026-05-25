"""
make_synthetic_era5.py
======================
Build a synthetic ERA5-like trial dataset that mimics what
`data_download_final.py` would produce, so the pipeline can be
exercised end-to-end on this machine.

Domain   : 30–60 °N, 75 °W–0 °, 0.25° resolution
Time     : 8 steps × 3 h on 2016-07-01
Levels   : 175, 200, 225 hPa
Variables: u, v, t, z, d, vo, pv  (named with CDS/cfgrib short names)

The fields are constructed to be *physically plausible* — a midlatitude
upper-level jet with realistic vertical shear, a hydrostatic geopotential,
a smooth temperature field, and divergence/vorticity computed *consistent*
with the wind field so the diagnostics produce sensible magnitudes.

This is NOT real ERA5 data.  Its purpose is solely to demonstrate that the
pipeline runs and gives outputs in the right ballpark.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr


def build_synthetic_era5(out_path: Path) -> xr.Dataset:
    # ---- grid -------------------------------------------------------------
    lat = np.arange(60.0, 29.99, -0.25)      # 121 lats, descending  (ERA5 convention)
    lon = np.arange(-75.0,  0.001, 0.25)     # 301 lons
    levels = np.array([175.0, 200.0, 225.0])
    times = pd.date_range("2016-07-01T00:00", periods=8, freq="3h")

    LON, LAT = np.meshgrid(lon, lat)         # shape (nlat, nlon)

    # ---- jet stream profile ----------------------------------------------
    # Gaussian zonal jet centred at 42 °N with peak ~50 m/s near 200 hPa
    jet_lat0 = 42.0
    jet_width = 8.0
    u_jet_profile = np.exp(-((LAT - jet_lat0) / jet_width) ** 2)  # 0–1 envelope

    # Some longitudinal wiggle (Rossby-wave-like)
    wave = 0.15 * np.cos(np.deg2rad(LON * 2.0)) + 0.10 * np.sin(np.deg2rad(LON * 3.0 + 30))
    u_jet_profile *= 1.0 + wave

    # Vertical structure: jet strongest at 200 hPa, weaker above and below.
    level_factor = {175: 0.85, 200: 1.00, 225: 0.78}

    # Time-varying intensity to give the dataset some time variation
    t_factor = 1.0 + 0.05 * np.sin(2 * np.pi * np.arange(8) / 8)

    rng = np.random.default_rng(20160701)

    u_all = np.empty((len(lat), len(lon), len(times), len(levels)), dtype=np.float32)
    v_all = np.empty_like(u_all)
    for k, lev in enumerate(levels):
        for t_idx, _ in enumerate(times):
            base_u = 55.0 * u_jet_profile * level_factor[int(lev)] * t_factor[t_idx]
            base_v = 5.0 * np.sin(np.deg2rad(LON * 4.0)) * level_factor[int(lev)]
            # small spatial noise so derivatives don't vanish
            base_u += rng.normal(0, 1.5, size=base_u.shape)
            base_v += rng.normal(0, 1.0, size=base_v.shape)
            u_all[:, :, t_idx, k] = base_u
            v_all[:, :, t_idx, k] = base_v

    # ---- temperature ------------------------------------------------------
    # ICAO-like tropopause-region temperature: drops with latitude, slight
    # vertical structure (lapses through upper troposphere).
    T_base = 220.0 - 0.3 * (LAT - 30.0)                          # warmer south
    level_T = {175: -3.0, 200: 0.0, 225: 4.0}                    # K offsets
    t_all = np.empty_like(u_all)
    for k, lev in enumerate(levels):
        for t_idx in range(len(times)):
            t_all[:, :, t_idx, k] = (
                T_base + level_T[int(lev)] + rng.normal(0, 0.3, size=T_base.shape)
            )

    # ---- geopotential (hydrostatic, gives sensible delta_Z ~750 m) -------
    # Approx mid-tropospheric heights: 175 hPa ~ 12500 m, 200 hPa ~ 11800 m, 225 hPa ~ 11150 m
    h_at_lev = {175: 12500.0, 200: 11800.0, 225: 11150.0}
    g = 9.80665
    z_all = np.empty_like(u_all)
    for k, lev in enumerate(levels):
        H = h_at_lev[int(lev)]
        # slight latitudinal tilt (ridge southward) + small horizontal noise
        H_field = H - 80.0 * (LAT - 45.0)
        for t_idx in range(len(times)):
            z_all[:, :, t_idx, k] = (H_field + rng.normal(0, 5.0, size=H_field.shape)) * g

    # ---- divergence / vorticity from the wind field, geospatially --------
    # Compute on a sphere so the magnitudes are realistic.
    R = 6.371e6
    lat_r = np.deg2rad(lat)
    coslat = np.cos(lat_r)
    dlat = np.deg2rad(np.gradient(lat))                                       # (nlat,)
    dlon = np.deg2rad(np.gradient(lon))                                       # (nlon,)
    dy = R * dlat                                                             # (nlat,)
    dx = R * coslat[:, None] * dlon[None, :]                                  # (nlat, nlon)

    d_all = np.empty_like(u_all)
    vo_all = np.empty_like(u_all)
    for k in range(len(levels)):
        for t_idx in range(len(times)):
            u = u_all[:, :, t_idx, k]
            v = v_all[:, :, t_idx, k]
            du_dx = np.gradient(u, axis=1) / dx
            dv_dy = np.gradient(v, axis=0) / dy[:, None]
            dv_dx = np.gradient(v, axis=1) / dx
            du_dy = np.gradient(u, axis=0) / dy[:, None]
            d_all[:, :, t_idx, k]  = du_dx + dv_dy
            vo_all[:, :, t_idx, k] = dv_dx - du_dy

    # ---- potential vorticity ---------------------------------------------
    # Realistic upper-trop / lower-strat PV: 1–10 PVU.  We don't compute it
    # from theta gradients; we set a smooth field with the right magnitude
    # plus dependence on absolute vorticity.
    f = 2 * 7.2921e-5 * np.sin(lat_r)[:, None]                                # (nlat, 1)
    pv_all = np.empty_like(u_all)
    pv_level = {175: 5.5e-6, 200: 4.5e-6, 225: 3.0e-6}                        # K m^2 / kg / s (SI)
    for k, lev in enumerate(levels):
        for t_idx in range(len(times)):
            pv_all[:, :, t_idx, k] = pv_level[int(lev)] + 5e-7 * vo_all[:, :, t_idx, k] / max(np.median(np.abs(vo_all)), 1e-6)

    # ---- assemble ---------------------------------------------------------
    coords = {
        "latitude": ("latitude", lat),
        "longitude": ("longitude", lon),
        "time": ("time", times),
        "isobaricInhPa": ("isobaricInhPa", levels),
    }
    dims = ("latitude", "longitude", "time", "isobaricInhPa")

    ds = xr.Dataset(
        data_vars={
            "u":  (dims, u_all.astype(np.float32)),
            "v":  (dims, v_all.astype(np.float32)),
            "t":  (dims, t_all.astype(np.float32)),
            "z":  (dims, z_all.astype(np.float32)),
            "d":  (dims, d_all.astype(np.float32)),
            "vo": (dims, vo_all.astype(np.float32)),
            "pv": (dims, pv_all.astype(np.float32)),
        },
        coords=coords,
    )
    ds["u"].attrs = {"units": "m s**-1", "long_name": "U component of wind"}
    ds["v"].attrs = {"units": "m s**-1", "long_name": "V component of wind"}
    ds["t"].attrs = {"units": "K", "long_name": "Temperature"}
    ds["z"].attrs = {"units": "m**2 s**-2", "long_name": "Geopotential"}
    ds["d"].attrs = {"units": "s**-1", "long_name": "Divergence"}
    ds["vo"].attrs = {"units": "s**-1", "long_name": "Vorticity (relative)"}
    ds["pv"].attrs = {"units": "K m**2 kg**-1 s**-1", "long_name": "Potential vorticity"}
    ds.attrs["title"] = "Synthetic ERA5-like trial data, 2016-07-01 (NOT REAL DATA)"
    ds.attrs["disclaimer"] = "Generated for pipeline validation. Do not use for any science."

    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out_path)
    return ds


if __name__ == "__main__":
    out = Path("/home/claude/work/synthetic_era5_2016-07-01.nc")
    ds = build_synthetic_era5(out)
    print(f"Wrote {out}")
    print(ds)
