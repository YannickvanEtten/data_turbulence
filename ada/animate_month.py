"""
animate_month.py -- watch one month of North Atlantic CAT move (or pop).

Two outputs for one month of the 21 stored diagnostics:

1. AN ANIMATION (GIF, plus MP4 if ffmpeg exists), one frame per 3-hourly step.
   Colour = the ensemble exceedance fraction Ebar(c,t): the share of the
   diagnostics defined at cell c that exceed the chosen severity threshold
   (exactly the quantity trends.py averages). Black contours = 200 hPa wind
   speed (the jet), taken from the `wind_speed` diagnostic itself.

2. "MOVING OR POPPING?" STATISTICS for a binary consensus event
       E(c,t) = 1  if  Ebar(c,t) >= --consensus   (default 0.25: about 5 of the 21 agree)
   a) spatial lift   P(E(c+h,t) = 1 | E(c,t) = 1) / P(E = 1)   for zonal and meridional h
   b) lag-shift map  P(E(c+d,t+tau) = 1 | E(c,t) = 1)          over a grid of shifts d
      If events are ADVECTED, the maximum sits at an eastward shift and is much
      larger than at d = 0. If they POP up and die in place, the map is flat
      and the maximum stays near d = 0 and small.
   c) births         share of event cells at t+1 with no event cell at t within
      --birth-radius cells of the best (advected) position -> "popping" share.

Nothing is written to the derived/ or calibration/ folders. Reads one month only.

On ADA (from the repo root, in the pixi environment):
    pixi run python ada/animate_month.py --year 2020 --month 2
Locally on the one-day file in the repo (for a preview):
    python ada/animate_month.py --nc cat_outputs/diagnostics.nc \
        --thresholds calibration/thresholds_2026-09-11.json --out preview
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
DEFAULT_DERIVED = BASE / "derived" / "north_atlantic_audit"
DEFAULT_THRESHOLDS = BASE / "calibration" / "thresholds_2026-09-11.json"
DEFAULT_OUT = BASE / "explore" / "animations"

DIAGNOSTICS = [
    "magnitude_pv", "colson_panofsky", "brown1", "temperature_gradient",
    "horizontal_divergence", "vertical_wind_shear", "endlich", "deformation",
    "wind_speed", "ngm2", "negative_richardson", "rva_magnitude", "ubf", "nva",
    "brown2", "vorticity_squared", "ti1", "ngm1", "ti2", "f2d", "ncsu1",
]

# A few fixed points so the map can be read without coastlines (lon, lat).
LANDMARKS = {"NYC": (-74.0, 40.7), "Gander": (-54.6, 48.9), "Reykjavik": (-21.9, 64.1),
             "Shannon": (-8.9, 52.7), "London": (-0.1, 51.5), "Azores": (-27.0, 38.5),
             "Lisbon": (-9.1, 38.7)}


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def open_month(args) -> xr.Dataset:
    if args.nc:
        ds = xr.open_dataset(args.nc)
    else:
        store = Path(args.derived) / f"diagnostics_na_{args.year}-{args.month:02d}.zarr"
        if not store.exists():
            raise FileNotFoundError(store)
        ds = xr.open_zarr(store)
    return ds.transpose("time", "latitude", "longitude")


def ensemble_fraction(ds: xr.Dataset, thr: dict, severity: str):
    """Ebar(t, lat, lon): share of defined diagnostics at/above threshold. One diagnostic at a time."""
    shape = (ds.sizes["time"], ds.sizes["latitude"], ds.sizes["longitude"])
    n_exc = np.zeros(shape, np.uint8)
    n_def = np.zeros(shape, np.uint8)
    used = []
    for key in DIAGNOSTICS:
        if key not in ds:
            continue
        x = ds[key].values.astype(np.float32)
        ok = np.isfinite(x)
        n_def += ok
        n_exc += ok & (x >= thr[key][severity])
        used.append(key)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(n_def > 0, n_exc / np.maximum(n_def, 1), np.nan).astype(np.float32)
    return frac, used


# ---------------------------------------------------------------------------
# Moving or popping?
# ---------------------------------------------------------------------------
def shifted_pair(a: np.ndarray, b: np.ndarray, dy: int, dx: int):
    """Overlapping parts of a(y, x) and b(y+dy, x+dx) -- no wrap-around."""
    ny, nx = a.shape[-2:]
    ya, yb = (slice(0, ny - dy), slice(dy, ny)) if dy >= 0 else (slice(-dy, ny), slice(0, ny + dy))
    xa, xb = (slice(0, nx - dx), slice(dx, nx)) if dx >= 0 else (slice(-dx, nx), slice(0, nx + dx))
    return a[..., ya, xa], b[..., yb, xb]


def cond_prob(E: np.ndarray, tau: int, dy: int, dx: int) -> float:
    """P(E(t+tau, y+dy, x+dx) = 1 | E(t, y, x) = 1), pooled over all t, y, x."""
    a, b = E[: E.shape[0] - tau], E[tau:]
    a, b = shifted_pair(a, b, dy, dx)
    n = a.sum()
    return float((a & b).sum() / n) if n else float("nan")


def moving_or_popping(E: np.ndarray, max_dx: int, max_dy: int, taus, birth_radius: int):
    p = float(E.mean())
    out = {"p_event": p, "n_event_cells": int(E.sum())}
    # (a) spatial lift at lag 0. Grid north-to-south: +dy = south, +dx = east.
    out["lift_zonal"] = {h: cond_prob(E, 0, 0, h) / p for h in range(1, max_dx + 1)}
    out["lift_meridional"] = {h: cond_prob(E, 0, h, 0) / p for h in range(1, max_dy + 1)}
    # (b) lag-shift maps
    maps = {}
    for tau in taus:
        m = np.full((2 * max_dy + 1, 2 * max_dx + 1), np.nan)
        for i, dy in enumerate(range(-max_dy, max_dy + 1)):
            for j, dx in enumerate(range(-max_dx, max_dx + 1)):
                m[i, j] = cond_prob(E, tau, dy, dx)
        i, j = np.unravel_index(np.nanargmax(m), m.shape)
        maps[tau] = m
        out[f"tau{tau}"] = {"p_same_cell": m[max_dy, max_dx], "p_best_shift": m[i, j],
                            "best_shift_cells_(east,south)": (int(j - max_dx), int(i - max_dy))}
    # (c) births: event cells at t+1 with nothing within birth_radius of the advected origin at t
    if 1 in maps:
        bx, by = out["tau1"]["best_shift_cells_(east,south)"]
        from scipy.ndimage import maximum_filter
        near = maximum_filter(E[:-1].astype(np.uint8), size=(1, 2 * birth_radius + 1, 2 * birth_radius + 1))
        # an event at (t+1, y, x) "has a parent" if near[t, y-by, x-bx] is set
        src, dst = shifted_pair(near, E[1:], by, bx)
        out["birth_share_tau1"] = float((dst & (src == 0)).sum() / max(dst.sum(), 1))
    return out, maps


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _add_landmarks(ax):
    for name, (lo, la) in LANDMARKS.items():
        if -75 <= lo <= 0 and 30 <= la <= 60:
            ax.plot(lo, la, "k^", ms=3)
            ax.annotate(name, (lo, la), xytext=(3, 3), textcoords="offset points", fontsize=7)


def make_animation(frac, wind, lat, lon, times, title, out: Path, fps: int, use_cartopy: bool):
    fig = plt.figure(figsize=(11, 5.2))
    if use_cartopy:
        import cartopy.crs as ccrs
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.coastlines(linewidth=0.6)
        kw = {"transform": ccrs.PlateCarree()}
    else:
        ax = fig.add_subplot(1, 1, 1)
        kw = {}
        _add_landmarks(ax)
    cmap = plt.get_cmap("magma_r").copy()
    cmap.set_under("white")
    mesh = ax.pcolormesh(lon, lat, frac[0], cmap=cmap, vmin=1 / 21, vmax=1.0, shading="auto", **kw)
    cb = fig.colorbar(mesh, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("share of diagnostics exceeding threshold")
    ax.set_xlim(lon.min(), lon.max()); ax.set_ylim(lat.min(), lat.max())
    ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    ttl = ax.set_title("")
    contours = []

    def draw(k):
        mesh.set_array(frac[k].ravel())
        for c in contours:
            c.remove()
        contours.clear()
        if wind is not None:
            contours.append(ax.contour(lon, lat, wind[k], levels=[40, 60, 80],
                                       colors="k", linewidths=[0.5, 0.9, 1.4], **kw))
        ttl.set_text(f"{title}   {np.datetime_as_string(times[k], unit='h')} UTC")
        return [mesh, ttl]

    anim = animation.FuncAnimation(fig, draw, frames=len(times), blit=False)
    anim.save(out.with_suffix(".gif"), writer=animation.PillowWriter(fps=fps), dpi=90)
    if shutil.which("ffmpeg"):
        anim.save(out.with_suffix(".mp4"), writer=animation.FFMpegWriter(fps=fps, bitrate=2400), dpi=120)
    plt.close(fig)


def plot_moving_or_popping(stats, maps, max_dx, max_dy, title, out: Path):
    taus = sorted(maps)
    fig, axes = plt.subplots(1, len(taus) + 1, figsize=(4.2 * (len(taus) + 1), 3.8))
    ax = axes[0]
    hz = list(stats["lift_zonal"]); hm = list(stats["lift_meridional"])
    ax.plot(hz, list(stats["lift_zonal"].values()), "o-", label="zonal (E-W)")
    ax.plot(hm, list(stats["lift_meridional"].values()), "s-", label="meridional (N-S)")
    ax.axhline(1, color="grey", lw=0.8)
    ax.set_xlabel("distance (grid cells, 0.25°)"); ax.set_ylabel("lift = P(neighbour | hit) / P(hit)")
    ax.set_title("same time"); ax.legend(fontsize=8); ax.set_yscale("log")
    vmax = np.nanmax([np.nanmax(m) for m in maps.values()])
    for ax, tau in zip(axes[1:], taus):
        im = ax.imshow(maps[tau], origin="upper", vmin=0, vmax=vmax, cmap="viridis",
                       extent=[-max_dx - .5, max_dx + .5, max_dy + .5, -max_dy - .5])
        bx, by = stats[f"tau{tau}"]["best_shift_cells_(east,south)"]
        ax.plot(bx, by, "r+", ms=12); ax.plot(0, 0, "wx", ms=8)
        ax.set_title(f"+{3 * tau} h: P(hit at shift | hit now)")
        ax.set_xlabel("shift east (cells)"); ax.set_ylabel("shift south (cells)")
    fig.colorbar(im, ax=axes[1:], fraction=0.02)
    fig.suptitle(title, fontsize=10)
    fig.savefig(out.with_suffix(".png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, default=2020)
    ap.add_argument("--month", type=int, default=2)
    ap.add_argument("--severity", default="moderate",
                    choices=["light", "light_to_moderate", "moderate", "moderate_to_severe", "severe"])
    ap.add_argument("--consensus", type=float, default=0.25,
                    help="Ebar >= this counts as an event (0.25 ~ 5 of 21 agree; 0.5 is very rare)")
    ap.add_argument("--derived", default=str(DEFAULT_DERIVED))
    ap.add_argument("--thresholds", default=str(DEFAULT_THRESHOLDS))
    ap.add_argument("--nc", default=None, help="read a single netCDF file instead of a zarr month")
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--fps", type=int, default=6)
    ap.add_argument("--max-dx", type=int, default=48, help="search window east-west (cells); 12 h needs ~40")
    ap.add_argument("--max-dy", type=int, default=16)
    ap.add_argument("--birth-radius", type=int, default=4)
    args = ap.parse_args()

    thr = json.loads(Path(args.thresholds).read_text())["thresholds"]
    ds = open_month(args)
    lat, lon, times = ds.latitude.values, ds.longitude.values, ds.time.values
    frac, used = ensemble_fraction(ds, thr, args.severity)
    wind = ds["wind_speed"].values if "wind_speed" in ds else None
    tag = np.datetime_as_string(times[0], unit="M")
    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / f"cat_{tag}_{args.severity}"
    print(f"{tag}: {len(times)} steps, {len(used)} diagnostics, severity {args.severity}")

    try:
        import cartopy  # noqa: F401
        use_cartopy = True
    except ImportError:
        use_cartopy = False

    E = np.nan_to_num(frac, nan=0.0) >= args.consensus
    taus = [t for t in (1, 2, 4) if t < len(times)]
    stats, maps = moving_or_popping(E, args.max_dx, args.max_dy, taus, args.birth_radius)
    plot_moving_or_popping(stats, maps, args.max_dx, args.max_dy,
                           f"{tag}, {args.severity}, event = Ebar ≥ {args.consensus}", Path(f"{stem}_moving_or_popping"))
    with open(f"{stem}_moving_or_popping.json", "w") as f:
        json.dump(stats, f, indent=1, default=float)
    for k, v in stats.items():
        print(f"  {k}: {v}")

    make_animation(frac, wind, lat, lon, times, f"CAT ({args.severity}), North Atlantic 200 hPa",
                   stem, args.fps, use_cartopy)
    print(f"written: {stem}.gif (+ .mp4 if ffmpeg), {stem}_moving_or_popping.png/.json")


if __name__ == "__main__":
    main()
