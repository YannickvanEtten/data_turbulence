"""
cat_month_render.py -- STEP 2 of 2 (any machine with matplotlib): animate a
frames file written by cat_month_frames.py.

Modes
  single       one map per severity (default: light, moderate, severe),
               colour = share of the 21 diagnostics above threshold (Ebar),
               with a strip underneath: box-mean Ebar over the month + cursor
  severities   light / moderate / severe side by side, same time step
  diagnostics  21 small maps, one per diagnostic (exceedance at the severity
               stored in the frames file, default moderate)
  stats        no animation: the "moving or popping?" figure + JSON

Examples
    python ada/cat_month_render.py frames/cat_frames_2020-02.npz
    python ada/cat_month_render.py frames/cat_frames_2020-02.npz --mode severities --days 7-14
    python ada/cat_month_render.py frames/cat_frames_2020-02.npz --mode diagnostics --days 8-10
    python ada/cat_month_render.py frames/cat_frames_2020-02.npz --mode stats
    add --jet to any animation to overlay the jet stream (dashed blue, 40/60/80 m/s)

Writes .gif always (Pillow) and .mp4 when ffmpeg is on the PATH.
Coastlines come from na_coastline.npz next to this file (GSHHS, no download).
"""
from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.ticker import FuncFormatter, PercentFormatter

HERE = Path(__file__).resolve().parent
COAST = HERE / "na_coastline.npz"

SEV_LABEL = {"light": "Light-or-greater", "light_to_moderate": "Light-to-moderate-or-greater",
             "moderate": "Moderate-or-greater", "moderate_to_severe": "Moderate-to-severe-or-greater",
             "severe": "Severe-or-greater"}
SEV_SHORT = {"light": "LOG", "light_to_moderate": "LMOG", "moderate": "MOG",
             "moderate_to_severe": "MSOG", "severe": "SOG"}
DIAG_LABEL = {"magnitude_pv": "|PV|", "colson_panofsky": "Colson-Panofsky", "brown1": "Brown Φ",
              "temperature_gradient": "|∇T|", "horizontal_divergence": "|div|",
              "vertical_wind_shear": "Vert. wind shear", "endlich": "Endlich", "deformation": "Deformation",
              "wind_speed": "Wind speed", "ngm2": "NGM2", "negative_richardson": "−Ri",
              "rva_magnitude": "|RVA|", "ubf": "UBF", "nva": "NVA", "brown2": "Brown ε",
              "vorticity_squared": "ζ²", "ti1": "TI1", "ngm1": "NGM1", "ti2": "TI2",
              "f2d": "Frontogenesis", "ncsu1": "NCSU1"}
CITIES = {"New York": (-74.0, 40.7), "Boston": (-71.1, 42.4), "Halifax": (-63.6, 44.6),
          "Gander": (-54.6, 48.9), "Bermuda": (-64.8, 32.3), "Azores": (-25.7, 37.7),
          "Lisbon": (-9.1, 38.7), "Madrid": (-3.7, 40.4), "Shannon": (-8.9, 52.7),
          "London": (-0.1, 51.5), "Glasgow": (-4.3, 55.9)}
PROSSER_BOX = (-55.0, -10.0, 36.0, 60.0)        # lon0, lon1, lat0, lat1 (Prosser 2023 Fig. 2)
JET_LEVELS = [40, 60, 80]                       # m/s
JET_STYLE = dict(colors="#1f5fa8", linestyles="--", alpha=0.8, zorder=2)


def draw_jet(ax, lon, lat, w, widths=(0.5, 0.8, 1.2), levels=JET_LEVELS):
    """Optional jet-stream contours: dashed blue, so they never read as coastlines."""
    return ax.contour(lon, lat, w, levels=levels, linewidths=list(widths)[:len(levels)], **JET_STYLE)

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "font.family": "DejaVu Sans"})


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def load(path: Path, days: str | None) -> dict:
    z = np.load(path, allow_pickle=False)
    d = {k: z[k] for k in z.files}
    d["time"] = d["time"].astype("datetime64[ns]")
    if days:
        a, b = (int(x) for x in days.split("-")) if "-" in days else (int(days), int(days))
        dom = (d["time"].astype("datetime64[D]") - d["time"].astype("datetime64[M]")).astype(int) + 1
        keep = (dom >= a) & (dom <= b)
        for k, v in list(d.items()):
            if isinstance(v, np.ndarray) and v.ndim >= 3 and v.shape[-3] == keep.size:
                d[k] = v[..., keep, :, :]
            elif k == "time":
                d[k] = v[keep]
    return d


def ebar(d: dict, sev: str) -> np.ndarray:
    n = d["n_defined"].astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, d[f"count_{sev}"] / np.maximum(n, 1), np.nan).astype(np.float32)


def single_diag(d: dict, sev_bits: str, i: int) -> np.ndarray:
    nx = int(d["bits_n_longitude"])
    return np.unpackbits(d[f"bits_{sev_bits}"][i], axis=-1, count=nx).astype(np.float32)


def box_mean(field: np.ndarray, lat, lon) -> np.ndarray:
    """cos(lat)-weighted mean over Prosser's box, per time step (as trends.py)."""
    lo0, lo1, la0, la1 = PROSSER_BOX
    jy = (lat >= la0) & (lat <= la1)
    jx = (lon >= lo0) & (lon <= lo1)
    f = field[:, jy][:, :, jx]
    w = np.cos(np.deg2rad(lat[jy]))[None, :, None] * np.ones_like(f)
    w = np.where(np.isfinite(f), w, 0)
    return np.nansum(np.nan_to_num(f) * w, axis=(1, 2)) / w.sum(axis=(1, 2))


def when(t: np.datetime64) -> str:
    dt = datetime.fromtimestamp(t.astype("datetime64[s]").astype(int), tz=timezone.utc)
    return dt.strftime("%a %d %b %Y  ·  %H:%M UTC")


# ---------------------------------------------------------------------------
# Map furniture
# ---------------------------------------------------------------------------
def lon_fmt(x, _):
    return "0°" if abs(x) < 1e-9 else (f"{abs(x):.0f}°W" if x < 0 else f"{x:.0f}°E")


def lat_fmt(y, _):
    return f"{y:.0f}°N" if y >= 0 else f"{-y:.0f}°S"


def base_map(ax, lat, lon, cities=True, ticks=True, box=True, fs=7):
    ax.set_xlim(lon.min(), lon.max())
    ax.set_ylim(lat.min(), lat.max())
    ax.set_aspect(1 / np.cos(np.deg2rad(45)))      # degrees look like distances at 45°N
    if COAST.exists():
        c = np.load(COAST)
        ax.plot(c["lon"], c["lat"], color="#3a3a3a", lw=0.6, zorder=3)
    if box:
        lo0, lo1, la0, la1 = PROSSER_BOX
        ax.plot([lo0, lo1, lo1, lo0, lo0], [la0, la0, la1, la1, la0], ls="--", lw=0.8,
                color="#2b6cb0", zorder=4)
        ax.text(lo1 - 0.5, la0 + 0.6, "analysis box", ha="right", va="bottom", fontsize=fs,
                color="#2b6cb0", zorder=4)
    if cities:
        for name, (x, y) in CITIES.items():
            ax.plot(x, y, "o", ms=2.5, color="#111", zorder=5)
            ax.annotate(name, (x, y), xytext=(3, 2), textcoords="offset points", fontsize=fs,
                        zorder=5, color="#111")
    if ticks:
        ax.set_xticks(np.arange(-75, 1, 15)); ax.set_yticks(np.arange(30, 61, 10))
        ax.xaxis.set_major_formatter(FuncFormatter(lon_fmt))
        ax.yaxis.set_major_formatter(FuncFormatter(lat_fmt))
        ax.grid(color="#999", lw=0.3, alpha=0.6)
    else:
        ax.set_xticks([]); ax.set_yticks([])
    ax.tick_params(length=2)


def cmap_share():
    cmap = plt.get_cmap("magma_r").copy()
    cmap.set_under("white"); cmap.set_bad("white")
    return cmap


def save(anim, fig, stem: Path, fps: int):
    anim.save(stem.with_suffix(".gif"), writer=animation.PillowWriter(fps=fps), dpi=100)
    made = [stem.with_suffix(".gif").name]
    if shutil.which("ffmpeg"):
        anim.save(stem.with_suffix(".mp4"), writer=animation.FFMpegWriter(fps=fps, bitrate=3000), dpi=150)
        made.append(stem.with_suffix(".mp4").name)
    plt.close(fig)
    print("written:", ", ".join(made))


# ---------------------------------------------------------------------------
# Animations
# ---------------------------------------------------------------------------
def anim_single(d, sev, stem, fps, jet=False):
    lat, lon, times = d["latitude"], d["longitude"], d["time"]
    F = ebar(d, sev)
    series = box_mean(F, lat, lon)
    wind = d.get("wind") if jet else None
    fig = plt.figure(figsize=(9.4, 5.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[6, 0.5], width_ratios=[40, 1], hspace=0.22, wspace=0.02,
                          left=0.08, right=0.92, top=0.885, bottom=0.06)
    ax = fig.add_subplot(gs[0, 0]); cax = fig.add_subplot(gs[0, 1]); axs = fig.add_subplot(gs[1, 0])
    mesh = ax.pcolormesh(lon, lat, F[0], cmap=cmap_share(), vmin=1 / 21, vmax=1, shading="nearest", zorder=1)
    base_map(ax, lat, lon)
    cb = fig.colorbar(mesh, cax=cax)
    cb.ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    cb.set_label("share of the 21 diagnostics above threshold")
    fig.suptitle(f"Clear-air turbulence over the North Atlantic at 200 hPa (≈ 11.8 km)  ·  "
                 f"{SEV_LABEL[sev]} ({SEV_SHORT[sev]})", fontsize=11, fontweight="bold", x=0.5, y=0.975)
    sub = "Colour: share of the 21 CAT diagnostics above their threshold  ·  ERA5, 0.25°, 3-hourly"
    if wind is not None:
        sub += "  ·  dashed blue: jet stream 40/60/80 m/s"
    fig.text(0.5, 0.925, sub, ha="center", fontsize=8, color="#444")
    stamp = ax.text(0.012, 0.975, "", transform=ax.transAxes, ha="left", va="top", fontsize=10,
                    fontweight="bold", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#888", lw=0.5, alpha=0.9))
    hours = (times - times[0]) / np.timedelta64(1, "D")
    axs.fill_between(hours, series, color="#c0504d", alpha=0.35, lw=0)
    axs.plot(hours, series, color="#8b1a1a", lw=0.9)
    cursor = axs.axvline(0, color="k", lw=1)
    axs.set_xlim(hours[0], hours[-1] if hours[-1] > 0 else 1)
    axs.set_ylim(0, max(series.max() * 1.15, 1e-4))
    span = hours[-1]
    step = 7 if span > 14 else (1 if span > 2 else 0.25)
    ticks = np.arange(0, span + 1e-9, step)
    axs.set_xticks(ticks)
    fmt = "%d %b" if step >= 1 else "%H:%M"
    axs.set_xticklabels([datetime.fromtimestamp((times[0] + np.timedelta64(int(x * 86400), "s"))
                         .astype("datetime64[s]").astype(int), tz=timezone.utc).strftime(fmt)
                         for x in ticks], fontsize=7)
    axs.set_yticks([])
    axs.text(0.0, 1.06, "analysis-box mean", transform=axs.transAxes, fontsize=7, color="#555",
             ha="left", va="bottom")
    axs.spines[["top", "right", "left"]].set_visible(False)
    fig.canvas.draw()                                   # line the strip up with the map
    pm, ps = ax.get_position(), axs.get_position()
    axs.set_position([pm.x0, ps.y0, pm.width, ps.height])
    pc = cax.get_position()
    cax.set_position([pm.x1 + 0.012, pm.y0, pc.width, pm.height])
    contours = []

    def draw(k):
        mesh.set_array(F[k].ravel())
        for c in contours:
            c.remove()
        contours.clear()
        if wind is not None:
            contours.append(draw_jet(ax, lon, lat, wind[k]))
        stamp.set_text(when(times[k]))
        cursor.set_xdata([hours[k], hours[k]])
        return []

    save(animation.FuncAnimation(fig, draw, frames=len(times)), fig, stem, fps)


def anim_severities(d, sevs, stem, fps, jet=False):
    lat, lon, times = d["latitude"], d["longitude"], d["time"]
    F = {s: ebar(d, s) for s in sevs}
    wind = d.get("wind") if jet else None
    fig, axes = plt.subplots(len(sevs), 1, figsize=(7.4, 3.35 * len(sevs) + 0.8))
    fig.subplots_adjust(left=0.1, right=0.86, top=0.90, bottom=0.04, hspace=0.18)
    meshes = []
    for ax, s in zip(axes, sevs):
        meshes.append(ax.pcolormesh(lon, lat, F[s][0], cmap=cmap_share(), vmin=1 / 21, vmax=1,
                                    shading="nearest", zorder=1))
        base_map(ax, lat, lon, cities=(s == sevs[0]), fs=6)
        ax.set_title(f"{SEV_LABEL[s]} ({SEV_SHORT[s]})", loc="left", fontsize=9)
    cax = fig.add_axes([0.875, 0.25, 0.02, 0.5])
    cb = fig.colorbar(meshes[0], cax=cax)
    cb.ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    cb.set_label("share of the 21 diagnostics above threshold")
    fig.suptitle("Clear-air turbulence at 200 hPa by severity", fontsize=11, fontweight="bold", y=0.985)
    stamp = fig.text(0.5, 0.935, "", ha="center", fontsize=10, fontweight="bold")
    contours = []

    def draw(k):
        for c in contours:
            c.remove()
        contours.clear()
        for ax, m, s in zip(axes, meshes, sevs):
            m.set_array(F[s][k].ravel())
            if wind is not None:
                contours.append(draw_jet(ax, lon, lat, wind[k], widths=(0.4, 0.7, 1.1)))
        stamp.set_text(when(times[k]) + ("   ·   dashed blue: jet 40/60/80 m/s" if wind is not None else ""))
        return []

    save(animation.FuncAnimation(fig, draw, frames=len(times)), fig, stem, fps)


def anim_diagnostics(d, stem, fps, jet=False):
    lat, lon, times = d["latitude"], d["longitude"], d["time"]
    sev = next(k[5:] for k in d if k.startswith("bits_") and k != "bits_n_longitude")
    names = [str(x) for x in d["diagnostics"]]
    E = [single_diag(d, sev, i) for i in range(len(names))]
    wind = d.get("wind") if jet else None
    nrow, ncol = 7, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 12.6))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.01, wspace=0.04, hspace=0.22)
    col = matplotlib.colors.ListedColormap(["white", "#b5323d"])
    meshes = []
    for ax, name, e in zip(axes.ravel(), names, E):
        meshes.append(ax.pcolormesh(lon, lat, e[0], cmap=col, vmin=0, vmax=1, shading="nearest", zorder=1))
        base_map(ax, lat, lon, cities=False, ticks=False, box=False)
        ax.set_title(DIAG_LABEL.get(name, name), fontsize=8, pad=2)
    for ax in axes.ravel()[len(names):]:
        ax.axis("off")
    fig.suptitle(f"The 21 CAT diagnostics one by one  ·  red = above the {SEV_LABEL[sev].lower()} "
                 f"({SEV_SHORT[sev]}) threshold  ·  200 hPa", fontsize=11, fontweight="bold", y=0.985)
    stamp = fig.text(0.5, 0.952, "", ha="center", fontsize=10)
    contours = []

    def draw(k):
        for c in contours:
            c.remove()
        contours.clear()
        for ax, m, e in zip(axes.ravel(), meshes, E):
            m.set_array(e[k].ravel())
            if wind is not None:
                contours.append(draw_jet(ax, lon, lat, wind[k], widths=(0.5,), levels=[60]))
        stamp.set_text(when(times[k]) + ("   ·   dashed blue: jet 60 m/s" if wind is not None else ""))
        return []

    save(animation.FuncAnimation(fig, draw, frames=len(times)), fig, stem, fps)


# ---------------------------------------------------------------------------
# Moving or popping?
# ---------------------------------------------------------------------------
def _pair(a, b, dy, dx):
    ny, nx = a.shape[-2:]
    ya, yb = (slice(0, ny - dy), slice(dy, ny)) if dy >= 0 else (slice(-dy, ny), slice(0, ny + dy))
    xa, xb = (slice(0, nx - dx), slice(dx, nx)) if dx >= 0 else (slice(-dx, nx), slice(0, nx + dx))
    return a[..., ya, xa], b[..., yb, xb]


def _cond(E, tau, dy, dx):
    a, b = _pair(E[: E.shape[0] - tau], E[tau:], dy, dx)
    n = a.sum()
    return float((a & b).sum() / n) if n else float("nan")


def stats(d, sev, consensus, stem, max_dx=48, max_dy=16, birth_radius=4):
    from scipy.ndimage import maximum_filter
    E = np.nan_to_num(ebar(d, sev)) >= consensus
    p = float(E.mean())
    res = {"severity": sev, "event": f"Ebar >= {consensus}", "p_event": p,
           "lift_east_west": {h: _cond(E, 0, 0, h) / p for h in range(1, max_dx + 1)},
           "lift_north_south": {h: _cond(E, 0, h, 0) / p for h in range(1, max_dy + 1)}}
    taus = [t for t in (1, 2, 4) if t < E.shape[0]]
    maps = {}
    for tau in taus:
        m = np.array([[_cond(E, tau, dy, dx) for dx in range(-max_dx, max_dx + 1)]
                      for dy in range(-max_dy, max_dy + 1)])
        i, j = np.unravel_index(np.nanargmax(m), m.shape)
        maps[tau] = m
        res[f"plus_{3 * tau}h"] = {"same_cell": m[max_dy, max_dx], "best_shift": m[i, j],
                                  "best_shift_deg_east": (j - max_dx) * 0.25,
                                  "best_shift_deg_north": -(i - max_dy) * 0.25}
    if 1 in maps:
        bx = int(round(res["plus_3h"]["best_shift_deg_east"] / 0.25))
        by = int(round(-res["plus_3h"]["best_shift_deg_north"] / 0.25))
        near = maximum_filter(E[:-1].astype(np.uint8), size=(1, 2 * birth_radius + 1, 2 * birth_radius + 1))
        src, dst = _pair(near, E[1:], by, bx)
        res["new_event_share_plus_3h"] = float((dst & (src == 0)).sum() / max(dst.sum(), 1))

    fig, axes = plt.subplots(1, len(taus) + 1, figsize=(4.6 * (len(taus) + 1), 3.9),
                             gridspec_kw={"width_ratios": [1.1] + [1] * len(taus), "wspace": 0.35})
    ax = axes[0]
    hx = np.array(list(res["lift_east_west"])) * 0.25
    hy = np.array(list(res["lift_north_south"])) * 0.25
    ax.plot(hx, list(res["lift_east_west"].values()), "-", lw=1.6, label="east–west")
    ax.plot(hy, list(res["lift_north_south"].values()), "-", lw=1.6, label="north–south")
    ax.axhline(1, color="grey", lw=0.8, ls=":")
    ax.set_yscale("log"); ax.set_xlabel("distance (degrees)")
    ax.set_ylabel("× more likely than average")
    ax.set_title("Neighbour of a hit, same time", loc="left"); ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    vmax = np.nanmax([np.nanmax(m) for m in maps.values()])
    ext = [(-max_dx - .5) * .25, (max_dx + .5) * .25, (-max_dy - .5) * .25, (max_dy + .5) * .25]
    for ax, tau in zip(axes[1:], taus):
        im = ax.imshow(maps[tau][::-1], origin="lower", extent=ext, vmin=0, vmax=vmax, cmap="viridis")
        r = res[f"plus_{3 * tau}h"]
        ax.plot(r["best_shift_deg_east"], r["best_shift_deg_north"], "r+", ms=12, mew=1.5)
        ax.plot(0, 0, "wx", ms=7)
        ax.set_title(f"{3 * tau} h later: same cell {r['same_cell']:.0%},\nbest shift {r['best_shift']:.0%}",
                     loc="left", fontsize=9)
        ax.set_xlabel("shift east (degrees)"); ax.set_ylabel("shift north (degrees)")
    cb = fig.colorbar(im, ax=axes[1:], fraction=0.025, pad=0.02)
    cb.set_label("P(hit at shifted cell | hit now)")
    t0 = np.datetime_as_string(d["time"][0], unit="D"); t1 = np.datetime_as_string(d["time"][-1], unit="D")
    fig.suptitle(f"Moving or popping?  {SEV_LABEL[sev]}, event = ≥ {consensus:.0%} of diagnostics agree,  "
                 f"{t0} to {t1}   ·   new events (no parent within 1°): "
                 f"{res.get('new_event_share_plus_3h', float('nan')):.0%}", fontsize=10, y=1.02)
    fig.savefig(stem.with_suffix(".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    stem.with_suffix(".json").write_text(json.dumps(res, indent=1, default=float))
    print("written:", stem.with_suffix(".png").name, stem.with_suffix(".json").name)


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("frames", help="cat_frames_YYYY-MM.npz from cat_month_frames.py")
    ap.add_argument("--mode", default="single", choices=["single", "severities", "diagnostics", "stats"])
    ap.add_argument("--severity", nargs="+", default=["light", "moderate", "severe"],
                    choices=list(SEV_LABEL))
    ap.add_argument("--days", default=None, help="day range inside the month, e.g. 7-14")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--jet", action="store_true", help="overlay jet-stream contours (dashed blue)")
    ap.add_argument("--consensus", type=float, default=0.25, help="stats mode: Ebar >= this is an event")
    ap.add_argument("--out", default=None, help="output folder (default: next to the frames file)")
    args = ap.parse_args()

    src = Path(args.frames)
    d = load(src, args.days)
    out = Path(args.out) if args.out else src.parent
    out.mkdir(parents=True, exist_ok=True)
    tag = src.stem.replace("cat_frames_", "") + (f"_d{args.days}" if args.days else "")
    print(f"{src.name}: {len(d['time'])} steps")

    if args.mode == "single":
        for s in args.severity:
            anim_single(d, s, out / f"cat_{tag}_{SEV_SHORT[s]}", args.fps, args.jet)
    elif args.mode == "severities":
        anim_severities(d, args.severity, out / f"cat_{tag}_severities", args.fps, args.jet)
    elif args.mode == "diagnostics":
        anim_diagnostics(d, out / f"cat_{tag}_21diagnostics", args.fps, args.jet)
    else:
        for s in args.severity:
            stats(d, s, args.consensus, out / f"cat_{tag}_{SEV_SHORT[s]}_moving_or_popping")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
