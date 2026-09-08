#!/usr/bin/env python
"""
Prosser (2023) figure replication — the RENDER half.

This script computes nothing from ERA5. It reads the small CSV reductions that
`ada/per_diagnostic_trend.py` already writes and turns them into the figure
panels of Prosser et al. (2023), so that a plot can be iterated on in seconds
instead of in a queue.

    reduce (ADA, hours, jobs/16 and jobs/21)  ->  cat_outputs/*.csv  (~100 KB)
    render (anywhere, seconds, this script)   ->  cat_outputs/figures/*.png

That split is deliberate. The CSVs are the archived result; the figures are
disposable and must be reproducible from one command with no notebook state.

WHAT IT REPRODUCES
------------------
  Figure 3(a)   annual diagnostic-mean MOG CAT probability for the North
                Atlantic box, 42 points, with the fitted trend.
                Prosser's panel (b) is the USA box and needs a download this
                project has not done — see PLAN_prosser_replication.md.
  Figure 4      the same fit for each of the 21 diagnostics separately.
  Figure S3     Figure 3 at LOG and at SOG.
  Figure S4     Figure 4 at LOG and at SOG.

Prosser's conventions, followed deliberately so the panels can be laid side by
side with his: blue crosses for the 42 yearly values, two red crosses for the
fitted 1979 and 2020 values, a green trend line that is SOLID when the trend is
significant at p = 0.05 (two-sided) and DASHED when it is not, and the fitted
change stated at the top of each panel.

INPUT
-----
`per_diagnostic_trend.py --csv cat_outputs/per_diagnostic_<season>_<severity>.csv`
writes two files. This script reads both:

  ..._series.csv   year, <21 diagnostic columns>, ensemble      <- the points
  ....csv          diagnostic, level_1979, level_2020, change,  <- the fit
                   ci_lo, ci_hi, t, r2, ...

The fit is RECOMPUTED here from the series rather than read, so the figure is
self-contained. The fit CSV is then used as a cross-check: any disagreement
above 0.5 percentage points is reported loudly, because it would mean the
picture and the table in STATUS disagree about the same run.

USAGE
-----
    python ada/prosser_figures.py                        # everything it finds
    python ada/prosser_figures.py --season annual --severity moderate
    python ada/prosser_figures.py --outputs cat_outputs --figdir /tmp/figs

Runs on ADA or on Windows — it needs only numpy and matplotlib, not the
pipeline. On ADA the pixi environment does not currently declare matplotlib;
either add `matplotlib-base = "*"` to pixi.toml, or (simpler) commit the CSVs
and render locally, which is the intended workflow.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUTS = REPO / "cat_outputs"

SEVERITY_LABEL = {
    "light": "LOG",
    "light_to_moderate": "LMOG",
    "moderate": "MOG",
    "moderate_to_severe": "MSOG",
    "severe": "SOG",
}
SEASON_TITLE = {
    "djf": "DJF", "mam": "MAM", "jja": "JJA", "son": "SON", "annual": "Annual",
}

# Which main-paper / SI figure each (season, severity) pair corresponds to.
# Prosser's Figures 3 and 4 are ANNUAL and MOG; the LOG and SOG versions are
# S3 and S4. Seasonal cuts have no published counterpart — they are ours, and
# are labelled as such rather than passed off as a replication.
FIGURE_ID = {
    ("annual", "moderate"): ("Figure 3a", "Figure 4"),
    ("annual", "light"): ("Figure S3-LOG a", "Figure S4-LOG"),
    ("annual", "severe"): ("Figure S3-SOG a", "Figure S4-SOG"),
}

# Prosser (2023) Table 1, annual row: fitted relative increase 1979->2020,
# diagnostic-mean, North Atlantic box. Used only to annotate the ensemble
# panel with the published value. CALIBRATION_REFERENCE.md §4.4.
PROSSER_ANNUAL_REL = {"light": 0.17, "light_to_moderate": 0.28,
                      "moderate": 0.37, "moderate_to_severe": 0.46,
                      "severe": 0.55}
# and the published 1979 level in hours per year, same source.
PROSSER_ANNUAL_HOURS_1979 = {"light": 466.5, "light_to_moderate": 152.7,
                             "moderate": 70.0, "moderate_to_severe": 35.5,
                             "severe": 17.7}

BLUE = "#1f4e9c"
RED = "#c1272d"
GREEN = "#2e7d32"


# --------------------------------------------------------------------------
# statistics — kept identical in form to ada/full_trend_check.py's ols()
# --------------------------------------------------------------------------

def tcrit(df: int) -> float:
    """Two-sided 5% critical value. scipy if present, else a table."""
    try:
        from scipy import stats  # noqa: PLC0415
        return float(stats.t.ppf(0.975, df))
    except Exception:
        table = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
                 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 12: 2.179,
                 15: 2.131, 20: 2.086, 25: 2.060, 30: 2.042, 40: 2.021,
                 50: 2.009, 60: 2.000, 80: 1.990, 100: 1.984, 1000: 1.962}
        if df in table:
            return table[df]
        keys = sorted(table)
        if df < keys[0]:
            return table[keys[0]]
        if df > keys[-1]:
            return 1.960
        hi = next(k for k in keys if k > df)
        lo = max(k for k in keys if k < df)
        f = (df - lo) / (hi - lo)
        return table[lo] + f * (table[hi] - table[lo])


def ols(x: np.ndarray, y: np.ndarray) -> dict:
    """Least squares with the slope's standard error and 95% half-width."""
    n = len(x)
    df = n - 2
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    resid = y - pred
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    sxx = float(np.sum((x - x.mean()) ** 2))
    se = math.sqrt(ss_res / df / sxx) if df > 0 and sxx > 0 else float("nan")
    t = slope / se if se and np.isfinite(se) and se > 0 else float("nan")
    half = tcrit(df) * se if np.isfinite(se) else float("nan")
    return dict(slope=float(slope), intercept=float(intercept), r2=r2,
                se=se, t=float(t), half=float(half), n=n, df=df,
                tcrit=tcrit(df))


def fit_summary(years: np.ndarray, rate: np.ndarray) -> dict:
    """Fitted endpoints, relative change and its interval — as STATUS uses."""
    f = ols(years, rate)
    y0 = f["slope"] * years[0] + f["intercept"]
    y1 = f["slope"] * years[-1] + f["intercept"]
    span = float(years[-1] - years[0])
    rel = (y1 - y0) / y0 if y0 > 0 else float("nan")
    rel_half = span * f["half"] / y0 if y0 > 0 else float("nan")
    return dict(fit=f, y0=y0, y1=y1, rel=rel,
                lo=rel - rel_half, hi=rel + rel_half,
                significant=abs(f["t"]) >= f["tcrit"])


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------

def read_series(path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Read a *_series.csv into (years, {column: values})."""
    path = Path(path)
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path} has a header but no rows")
    cols = [c for c in rows[0] if c and c != "year"]
    years = np.array([int(r["year"]) for r in rows], float)
    data = {c: np.array([float(r[c]) for r in rows], float) for c in cols}
    if "ensemble" not in data:
        raise SystemExit(f"{path} has no 'ensemble' column — is this a "
                         f"per_diagnostic_*_series.csv?")
    return years, data


def read_fit(path) -> dict[str, dict]:
    """Read the companion fit CSV, if it is there. Returns {} if not."""
    path = Path(path)
    if not path.exists():
        return {}
    with path.open() as fh:
        return {r["diagnostic"]: r for r in csv.DictReader(fh)}


def discover(outputs: Path) -> list[tuple[str, str, Path]]:
    """Find every per_diagnostic_<season>_<severity>_series.csv present."""
    found = []
    for p in sorted(outputs.glob("per_diagnostic_*_series.csv")):
        stem = p.name[len("per_diagnostic_"):-len("_series.csv")]
        for sev in sorted(SEVERITY_LABEL, key=len, reverse=True):
            if stem.endswith("_" + sev):
                season = stem[: -len(sev) - 1]
                if season in SEASON_TITLE:
                    found.append((season, sev, p))
                break
    return found


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------

def _panel(ax, years, rate, s, title, *, show_prosser=None, compact=False):
    """One Prosser-style regression panel."""
    ax.plot(years, rate * 100.0, "+", color=BLUE, markersize=5 if compact else 7,
            markeredgewidth=1.1, zorder=3, label="annual value")

    xs = np.array([years[0], years[-1]])
    ys = np.array([s["y0"], s["y1"]]) * 100.0
    ax.plot(xs, ys, "-" if s["significant"] else "--", color=GREEN,
            linewidth=1.6, zorder=2,
            label="fit (p<0.05)" if s["significant"] else "fit (n.s.)")
    ax.plot(xs, ys, "+", color=RED, markersize=8 if compact else 11,
            markeredgewidth=1.8, zorder=4, label="fitted endpoints")

    txt = f"{s['rel']:+.0%}"
    if show_prosser is not None:
        txt += f"   (Prosser {show_prosser:+.0%})"
    if not s["significant"]:
        txt += "  n.s."
    ax.set_title(f"{title}\n{txt}", fontsize=8 if compact else 11,
                 linespacing=1.35)
    ax.tick_params(labelsize=7 if compact else 9)
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.margins(x=0.03)


def figure_ensemble(years, rate, season, severity, figdir, fig_id):
    import matplotlib.pyplot as plt

    s = fit_summary(years, rate)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    show = (PROSSER_ANNUAL_REL.get(severity) if season == "annual" else None)
    _panel(ax, years, rate, s,
           f"{SEASON_TITLE[season]} {SEVERITY_LABEL[severity]} CAT probability, "
           f"North Atlantic box\n(36–60°N, 55–10°W), diagnostic mean over 21",
           show_prosser=show)
    ax.set_xlabel("year")
    ax.set_ylabel(f"{SEVERITY_LABEL[severity]} CAT probability (%)")
    ax.legend(fontsize=8, framealpha=0.9)
    sub = (f"ours {s['y0']:.3%} → {s['y1']:.3%}   "
           f"95% CI on the change [{s['lo']:+.0%}, {s['hi']:+.0%}]   "
           f"t={s['fit']['t']:.2f}, R²={s['fit']['r2']:.2f}, n={s['fit']['n']}")
    fig.text(0.5, 0.005, sub, ha="center", fontsize=8, color="#444444")
    fig.suptitle(f"replicating {fig_id}", fontsize=9, color="#666666", y=0.995)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))

    out = figdir / f"fig3_{season}_{severity}_ensemble.png"
    fig.savefig(out, dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    return out, s


def figure_per_diagnostic(years, data, season, severity, figdir, fig_id):
    import matplotlib.pyplot as plt

    names = [c for c in data if c != "ensemble"]
    n = len(names)
    ncol = 5 if n > 16 else 4
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.5 * nrow),
                             squeeze=False)

    summaries = {}
    order = sorted(names, key=lambda k: -float(np.mean(data[k])))
    for i, name in enumerate(order):
        ax = axes[i // ncol][i % ncol]
        s = fit_summary(years, data[name])
        summaries[name] = s
        _panel(ax, years, data[name], s, name, compact=True)
    for j in range(len(order), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    n_sig = sum(1 for s in summaries.values() if s["significant"])
    fig.suptitle(
        f"replicating {fig_id} — {SEASON_TITLE[season]} "
        f"{SEVERITY_LABEL[severity]} CAT probability by diagnostic, "
        f"North Atlantic box, {int(years[0])}–{int(years[-1])}\n"
        f"{n_sig} of {n} trends significant at p=0.05 (solid line); "
        f"panels ordered by mean level",
        fontsize=11, y=0.998)
    fig.supxlabel("year", fontsize=10)
    fig.supylabel(f"{SEVERITY_LABEL[severity]} CAT probability (%)", fontsize=10)
    fig.tight_layout(rect=(0.012, 0.012, 1, 0.965))

    out = figdir / f"fig4_{season}_{severity}_per_diagnostic.png"
    fig.savefig(out, dpi=170)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    return out, summaries


# --------------------------------------------------------------------------
# the per-diagnostic comparison itself — shared by the script and the notebook
# so the two can never disagree about the same run
# --------------------------------------------------------------------------

def compare_to_prosser(years, data, published, stencil, ylims=None):
    """One row per diagnostic: ours, his, the ratios, and transparent flags.

    NO SCORE AND NO VERDICT. §12.8's lesson, applied to this project's own
    tooling: report the statistic and let a human read it. The flags below are
    plain thresholds, each one checkable by eye against the columns beside it.

      LEVEL  our 1979 exceedance level is more than 2x from his either way.
             This is the one that carries information about the FORMULA.
             Thresholds are percentiles of our own data, so a constant scale
             error - units, resolution, any uniform bias - cancels exactly out
             of the exceedance field (STATUS §5 point 5) and CANNOT move this.
      TREND  the fitted relative change is more than 2x from his, counted only
             where his own trend is big enough for a ratio to mean anything.
      SIG    we and he disagree on significance at p = 0.05.
      AXIS   our data does not fit on his published y axis, materially: the
             fitted line leaves it, or more than a tenth of the points do.
             His limits are read off a figure, so smaller overhangs are noise.
    """
    rows = []
    for name, (panel, title, rel, abs_, p, lvl) in published.items():
        f = fit_summary(years, data[name])
        our_rel = f["rel"] * 100.0
        our_lvl = f["y0"] * 100.0
        our_p = pvalue(f["fit"]["t"], f["fit"]["df"])
        level_ratio = our_lvl / lvl if lvl else float("nan")
        trend_ratio = our_rel / rel if rel else float("nan")

        flags = []
        if not (0.5 <= level_ratio <= 2.0):
            flags.append("LEVEL")
        if abs(rel) > 5 and not (0.5 <= trend_ratio <= 2.0):
            flags.append("TREND")
        if f["significant"] != (p < 0.05):
            flags.append("SIG")
        if ylims is not None:
            vals = data[name] * 100.0
            lo, hi = ylims[name]
            n_out = int(((vals < lo) | (vals > hi)).sum())
            fit_out = not (lo <= f["y0"] * 100 <= hi
                           and lo <= f["y1"] * 100 <= hi)
            if fit_out or n_out > 0.1 * len(vals):
                flags.append("AXIS")

        rows.append(dict(
            diagnostic=name, panel=panel, prosser_title=title,
            uses_vertical_stencil=stencil[name],
            ours_change_pct=our_rel, prosser_change_pct=rel,
            trend_ratio=trend_ratio,
            ours_ci_lo_pct=f["lo"] * 100, ours_ci_hi_pct=f["hi"] * 100,
            ours_level_1979_pct=our_lvl, prosser_level_1979_pct=lvl,
            level_ratio=level_ratio,
            ours_abs_pct_per_yr=f["fit"]["slope"] * 100,
            prosser_abs_pct_per_yr=abs_,
            ours_t=f["fit"]["t"], ours_p=our_p, prosser_p=p,
            ours_significant=f["significant"], prosser_significant=p < 0.05,
            flags=" ".join(flags), n_flags=len(flags)))
    rows.sort(key=lambda r: (-r["n_flags"], r["panel"]))
    return rows


# --------------------------------------------------------------------------
# side-by-side layout: his order, his titles, his axes, both sets of numbers
# --------------------------------------------------------------------------

def pvalue(t: float, df: int) -> float:
    """Two-sided p for a t statistic. scipy if present, else a normal approx."""
    try:
        from scipy import stats  # noqa: PLC0415
        return float(2 * stats.t.sf(abs(t), df))
    except Exception:
        z = abs(t)
        return float(math.erfc(z / math.sqrt(2)))


def _fmt_p(p: float) -> str:
    if p >= 0.01:
        return f"{p:.2g}"
    return f"{p:.0e}".replace("e-0", "e-")


def figure_prosser_layout(years, data, figdir, published, ylims,
                          match_axes: bool = True):
    """Our 21 panels in Prosser's Figure 4 layout, for direct comparison.

    His figure fills 7 rows x 3 columns COLUMN-MAJOR: (a)-(g) down the first
    column, (h)-(n) the second, (o)-(u) the third. Same order here, same
    titles, and — when match_axes — the same y limits, so the two figures can
    be laid side by side and read panel for panel.

    Under each panel, two lines: his printed rel/abs/p, then ours computed the
    same way. `abs` is the slope in percentage points per year, which is what
    his caption says his is.

    A panel whose data will not fit on his axis is flagged in red with our
    actual range. That is a finding, not a plotting failure: it means the
    LEVEL disagrees, and where the whole panel is off scale it belongs on the
    investigate list.
    """
    import matplotlib.pyplot as plt

    order = sorted(published, key=lambda k: published[k][0])   # a..u
    ncol, nrow = 3, 7
    fig, axes = plt.subplots(nrow, ncol, figsize=(16.5, 26), squeeze=False)

    offscale = {}
    for i, name in enumerate(order):
        ax = axes[i % nrow][i // nrow]                          # column-major
        panel, title, rel, abs_, p, lvl = published[name]
        s = fit_summary(years, data[name])
        _panel(ax, years, data[name], s, "", compact=True)
        ax.set_title(f"{title}\n({panel})", fontsize=9, linespacing=1.3)

        vals = data[name] * 100.0
        lo, hi = ylims[name]
        n_out = int(((vals < lo) | (vals > hi)).sum())
        fit_out = not (lo <= s["y0"] * 100 <= hi and lo <= s["y1"] * 100 <= hi)
        if match_axes:
            ax.set_ylim(lo, hi)
        if n_out:
            offscale[name] = dict(n_out=n_out, fit_out=fit_out,
                                  ours=(float(vals.min()), float(vals.max())),
                                  his=(lo, hi))
            # His y limits were read off a rendered figure, so a point or two
            # over the edge means nothing. Only flag a MATERIAL disagreement:
            # the fitted line off scale, or more than a tenth of the points.
            material = fit_out or n_out > 0.1 * len(vals)
            if match_axes and material:
                severity = ("ENTIRELY OFF SCALE" if n_out == len(vals)
                            else f"{n_out}/{len(vals)} points off scale")
                ax.text(0.5, 0.5, f"{severity}\nours {vals.min():.3f}–{vals.max():.3f}",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=8.5, color=RED, weight="bold",
                        bbox=dict(fc="white", ec=RED, alpha=0.85, lw=0.8))

        our_abs = s["fit"]["slope"] * 100.0
        our_p = pvalue(s["fit"]["t"], s["fit"]["df"])
        ax.text(0.0, -0.30,
                f"his   rel={rel:+.1f}   abs={abs_:+.3f}   p={_fmt_p(p)}",
                transform=ax.transAxes, fontsize=8.5, family="monospace",
                color="#333333")
        ax.text(0.0, -0.42,
                f"ours  rel={s['rel'] * 100:+.1f}   abs={our_abs:+.3f}   "
                f"p={_fmt_p(our_p)}",
                transform=ax.transAxes, fontsize=8.5, family="monospace",
                color=BLUE if abs(s["rel"] * 100 - rel) < 10 else RED)

    axis_note = ("drawn on HIS y axes" if match_axes else "drawn on OUR y axes")
    fig.suptitle(
        "Prosser (2023) Figure 4 — his layout, his titles, "
        f"{axis_note}\n"
        "blue crosses: our 42 annual values · red crosses: our fitted 1979 and "
        "2020 · solid green: significant at p=0.05\n"
        "under each panel: his printed rel/abs/p, then ours "
        "(rel = % change 1979→2020, abs = slope in %/yr)",
        fontsize=12, y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.972), h_pad=5.2)

    suffix = "his_axes" if match_axes else "our_axes"
    out = figdir / f"fig4_prosser_layout_{suffix}.png"
    fig.savefig(out, dpi=170)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    return out, offscale


# --------------------------------------------------------------------------

def crosscheck(summaries: dict[str, dict], fit_csv: dict[str, dict],
               tol_pp: float = 0.5) -> list[str]:
    """Recomputed fit vs the archived fit CSV. Silence means they agree."""
    problems = []
    for name, s in summaries.items():
        row = fit_csv.get(name)
        if not row:
            continue
        try:
            archived = float(row["change"])
        except (KeyError, ValueError):
            continue
        d = abs(archived - s["rel"]) * 100.0
        if d > tol_pp:
            problems.append(
                f"{name}: figure says {s['rel']:+.1%}, "
                f"CSV says {archived:+.1%}  (Δ {d:.2f} pp)")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Render Prosser (2023) Figures 3/4/S3/S4 from the "
                    "archived CSV reductions.")
    ap.add_argument("--outputs", default=str(DEFAULT_OUTPUTS),
                    help="directory holding per_diagnostic_*_series.csv")
    ap.add_argument("--figdir", default=None,
                    help="where to write figures (default <outputs>/figures)")
    ap.add_argument("--season", default=None, choices=sorted(SEASON_TITLE))
    ap.add_argument("--severity", default=None, choices=sorted(SEVERITY_LABEL))
    args = ap.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        print("matplotlib is not installed in this environment.\n"
              "  On ADA: add  matplotlib-base = \"*\"  to pixi.toml, or\n"
              "  commit cat_outputs/*.csv, pull on Windows, and run this there\n"
              "  (the intended workflow — the CSVs are the archived result).",
              file=sys.stderr)
        return 2

    outputs = Path(args.outputs)
    figdir = Path(args.figdir) if args.figdir else outputs / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    available = discover(outputs)
    if args.season:
        available = [a for a in available if a[0] == args.season]
    if args.severity:
        available = [a for a in available if a[1] == args.severity]
    if not available:
        print(f"no per_diagnostic_*_series.csv found under {outputs}\n"
              f"  produce them with jobs/21, e.g.\n"
              f"    pixi run python ada/per_diagnostic_trend.py "
              f"--season annual --severity moderate \\\n"
              f"        --csv cat_outputs/per_diagnostic_annual_moderate.csv",
              file=sys.stderr)
        return 1

    print(f"outputs : {outputs}")
    print(f"figures : {figdir}\n")

    all_problems = []
    for season, severity, series_path in available:
        years, data = read_series(series_path)
        fit_csv = read_fit(series_path.with_name(
            series_path.name.replace("_series.csv", ".csv")))
        ens_id, diag_id = FIGURE_ID.get(
            (season, severity),
            (f"[no published counterpart: {SEASON_TITLE[season]} "
             f"{SEVERITY_LABEL[severity]} ensemble]",
             f"[no published counterpart: {SEASON_TITLE[season]} "
             f"{SEVERITY_LABEL[severity]} per-diagnostic]"))

        print(f"=== {SEASON_TITLE[season]} {SEVERITY_LABEL[severity]} "
              f"({len(years)} years, {len(data) - 1} diagnostics)")

        p1, s_ens = figure_ensemble(years, data["ensemble"], season, severity,
                                    figdir, ens_id)
        print(f"    {p1.name:<48} {ens_id}")
        print(f"      ours {s_ens['rel']:+.0%} "
              f"[{s_ens['lo']:+.0%}, {s_ens['hi']:+.0%}], "
              f"t={s_ens['fit']['t']:.2f}, "
              f"{'significant' if s_ens['significant'] else 'NOT significant'}")
        if season == "annual" and severity in PROSSER_ANNUAL_REL:
            pr = PROSSER_ANNUAL_REL[severity]
            inside = s_ens["lo"] <= pr <= s_ens["hi"]
            print(f"      Prosser {pr:+.0%}  ratio {s_ens['rel'] / pr:.2f}  "
                  f"inside our 95% CI: {'yes' if inside else 'NO'}")

        p2, s_diag = figure_per_diagnostic(years, data, season, severity,
                                           figdir, diag_id)
        n_sig = sum(1 for s in s_diag.values() if s["significant"])
        print(f"    {p2.name:<48} {diag_id}")
        print(f"      {n_sig}/{len(s_diag)} per-diagnostic trends significant")

        # The direct comparison: only annual MOG has a published counterpart.
        if (season, severity) == ("annual", "moderate"):
            try:
                sys.path.insert(0, str(Path(__file__).resolve().parent))
                from prosser_published import (FIGURE4, FIGURE4_SUMMARY,
                                               FIGURE4_YLIM)
            except ImportError:
                print("      (ada/prosser_published.py not found — skipping "
                      "the side-by-side layout)")
            else:
                for match in (True, False):
                    p3, offscale = figure_prosser_layout(
                        years, data, figdir, FIGURE4, FIGURE4_YLIM,
                        match_axes=match)
                    print(f"    {p3.name:<48} "
                          f"{'his axes' if match else 'our axes'}")
                print(f"      Prosser states {FIGURE4_SUMMARY['n_significant']}"
                      f"/{FIGURE4_SUMMARY['n_total']} significant, max "
                      f"{FIGURE4_SUMMARY['max_change']:+.1f}%; "
                      f"ours {n_sig}/{len(s_diag)}, max "
                      f"{max(v['rel'] for v in s_diag.values()) * 100:+.1f}%")
                if offscale:
                    material = {k: v for k, v in offscale.items()
                                if v["fit_out"] or v["n_out"] > 0.1 * len(years)}
                    print(f"      PANELS THAT DO NOT FIT HIS AXIS: "
                          f"{len(material)} material, {len(offscale)} including "
                          f"one-or-two-point overhangs (his limits are read off "
                          f"a figure, so small overhangs mean nothing).")
                    for nm, d in sorted(
                            material.items(),
                            key=lambda kv: -kv[1]["n_out"]):
                        tag = ("ENTIRELY off" if d["n_out"] == len(years)
                               else f"{d['n_out']}/{len(years)} off")
                        fit = "  fitted line off too" if d["fit_out"] else ""
                        print(f"        {nm:<23}{tag:<16}"
                              f"ours {d['ours'][0]:.3f}-{d['ours'][1]:.3f}  "
                              f"his axis {d['his'][0]:.3f}-{d['his'][1]:.3f}"
                              f"{fit}")

        probs = crosscheck(s_diag, fit_csv)
        if probs:
            all_problems += [f"{season}/{severity}: {p}" for p in probs]
        elif fit_csv:
            print("      cross-check against the archived fit CSV: agrees")
        print()

    if all_problems:
        print("!! FIGURE AND ARCHIVED CSV DISAGREE — do not use these figures")
        for p in all_problems:
            print("   " + p)
        return 3

    print("done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
