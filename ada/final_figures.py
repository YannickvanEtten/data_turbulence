#!/usr/bin/env python
"""
ada/final_figures.py
====================
(FINAL-CHECKS-2026-09-16)
The three summary figures that the 21-panel layouts cannot show: how well the
final dataset reproduces Prosser across ALL published per-diagnostic cells, how
far each diagnostic's LEVEL sits from his, and where each diagnostic's global
tail lives.

    python ada/final_figures.py --outputs cat_outputs/final --tag _audit

Renders only; computes nothing from ERA5. Inputs are the archived CSVs:

    per_diagnostic_annual_<severity><tag>.csv          the fits
    per_diagnostic_annual_<severity><tag>_series.csv   the 42 points
    tail_latitude_2000<tag>.csv                        optional (step 1b)

The 21-panel Prosser-layout figures stay in ada/prosser_figures.py; this module
imports its fit so the two can never disagree, and imports the published values
from ada/prosser_published.py and ada/prosser_published_s4.py.

COLOUR
    The panel figures keep Prosser's own convention (blue crosses, red fitted
    endpoints, green trend line) because their job is to be laid beside his.
    These three figures are ours, so they use a validated categorical palette:
    blue / orange / aqua for LOG / MOG / SOG (all-pairs CVD dE 9.2, normal-vision
    24.0 on a light surface), a blue-red diverging pair with a neutral grey
    midpoint for the tail-density ratio, and text in ink rather than in series
    colour. Every series is also direct-labelled, so identity is never carried
    by colour alone.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))

from final_scorecard import GLOBAL_PCT, LABEL, TCRIT_DF40, his, read  # noqa: E402
from prosser_figures import fit_summary, read_series                  # noqa: E402
from prosser_published import FIGURE4, USES_VERTICAL_STENCIL          # noqa: E402

SEVS = ["light", "moderate", "severe"]
SEV_COLOUR = {"light": "#2a78d6", "moderate": "#eb6834", "severe": "#1baf7a"}
SEV_MARKER = {"light": "o", "moderate": "s", "severe": "^"}
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
DIVERGE_LOW, DIVERGE_MID, DIVERGE_HIGH = "#2a78d6", "#f0efec", "#e34948"
STENCIL_AUDIT = dict(USES_VERTICAL_STENCIL, magnitude_pv=True)

# Class as FINAL_DATASET.md 5 assigns it. Only used to mark points, never to
# compute anything.
CLASS_C = ["negative_richardson", "colson_panofsky", "f2d", "ubf", "ncsu1",
           "endlich"]


def load(outputs: Path, tag: str) -> dict[str, dict]:
    """{severity: {name: {ours..., his...}}} for every severity present.

    The stencil grouping depends on the series: under four-variable provenance
    `magnitude_pv` is Sharman A18, which contains d(theta)/dp and therefore
    uses the 175/225 stencil; on the baseline series it is ERA5's archived
    Ertel PV, computed on a single level. The tag decides which map is used,
    so a baseline render is never labelled with the audit grouping.
    """
    stencil = STENCIL_AUDIT if "audit" in tag else USES_VERTICAL_STENCIL
    out = {}
    for sev in SEVS:
        fit = read(outputs / f"per_diagnostic_annual_{sev}{tag}.csv")
        if not fit:
            continue
        rows = {}
        for n in FIGURE4:
            o, h = fit[n], his(n, sev)
            rows[n] = dict(
                level=o["level_1979"] * 100, change=o["change"],
                lo=o["ci_lo"], hi=o["ci_hi"], t=o["t"],
                his_level=h["level"], his_rel=h["rel"], his_p=h["p"],
                ratio=o["level_1979"] * 100 / h["level"],
                share=o["level_1979"] * 100 / GLOBAL_PCT[sev],
                his_share=h["level"] / GLOBAL_PCT[sev],
                inside=o["ci_lo"] <= h["rel"] <= o["ci_hi"],
                sig=abs(o["t"]) > TCRIT_DF40, his_sig=h["p"] < 0.05,
                stencil=stencil[n])
        out[sev] = rows
    return out


def _style(ax):
    ax.grid(alpha=0.22, linewidth=0.5, color=MUTED)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)


def figure_fit_vs_published(data, figdir: Path, tag: str):
    """Our fitted 1979-2020 change against Prosser's printed one, all cells.

    One point per diagnostic per severity, with our 95 % interval as a bar.
    A point on the 1:1 line reproduces his number exactly; a bar crossing the
    line contains it. The cells that do not are labelled, because those are
    the only ones a reader needs to look up.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.6, 7.0))
    lim = [-0.25, 1.25]
    ax.plot(lim, lim, "-", color=MUTED, linewidth=1, zorder=1)
    ax.text(1.18, 1.21, "1:1", color=MUTED, fontsize=9, ha="right")
    n_in = n_tot = 0
    for sev in SEVS:
        if sev not in data:
            continue
        rows = data[sev]
        x = np.array([r["his_rel"] for r in rows.values()])
        y = np.array([r["change"] for r in rows.values()])
        lo = np.array([r["lo"] for r in rows.values()])
        hi = np.array([r["hi"] for r in rows.values()])
        c = SEV_COLOUR[sev]
        ax.vlines(x, lo, hi, color=c, linewidth=1.1, alpha=0.45, zorder=2)
        ax.scatter(x, y, s=42, marker=SEV_MARKER[sev], facecolor=c,
                   edgecolor="white", linewidth=0.8, zorder=3,
                   label=f"{LABEL[sev]} ({sum(r['inside'] for r in rows.values())}/21 contain his value)")
        n_in += sum(r["inside"] for r in rows.values())
        n_tot += len(rows)
        for name, r in rows.items():
            if not r["inside"]:
                ax.annotate(f"{name} ({LABEL[sev]})", (r["his_rel"], r["change"]),
                            textcoords="offset points", xytext=(9, -3),
                            fontsize=8.5, color=INK)
    _style(ax)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:+.0%}")
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:+.0%}")
    ax.set_xlabel("Prosser's printed change, 1979 → 2020", color=INK2, fontsize=10)
    ax.set_ylabel("our fitted change, with 95 % interval", color=INK2, fontsize=10)
    ax.set_title(f"Per-diagnostic replication, annual, North Atlantic box\n"
                 f"{n_in} of {n_tot} published cells lie inside our 95 % interval",
                 fontsize=12, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)
    fig.tight_layout()
    out = figdir / f"final_fit_vs_published{tag}.png"
    fig.savefig(out, dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    return out


def figure_level_ratio(data, figdir: Path, tag: str):
    """Level ratio (ours / his 1979 level) per diagnostic, three severities.

    Sorted by the MOG ratio. The vertical rule at 1.0 is agreement; the shaded
    band is 0.75-1.33, the band FINAL_DATASET.md 4.3 uses when deciding
    whether a disagreement is material. Diagnostics that use the 175/225 hPa
    vertical stencil are named in bold-ish ink, the single-level ones in grey,
    so the stencil grouping is visible without a second chart.
    """
    import matplotlib.pyplot as plt

    ref = data.get("moderate") or data[next(iter(data))]
    order = sorted(ref, key=lambda k: ref[k]["ratio"])
    fig, ax = plt.subplots(figsize=(8.4, 8.2))
    ax.axvspan(0.75, 1.33, color="#eceae4", zorder=0)
    ax.axvline(1.0, color=MUTED, linewidth=1, zorder=1)
    ys = np.arange(len(order))
    for sev in SEVS:
        if sev not in data:
            continue
        r = [data[sev][n]["ratio"] for n in order]
        ax.scatter(r, ys, s=44, marker=SEV_MARKER[sev], zorder=3,
                   facecolor=SEV_COLOUR[sev], edgecolor="white", linewidth=0.8,
                   label=LABEL[sev])
    ax.set_yticks(ys)
    ax.set_yticklabels(
        [f"{n}{'  ·S' if ref[n]['stencil'] else ''}" for n in order], fontsize=9)
    for lbl, n in zip(ax.get_yticklabels(), order):
        lbl.set_color(INK if ref[n]["stencil"] else MUTED)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4, 8])
    ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4", "8"])
    _style(ax)
    ax.set_xlabel("our 1979 exceedance level ÷ Prosser's  (log scale)",
                  color=INK2, fontsize=10)
    ax.set_title("How far each diagnostic's LEVEL sits from Prosser's\n"
                 "·S = uses the 175/225 hPa vertical stencil; shaded band 0.75–1.33",
                 fontsize=12, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)
    fig.tight_layout()
    out = figdir / f"final_level_ratio{tag}.png"
    fig.savefig(out, dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    return out


def figure_tail_latitude(csv_path: Path, figdir: Path, tag: str,
                         severity: str = "severe"):
    """Where each diagnostic's global tail lives, at one severity.

    Density = band share of the cos-weighted exceedance / band share of the
    area. 1.00 is an evenly spread tail, so the scale diverges about 1.00 on a
    log axis: blue below, red above, neutral at agreement.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    with Path(csv_path).open() as fh:
        rows = [r for r in csv.DictReader(fh) if r["severity"] == severity]
    if not rows:
        raise SystemExit(f"{csv_path} has no rows for severity {severity!r}")
    bands = [k[len("density_"):] for k in rows[0] if k.startswith("density_")]
    names = [r["diagnostic"] for r in rows]
    m = np.array([[float(r[f"density_{b}"]) for b in bands] for r in rows])
    order = np.argsort(-m[:, [bands.index("30N-60N")]].ravel())
    m, names = m[order], [names[i] for i in order]

    cmap = LinearSegmentedColormap.from_list(
        "tail", [DIVERGE_LOW, DIVERGE_MID, DIVERGE_HIGH])
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-3.0, vmax=3.0)
    fig, ax = plt.subplots(figsize=(8.0, 8.6))
    im = ax.imshow(np.log2(np.clip(m, 1e-3, None)), cmap=cmap, norm=norm,
                   aspect="auto")
    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels(bands, fontsize=9, rotation=30, ha="right")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    for lbl, n in zip(ax.get_yticklabels(), names):
        lbl.set_color(INK if n in CLASS_C else INK2)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7.5,
                    color=INK if 0.35 < v < 2.8 else "white")
    ax.tick_params(colors=INK2, length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02,
                      ticks=[-2, -1, 0, 1, 2])
    cb.ax.set_yticklabels(["0.25×", "0.5×", "even", "2×", "4×"], fontsize=9,
                          color=INK2)
    cb.outline.set_visible(False)
    ax.set_title(f"Where each diagnostic's {LABEL[severity]} tail lives — "
                 f"reference year 2000\n"
                 f"density = share of the cos-weighted exceedance ÷ share of "
                 f"the area; rows sorted by the 30–60°N band",
                 fontsize=12, color=INK, loc="left")
    fig.tight_layout()
    out = figdir / f"final_tail_latitude_{severity}{tag}.png"
    fig.savefig(out, dpi=200)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    return out


def table_rows(data) -> list[dict]:
    """The same numbers as the figures, for a table view (the relief rule)."""
    rows = []
    for sev in SEVS:
        for n, r in data.get(sev, {}).items():
            rows.append(dict(severity=LABEL[sev], diagnostic=n,
                             stencil="S" if r["stencil"] else "",
                             level_pct=round(r["level"], 4),
                             his_level_pct=r["his_level"],
                             ratio=round(r["ratio"], 2),
                             na_share=round(r["share"], 2),
                             his_na_share=round(r["his_share"], 2),
                             change=round(r["change"], 4),
                             ci=f"[{r['lo']:+.0%}, {r['hi']:+.0%}]",
                             t=round(r["t"], 2), his_rel=round(r["his_rel"], 4),
                             his_p=r["his_p"], inside=r["inside"],
                             sig_agree=r["sig"] == r["his_sig"]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outputs", default=str(REPO / "cat_outputs" / "final"))
    ap.add_argument("--tag", default="_audit")
    ap.add_argument("--figdir", default=None)
    ap.add_argument("--tail-csv", default=None)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")

    outputs = Path(args.outputs)
    figdir = Path(args.figdir) if args.figdir else outputs / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    data = load(outputs, args.tag)
    if not data:
        print(f"no per_diagnostic_annual_*{args.tag}.csv under {outputs}")
        return 1
    print(f"severities found: {', '.join(LABEL[s] for s in data)}")
    print(figure_fit_vs_published(data, figdir, args.tag).name)
    print(figure_level_ratio(data, figdir, args.tag).name)
    tail = Path(args.tail_csv) if args.tail_csv else \
        outputs / f"tail_latitude_2000{args.tag}.csv"
    if tail.exists():
        for sev in ("light", "severe"):
            print(figure_tail_latitude(tail, figdir, args.tag, sev).name)
    else:
        print(f"(no {tail.name} — skipping the tail-latitude figure)")
    with (figdir / f"final_summary_table{args.tag}.csv").open("w", newline="") as fh:
        rows = table_rows(data)
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print(f"final_summary_table{args.tag}.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
