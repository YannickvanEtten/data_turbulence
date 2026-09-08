#!/usr/bin/env python
"""
Prosser (2023) replication scorecard.

One row per published item, with our value beside theirs and an explicit
verdict — including for the items we CANNOT compute, which are named rather
than omitted, because a scorecard that silently drops what it cannot do is
worse than no scorecard.

    python ada/prosser_scorecard.py --log logs/full-trend-1101414.out

Three sources, all cheap, none of which touch ERA5:

  1. Prosser's published values, hardcoded below with their source line.
     Table 1 is the only thing he tabulates; everything else in the paper is a
     figure, so those comparisons are read off an image and are sign-and-
     magnitude checks, not precision ones. The scorecard says so per row.
  2. `logs/full-trend-*.out` — the MASTER SUMMARY block written by
     ada/full_trend_check.py, which is where the 25 season x severity cells
     live. Parsed rather than recomputed: the run costs 6 h 20 m.
  3. `cat_outputs/per_diagnostic_*.csv` — the per-diagnostic fits from
     ada/per_diagnostic_trend.py.

Writes cat_outputs/scorecard.csv and cat_outputs/SCORECARD.md.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

SEV_ORDER = ["LOG", "LMOG", "MOG", "MSOG", "SOG"]
SEASON_ORDER = ["DJF", "MAM", "JJA", "SON", "Annual"]

# Prosser (2023) Table 1 — fitted relative increase 1979->2020, diagnostic
# mean, North Atlantic box (36-60N, 55-10W), ERA5 at 197 hPa.
# Source: CALIBRATION_REFERENCE.md §4.4, read from the paper.
PROSSER_TABLE1 = {
    "DJF":    {"LOG": 0.21, "LMOG": 0.30, "MOG": 0.37, "MSOG": 0.43, "SOG": 0.49},
    "MAM":    {"LOG": 0.26, "LMOG": 0.43, "MOG": 0.57, "MSOG": 0.71, "SOG": 0.85},
    "JJA":    {"LOG": 0.09, "LMOG": 0.20, "MOG": 0.31, "MSOG": 0.41, "SOG": 0.52},
    "SON":    {"LOG": 0.15, "LMOG": 0.23, "MOG": 0.31, "MSOG": 0.39, "SOG": 0.47},
    "Annual": {"LOG": 0.17, "LMOG": 0.28, "MOG": 0.37, "MSOG": 0.46, "SOG": 0.55},
}
PROSSER_TABLE1_HOURS = {  # (1979 h, 2020 h)
    "DJF":    {"LOG": (128.9, 155.6), "LMOG": (45.6, 59.3), "MOG": (22.3, 30.6),
               "MSOG": (12.1, 17.2), "SOG": (6.4, 9.6)},
    "MAM":    {"LOG": (90.4, 113.4), "LMOG": (27.2, 38.9), "MOG": (11.8, 18.6),
               "MSOG": (5.7, 9.7), "SOG": (2.7, 5.0)},
    "JJA":    {"LOG": (114.1, 124.5), "LMOG": (36.5, 43.8), "MOG": (16.1, 21.1),
               "MSOG": (7.7, 10.9), "SOG": (3.6, 5.5)},
    "SON":    {"LOG": (133.1, 153.2), "LMOG": (43.4, 53.4), "MOG": (19.8, 25.8),
               "MSOG": (10.0, 13.9), "SOG": (5.0, 7.4)},
    "Annual": {"LOG": (466.5, 546.8), "LMOG": (152.7, 195.4), "MOG": (70.0, 96.1),
               "MSOG": (35.5, 51.8), "SOG": (17.7, 27.4)},
}

# Prosser's own statement about Figure 4 (main text, MOG, North Atlantic):
# "...significant (p < 0.05) upward trends, with relative changes of up to
#  75.6%. The remaining 4 diagnostics show no significant trend, and none of
#  the diagnostics shows a significant downward trend."
PROSSER_FIG4 = dict(n_significant=17, n_total=21, max_change=0.756,
                    n_significant_downward=0)

# Everything in the paper, and what each needs. The items marked
# needs="global" are the reason full replication is a two-week job rather than
# an afternoon: they are global maps over all 42 years.
INVENTORY = [
    ("Table 1", "Fitted 1979-2020 hours, 5 seasons x 5 severities, NA box",
     "NA box", "table"),
    ("Figure 1a-d", "Global maps of annual-mean MOG probability: 1979, 2020, "
     "and both inferred from the regression", "global", "map"),
    ("Figure 2a,b", "Global maps of absolute and relative MOG change, with "
     "Wald-test stippling", "global", "map"),
    ("Figure 3a", "Annual diagnostic-mean MOG probability, NA box, 42 points",
     "NA box", "series"),
    ("Figure 3b", "Same, USA box (30-55N, 124-60W)", "USA box", "series"),
    ("Figure 4", "Per-diagnostic annual MOG probability, NA box, 21 panels",
     "NA box", "series"),
    ("Figure S1-LOG/SOG", "Figure 1 at LOG and SOG", "global", "map"),
    ("Figure S2-LOG/SOG", "Figure 2 at LOG and SOG", "global", "map"),
    ("Figure S3-LOG/SOG", "Figure 3 at LOG and SOG (both boxes)",
     "NA box + USA box", "series"),
    ("Figure S4-LOG/SOG", "Figure 4 at LOG and SOG, NA box", "NA box", "series"),
    ("Figure S5", "Figure 2a broken down by all 21 diagnostics", "global", "map"),
]

MASTER = re.compile(
    r"^\s+(DJF|MAM|JJA|SON|Annual)\s+(LOG|LMOG|MOG|MSOG|SOG)\s+"
    r"(-?\d+)%\s+(-?\d+)%\s+(-?[\d.]+|nan)\s+(-?[\d.]+|nan)\s+"
    r"(yes|no)\s+(yes|no)\s*$")


def parse_master_summary(path: Path) -> dict[tuple[str, str], dict]:
    """Pull the 25 season x severity cells out of a jobs/16 log."""
    if not path.exists():
        return {}
    text = path.read_text(errors="replace").splitlines()
    try:
        start = next(i for i, l in enumerate(text) if "MASTER SUMMARY" in l)
    except StopIteration:
        print(f"   !! {path.name} has no MASTER SUMMARY block — is this a "
              f"jobs/16 log?")
        return {}
    cells = {}
    for line in text[start:]:
        m = MASTER.match(line)
        if not m:
            continue
        season, sev, ours, theirs, ratio, t, sig, inside = m.groups()
        cells[(season, sev)] = dict(
            ours=int(ours) / 100.0, theirs=int(theirs) / 100.0,
            ratio=float(ratio) if ratio != "nan" else float("nan"),
            t=float(t) if t != "nan" else float("nan"),
            significant=(sig == "yes"), inside=(inside == "yes"))
    return cells


def read_per_diagnostic(outputs: Path, season: str, severity: str):
    p = outputs / f"per_diagnostic_{season}_{severity}.csv"
    if not p.exists():
        return None
    with p.open() as fh:
        return list(csv.DictReader(fh))


def fmt_pct(v) -> str:
    return "—" if v is None else f"{v:+.0%}"


def _ensemble_row(ens_fig: str, sev_key: str, sev_label: str,
                  cells: dict, figdir: Path) -> dict:
    """The Figure 3a / S3 row: our annual ensemble trend against Table 1."""
    c = cells.get(("Annual", sev_label))
    rendered = (figdir / f"fig3_annual_{sev_key}_ensemble.png").exists()
    return dict(
        item=ens_fig, quantity="annual diagnostic-mean series + fit",
        domain="NA box",
        ours=fmt_pct(c["ours"]) if c else "see the rendered panel",
        prosser=f"{PROSSER_TABLE1['Annual'][sev_label]:+.0%}",
        ratio=f"{c['ratio']:.2f}" if c else "",
        precision="the trend is Table 1 (exact); the 42 plotted points are "
                  "read off the figure",
        verdict=("MATCH — inside our 95% CI" if c and c["inside"]
                 else ("OUTSIDE our 95% CI" if c else "")),
        note=f"rendered: {'yes' if rendered else 'no'}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Prosser replication scorecard.")
    ap.add_argument("--log", default=None,
                    help="a jobs/16 log holding the MASTER SUMMARY block "
                         "(e.g. logs/full-trend-1101414.out)")
    ap.add_argument("--outputs", default=str(REPO / "cat_outputs"))
    ap.add_argument("--figdir", default=None,
                    help="where prosser_figures.py wrote its panels, to note "
                         "which rows have a rendered figure")
    args = ap.parse_args()

    outputs = Path(args.outputs)
    figdir = Path(args.figdir) if args.figdir else outputs / "figures"
    outputs.mkdir(parents=True, exist_ok=True)

    cells = parse_master_summary(Path(args.log)) if args.log else {}
    if args.log:
        print(f"Table 1 cells parsed from {args.log}: {len(cells)}/25")
    else:
        print("no --log given: Table 1 rows will be listed as not loaded")

    rows: list[dict] = []

    # ---------------------------------------------------------- Table 1
    for season in SEASON_ORDER:
        for sev in SEV_ORDER:
            c = cells.get((season, sev))
            theirs = PROSSER_TABLE1[season][sev]
            h79, h20 = PROSSER_TABLE1_HOURS[season][sev]
            rows.append(dict(
                item=f"Table 1 [{season}/{sev}]",
                quantity="fitted relative change 1979-2020",
                domain="NA box",
                ours=fmt_pct(c["ours"]) if c else "not loaded",
                prosser=f"{theirs:+.0%}",
                ratio=f"{c['ratio']:.2f}" if c else "",
                precision="exact (published table)",
                verdict=("" if not c else
                         ("MATCH — inside our 95% CI"
                          + ("" if c["significant"]
                             else ", but our trend is not significant"))
                         if c["inside"] else "OUTSIDE our 95% CI"),
                note=f"Prosser {h79:.1f} h -> {h20:.1f} h"))

    # ---------------------------------------------------- Figures 3a / 4 / S
    for sev_key, sev_label, ens_fig, diag_fig in [
            ("light", "LOG", "Figure S3-LOG a", "Figure S4-LOG"),
            ("moderate", "MOG", "Figure 3a", "Figure 4"),
            ("severe", "SOG", "Figure S3-SOG a", "Figure S4-SOG")]:
        recs = read_per_diagnostic(outputs, "annual", sev_key)
        have_fig = (figdir / f"fig4_annual_{sev_key}_per_diagnostic.png").exists()

        if recs is None:
            rows.append(dict(
                item=diag_fig, quantity="21 per-diagnostic annual fits",
                domain="NA box", ours="NOT RUN", prosser="", ratio="",
                precision="read off a figure (sign and magnitude only)",
                verdict="run jobs/21 at --season annual "
                        f"--severity {sev_key}",
                note=""))
            rows.append(_ensemble_row(ens_fig, sev_key, sev_label, cells, figdir))
            continue

        n = len(recs)
        n_sig = sum(1 for r in recs
                    if abs(float(r["t"])) >= 2.021)          # df=40
        changes = [float(r["change"]) for r in recs]
        n_neg_sig = sum(1 for r in recs
                        if float(r["change"]) < 0 and abs(float(r["t"])) >= 2.021)
        note = (f"range {min(changes):+.0%} to {max(changes):+.0%}; "
                f"rendered: {'yes' if have_fig else 'no'}")

        if sev_label == "MOG":
            p = PROSSER_FIG4
            verdict = []
            verdict.append(
                f"significant {n_sig}/{n} vs Prosser {p['n_significant']}/"
                f"{p['n_total']}")
            verdict.append(
                f"max change {max(changes):+.0%} vs his {p['max_change']:+.1%}")
            verdict.append(
                f"significant downward {n_neg_sig} vs his "
                f"{p['n_significant_downward']}")
            rows.append(dict(
                item=diag_fig, quantity="21 per-diagnostic annual fits",
                domain="NA box", ours=f"{n_sig}/{n} significant",
                prosser=f"{p['n_significant']}/{p['n_total']} significant",
                ratio="",
                precision="his counts are stated in the text — exact; "
                          "per-diagnostic values are read off Figure 4",
                verdict="; ".join(verdict), note=note))
        else:
            rows.append(dict(
                item=diag_fig, quantity="21 per-diagnostic annual fits",
                domain="NA box", ours=f"{n_sig}/{n} significant",
                prosser="figure only, no stated counts", ratio="",
                precision="read off a figure (sign and magnitude only)",
                verdict="computed; compare panel by panel against the SI",
                note=note))

        # the ensemble panel for the same severity
        rows.append(_ensemble_row(ens_fig, sev_key, sev_label, cells, figdir))

    # ------------------------------------------------- what we cannot do yet
    for item, what, needs, kind in INVENTORY:
        if needs == "global":
            rows.append(dict(
                item=item, quantity=what, domain="global",
                ours="NOT COMPUTED", prosser="published", ratio="",
                precision="n/a",
                verdict="needs global 42-year fields — see "
                        "PLAN_prosser_replication.md",
                note="storage is not the obstacle (~17 GB of annual "
                     "aggregates); CDS queue and compute are"))
        elif needs == "USA box":
            rows.append(dict(
                item=item, quantity=what, domain="USA box",
                ours="NOT COMPUTED", prosser="published", ratio="",
                precision="n/a",
                verdict="needs a USA-box download (~470 GB, ~1.5 days)",
                note="30-55N, 124-60W"))

    # ------------------------------------------------------------- write out
    csv_path = outputs / "scorecard.csv"
    fields = ["item", "quantity", "domain", "ours", "prosser", "ratio",
              "precision", "verdict", "note"]
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

    md_path = outputs / "SCORECARD.md"
    with md_path.open("w") as fh:
        fh.write("# Prosser (2023) replication scorecard\n\n")
        fh.write("Generated by `ada/prosser_scorecard.py`. One row per "
                 "published item.\n\n")
        fh.write("**Precision matters per row.** Table 1 is the only thing "
                 "Prosser tabulates, so those comparisons are exact. Every "
                 "figure comparison is read off an image and is a sign-and-"
                 "magnitude check.\n\n")
        n_match = sum(1 for r in rows if r["verdict"].startswith("MATCH"))
        n_out = sum(1 for r in rows if "OUTSIDE" in r["verdict"])
        n_todo = sum(1 for r in rows if r["ours"] in ("NOT COMPUTED", "NOT RUN"))
        fh.write(f"- inside our 95% CI: **{n_match}**\n")
        fh.write(f"- outside our 95% CI: **{n_out}**\n")
        fh.write(f"- not computed yet: **{n_todo}**\n\n")
        fh.write("| item | domain | ours | Prosser | ratio | verdict |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for r in rows:
            fh.write(f"| {r['item']} | {r['domain']} | {r['ours']} | "
                     f"{r['prosser']} | {r['ratio']} | {r['verdict']} |\n")
        fh.write("\n## Precision and notes\n\n")
        fh.write("| item | precision | note |\n|---|---|---|\n")
        for r in rows:
            if r["precision"] != "n/a" or r["note"]:
                fh.write(f"| {r['item']} | {r['precision']} | {r['note']} |\n")

    print(f"\nwrote {csv_path}")
    print(f"wrote {md_path}")
    print(f"\n   inside 95% CI {n_match}   outside {n_out}   "
          f"not computed {n_todo}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
