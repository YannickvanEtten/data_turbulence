#!/usr/bin/env python
"""
ada/final_scorecard.py
======================
(FINAL-CHECKS-2026-09-16)
The per-diagnostic scorecard on the FINAL series, at LOG, MOG and SOG, against
every exact per-diagnostic number Prosser (2023) prints: main-text Figure 4
(MOG) and Supporting Information Figures S4-LOG and S4-SOG.

    python ada/final_scorecard.py                       # reads cat_outputs/final/
    python ada/final_scorecard.py --light X --moderate Y --severe Z

No cluster, no xarray: CSV in, CSV and Markdown out. Runs on ADA or Windows.

WHAT IT REPORTS, PER DIAGNOSTIC AND SEVERITY
    ours: fitted 1979 level, fitted change, 95 % CI, t
    his : read-off 1979 level (approximate), printed rel and p (exact)
    level ratio (ours / his) and NORTH-ATLANTIC SHARE (level / global
    percentile target: 3 %, 0.4 %, 0.1 %) for both
    whether his rel lies inside our CI; whether the two significance calls agree

WHY THE NA SHARE
    Thresholds are GLOBAL percentiles, so by construction the global
    exceedance is exactly 3 / 0.4 / 0.1 %. The North-Atlantic level divided by
    that number is the fraction of the global tail that falls in the box,
    relative to an even spread. A constant factor cannot move it; only a
    change in WHERE a diagnostic's tail lives can. That makes it the right
    quantity for comparing a diagnostic's behaviour across severities.

WHAT IT DOES NOT DO
    It prints no verdict. STATUS 12.8 / 15.6: automated pass/fail on n = 42 has
    flipped on rounding twice in this project. The two CONSISTENCY checks below
    are equality tests on identical inputs, not statistical judgements.

CONSISTENCY CHECKS (printed first)
    1. regression: the MOG table written by ada/per_diagnostic_all_severities.py
       must equal the existing cat_outputs/per_diagnostic_annual_moderate_audit.csv
       (job 1108148, ada/per_diagnostic_trend.py) to the printed precision.
    2. the five diagnostics that are bit-identical between the two conventions
       (STATUS 18.13 / 18.14) must give identical rows in the audit and the
       baseline tables, at every severity for which both exist.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))

from prosser_published import FIGURE4, USES_VERTICAL_STENCIL  # noqa: E402
from prosser_published_s4 import S4_LOG, S4_SOG, TABLE1_ANNUAL_1979_PCT  # noqa: E402

SEVS = ["light", "moderate", "severe"]
LABEL = {"light": "LOG", "moderate": "MOG", "severe": "SOG"}
GLOBAL_PCT = {"light": 3.0, "moderate": 0.4, "severe": 0.1}
TCRIT_DF40 = 2.021
IDENTICAL_FIVE = ["colson_panofsky", "endlich", "negative_richardson",
                  "vertical_wind_shear", "wind_speed"]
COLS = ["level_1979", "change", "ci_lo", "ci_hi", "t"]
# FINAL_DATASET.md section 5, class C: kept in the 21-member ensemble, not
# recommended as stand-alone indicators. The scorecard reports how much the
# ensemble trend depends on them.
CLASS_C = ["negative_richardson", "colson_panofsky", "f2d", "ubf", "ncsu1"]

# Under four-variable provenance, magnitude_pv is Sharman A18, which contains
# d(theta)/dp -- so on the AUDIT series it uses the 175/225 stencil too
# (STATUS 18.12, audit F3 warning). prosser_published.py's split was written
# for the archived-PV baseline.
STENCIL_AUDIT = dict(USES_VERTICAL_STENCIL, magnitude_pv=True)


def his(name: str, sev: str) -> dict:
    if sev == "moderate":
        panel, _t, rel, _abs, p, lvl = FIGURE4[name]
    else:
        panel, rel, p, lvl = (S4_LOG if sev == "light" else S4_SOG)[name]
    return dict(panel=panel, rel=rel / 100.0, p=p, level=lvl)


def read(path: Path | None) -> dict | None:
    if path is None or not path.exists():
        return None
    with path.open() as fh:
        return {r["diagnostic"]: {k: float(r[k]) for k in COLS}
                for r in csv.DictReader(fh)}


def ensemble_sensitivity(table_path: Path) -> str:
    """Fitted 1979->2020 relative change of the mean of the 21 annual series,
    and of the mean without class C. Uses the _series.csv beside the table."""
    sp = table_path.with_name(table_path.stem + "_series" + table_path.suffix)
    if not sp.exists():
        return "series file not found"
    with sp.open() as fh:
        rows = list(csv.DictReader(fh))
    names = [k for k in rows[0] if k not in ("year", "ensemble")]
    x = [float(r["year"]) for r in rows]

    def rel(keep):
        y = [statistics.fmean(float(r[n]) for n in keep) for r in rows]
        slope, icpt = statistics.linear_regression(x, y)
        y0, y1 = slope * 1979 + icpt, slope * 2020 + icpt
        return (y1 - y0) / y0

    all21 = rel(names)
    wo = rel([n for n in names if n not in CLASS_C])
    return (f"mean of 21 {all21:+.1%}, without class C ({len(CLASS_C)}) {wo:+.1%}, "
            f"difference {100 * (wo - all21):+.1f} pp")


def era51_block_test(table_path: Path) -> list[str]:
    """ERA5 (not ERA5.1) is used for 2000-2006. Fit y = a + b (t - tbar) + c D,
    D = 1 in 2000-2006, per series. Report c as a fraction of the mean level,
    its t, and the fitted 1979->2020 change with and without D. A block of
    seven years centred near the sample mean has low leverage on the slope
    (the change moves by 0.163 x the block bias), so this bounds how much an
    ERA5/ERA5.1 difference confined to those years could move the trend.
    Interannual variability (NAO) produces block anomalies too: this is a
    bound and a sign pattern, not an attribution."""
    import numpy as np
    sp = table_path.with_name(table_path.stem + "_series" + table_path.suffix)
    if not sp.exists():
        return ["series file not found"]
    with sp.open() as fh:
        rows = list(csv.DictReader(fh))
    x = np.array([float(r["year"]) for r in rows])
    d = ((x >= 2000) & (x <= 2006)).astype(float)
    if d.sum() == 0 or d.sum() == len(d) or len(d) < 5:
        return ["not applicable: the series must contain 2000-2006 and other years"]
    X3 = np.column_stack([np.ones_like(x), x - x.mean(), d])
    out = []
    for n in [k for k in rows[0] if k != "year"]:
        y = np.array([float(r[n]) for r in rows])
        b, *_ = np.linalg.lstsq(X3, y, rcond=None)
        e = y - X3 @ b
        cov = (e @ e / (len(y) - 3)) * np.linalg.inv(X3.T @ X3)
        t_c = b[2] / math.sqrt(cov[2, 2]) if cov[2, 2] > 0 else float("nan")
        # change 1979 -> 2020 from the dummy model (D = 0 at both ends)
        y0 = b[0] + b[1] * (1979 - x.mean())
        y1 = b[0] + b[1] * (2020 - x.mean())
        c1 = np.polyfit(x, y, 1)
        z0, z1 = c1[0] * 1979 + c1[1], c1[0] * 2020 + c1[1]
        out.append(f"{n:<22} block {b[2] / b[0]:+7.1%} (t {t_c:+5.2f})   change "
                   f"{(z1 - z0) / z0:+6.1%} -> {(y1 - y0) / y0:+6.1%} with the 2000-06 dummy")
    return out


def compare_tables(a: dict, b: dict, names) -> float:
    """Largest absolute difference over the printed columns. The CSVs print
    level with 8 decimals, change/ci with 6, t with 4, so a difference below
    half the last printed unit is 'identical to printed precision'."""
    tol = {"level_1979": 5e-9, "change": 5e-7, "ci_lo": 5e-7, "ci_hi": 5e-7,
           "t": 5e-5}
    worst = 0.0
    for n in names:
        for c in COLS:
            d = abs(a[n][c] - b[n][c])
            worst = max(worst, d / tol[c])
    return worst          # in units of "half the last printed digit"


def main() -> int:
    ap = argparse.ArgumentParser()
    final = REPO / "cat_outputs" / "final"
    ap.add_argument("--light", type=Path,
                    default=final / "per_diagnostic_annual_light_audit.csv")
    ap.add_argument("--moderate", type=Path,
                    default=final / "per_diagnostic_annual_moderate_audit.csv")
    ap.add_argument("--severe", type=Path,
                    default=final / "per_diagnostic_annual_severe_audit.csv")
    ap.add_argument("--existing-audit-moderate", type=Path,
                    default=REPO / "cat_outputs" / "per_diagnostic_annual_moderate_audit.csv")
    ap.add_argument("--baseline-moderate", type=Path,
                    default=REPO / "data_prosser" / "per_diagnostic_annual_moderate.csv")
    ap.add_argument("--baseline-severe", type=Path,
                    default=REPO / "data_prosser" / "per_diagnostic_annual_severe.csv")
    ap.add_argument("--series-label", default="audit",
                    help="'audit' (default) assigns magnitude_pv to the stencil "
                         "group; anything else uses prosser_published's split")
    ap.add_argument("--out", type=Path, default=final / "FINAL_SCORECARD")
    args = ap.parse_args()

    ours = {s: read(getattr(args, s)) for s in SEVS}
    base = {"moderate": read(args.baseline_moderate),
            "severe": read(args.baseline_severe)}
    stencil = STENCIL_AUDIT if args.series_label == "audit" else USES_VERTICAL_STENCIL

    print("=" * 100)
    print(f"FINAL PER-DIAGNOSTIC SCORECARD — series label '{args.series_label}'")
    for s in SEVS:
        p = getattr(args, s)
        print(f"   {LABEL[s]}  {p}  {'OK' if ours[s] else 'MISSING'}")
    print("=" * 100)

    # ---------------------------------------------------------- consistency
    print("\nCONSISTENCY 1 — new MOG table vs the existing audit MOG table")
    old = read(args.existing_audit_moderate)
    if ours["moderate"] and old:
        w = compare_tables(ours["moderate"], old, FIGURE4)
        print(f"   largest difference: {w:.2f} x half-last-printed-digit "
              f"({'identical to printed precision' if w <= 1 else 'DIFFERENT — read before using anything below'})")
    else:
        print("   not available (one of the two files is missing)")

    print("\nCONSISTENCY 2 — the five convention-invariant diagnostics, audit vs baseline")
    for s in ("moderate", "severe"):
        if ours[s] and base[s]:
            w = compare_tables(ours[s], base[s], IDENTICAL_FIVE)
            print(f"   {LABEL[s]}: largest difference {w:.2f} x half-last-printed-digit "
                  f"({'identical to printed precision' if w <= 1 else 'DIFFERENT'})")
        else:
            print(f"   {LABEL[s]}: not available")

    # ----------------------------------------------------- read-off sanity
    print("\nREAD-OFF VALIDATION — mean of Prosser's 21 read-off 1979 levels vs his exact Table 1")
    for s in SEVS:
        m = statistics.fmean(his(n, s)["level"] for n in FIGURE4)
        print(f"   {LABEL[s]}: {m:.4f} % vs {TABLE1_ANNUAL_1979_PCT[s]:.4f} %  "
              f"({(m / TABLE1_ANNUAL_1979_PCT[s] - 1):+.1%})")

    rows = []
    for s in SEVS:
        tab = ours[s]
        if not tab:
            continue
        print("\n" + "=" * 100)
        print(f"{LABEL[s]} — annual, n = 42")
        print("=" * 100)
        print(f"   {'diagnostic':<22}{'st':>3}{'lvl ours':>9}{'his':>8}{'ratio':>7}"
              f"{'share o':>8}{'his':>6}{'chg ours':>9}{'95% CI':>15}{'t':>6}"
              f"{'his rel':>8}{'p':>8}{'in':>4}{'sig':>6}")
        for n in sorted(FIGURE4, key=lambda k: FIGURE4[k][0]):
            o, h = tab[n], his(n, s)
            lvl = o["level_1979"] * 100
            ratio = lvl / h["level"]
            sig_o = abs(o["t"]) > TCRIT_DF40
            sig_h = h["p"] < 0.05
            inside = o["ci_lo"] <= h["rel"] <= o["ci_hi"]
            rows.append(dict(severity=LABEL[s], diagnostic=n, panel=h["panel"],
                             uses_stencil=int(stencil[n]),
                             level_ours_pct=lvl, level_his_pct=h["level"],
                             level_ratio=ratio,
                             na_share_ours=lvl / GLOBAL_PCT[s],
                             na_share_his=h["level"] / GLOBAL_PCT[s],
                             change_ours=o["change"], ci_lo=o["ci_lo"],
                             ci_hi=o["ci_hi"], t_ours=o["t"],
                             rel_his=h["rel"], p_his=h["p"],
                             his_inside_our_ci=int(inside),
                             sig_ours=int(sig_o), sig_his=int(sig_h),
                             sig_agree=int(sig_o == sig_h)))
            print(f"   {n:<22}{'S' if stencil[n] else '-':>3}{lvl:>8.3f}%{h['level']:>7.3f}%"
                  f"{ratio:>7.2f}{lvl / GLOBAL_PCT[s]:>8.2f}{h['level'] / GLOBAL_PCT[s]:>6.2f}"
                  f"{o['change']:>+9.0%}  [{o['ci_lo']:>+4.0%},{o['ci_hi']:>+5.0%}]"
                  f"{o['t']:>6.2f}{h['rel']:>+8.0%}{h['p']:>8.0e}"
                  f"{'y' if inside else 'N':>4}"
                  f"{('=' if sig_o == sig_h else ('o' if sig_o else 'h')):>6}")
        sub = [r for r in rows if r["severity"] == LABEL[s]]
        n_o = sum(r["sig_ours"] for r in sub)
        n_h = sum(r["sig_his"] for r in sub)
        n_in = sum(r["his_inside_our_ci"] for r in sub)
        n_ag = sum(r["sig_agree"] for r in sub)
        diff = [r["diagnostic"] for r in sub if not r["sig_agree"]]
        outside = [r["diagnostic"] for r in sub if not r["his_inside_our_ci"]]
        print(f"\n   significant: ours {n_o}/21, his {n_h}/21; same call on {n_ag}/21"
              f"{'  (differ: ' + ', '.join(diff) + ')' if diff else ''}")
        print(f"   his rel inside our 95% CI: {n_in}/21"
              f"{'  (outside: ' + ', '.join(outside) + ')' if outside else ''}")
        st = [r["level_ratio"] for r in sub if r["uses_stencil"]]
        sl = [r["level_ratio"] for r in sub if not r["uses_stencil"]]
        print(f"   median level ratio: stencil group ({len(st)}) {statistics.median(st):.2f}, "
              f"single-level ({len(sl)}) {statistics.median(sl):.2f}")
        his_mean = statistics.fmean(r["level_his_pct"] for r in sub)
        gap = {r["diagnostic"]: (r["level_ours_pct"] - r["level_his_pct"]) / 21 / his_mean
               for r in sub}
        g_st = sum(v for k, v in gap.items() if stencil[k])
        g_sl = sum(v for k, v in gap.items() if not stencil[k])
        print(f"   ensemble level gap (mean of 21 / his mean - 1): {sum(gap.values()):+.3f}"
              f"  = stencil group {g_st:+.3f} + single-level group {g_sl:+.3f}")
        try:
            sens = ensemble_sensitivity(getattr(args, s))
        except Exception as exc:
            sens = f"failed: {type(exc).__name__}: {exc}"
        print(f"   ensemble trend sensitivity: {sens}")
        print("   ERA5-vs-ERA5.1 block test (2000-2006 dummy):")
        try:
            lines = era51_block_test(getattr(args, s))
        except Exception as exc:                  # an extra, never fatal
            lines = [f"failed: {type(exc).__name__}: {exc}"]
        for line in lines:
            print("      " + line)
        top = sorted(gap.items(), key=lambda kv: kv[1])
        print("   largest contributors: " + ", ".join(f"{k} {v:+.3f}" for k, v in top[:4])
              + " ... " + ", ".join(f"{k} {v:+.3f}" for k, v in top[-2:]))

    if not rows:
        print("\nno input tables found — nothing written")
        return 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.with_suffix(".csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                         for k, v in r.items()})
    with args.out.with_suffix(".md").open("w") as fh:
        fh.write("# Final per-diagnostic scorecard\n\nGenerated by "
                 "`ada/final_scorecard.py`. Levels of Prosser are read off his "
                 "figures (approximate); his rel and p are printed (exact). "
                 "No verdicts.\n\n")
        fh.write("| sev | diagnostic | S | level ours % | his % | ratio | NA share ours | his "
                 "| change ours | 95% CI | t | his rel | his p | inside | same sig |\n")
        fh.write("|" + "---|" * 15 + "\n")
        for r in rows:
            fh.write(f"| {r['severity']} | {r['diagnostic']} | {'S' if r['uses_stencil'] else ''} "
                     f"| {r['level_ours_pct']:.3f} | {r['level_his_pct']:.3f} "
                     f"| {r['level_ratio']:.2f} | {r['na_share_ours']:.2f} "
                     f"| {r['na_share_his']:.2f} | {r['change_ours']:+.0%} "
                     f"| [{r['ci_lo']:+.0%}, {r['ci_hi']:+.0%}] | {r['t_ours']:.2f} "
                     f"| {r['rel_his']:+.0%} | {r['p_his']:.0e} "
                     f"| {'yes' if r['his_inside_our_ci'] else 'NO'} "
                     f"| {'yes' if r['sig_agree'] else 'NO'} |\n")
    print(f"\nwrote {args.out.with_suffix('.csv')}\nwrote {args.out.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
