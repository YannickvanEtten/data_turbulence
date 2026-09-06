"""
ada/compare_thresholds.py
=========================
Diff two threshold sets, and translate the difference into what it actually
costs downstream.

WHY THIS EXISTS
---------------
The year-2000 calibration has been computed twice on different samples:

    thresholds_2026-09-03.json   48 days (1/9/17/25 of each month)
    thresholds_2026-09-??.json   the full reference year, all 2928 timesteps

Everything the project has published so far — including the 42-year
all-season trend in RESULT_full_trend_42yr.md — used the 48-day set. The
honest way to retire that caveat is not to quietly switch and re-run, but to
MEASURE what switching costs, so the sub-sample becomes a reported
sensitivity rather than an unexamined choice. That is the number a referee
will ask for the moment they read "days 1, 9, 17 and 25".

WHAT IT REPORTS, AND WHY NOT JUST THE THRESHOLD DIFFERENCE
----------------------------------------------------------
A threshold moving by 1% does NOT move the exceedance frequency by 1%. These
distributions are heavy-tailed, so the frequency responds with a local
elasticity

    alpha = -dP/dthreshold * threshold/P

which for CAT diagnostics runs to 8 or so around MOG. A 1% threshold error is
therefore an ~8% frequency error, and that is the quantity that propagates
into every exceedance number in the project.

alpha is estimated PER DIAGNOSTIC from its own ladder rather than assumed,
using the five (threshold, exceedance) pairs the severity ladder already
provides:

    LOG p97.0 -> 3.0%   LMOG p99.1 -> 0.9%   MOG p99.6 -> 0.4%
    MSOG p99.8 -> 0.2%  SOG p99.9 -> 0.1%

The local slope dP/dthr is taken from the neighbouring rungs by central
difference. This is deliberately done on the RAW scale, not in logs: `ubf`,
`colson_panofsky` and `negative_richardson` have negative thresholds, and a
log-log elasticity is undefined for them. A central difference is not.

WHAT A LARGE MOVE MEANS
-----------------------
Not necessarily an error. Two diagnostics are expected to move more than the
rest, for known reasons:

  * `f2d` — its 48-day threshold was computed with 25% of timesteps on a
    corrupted time stencil (PLAN_full_year_calibration.md §5): the
    non-contiguous sub-sample straddles 8-day gaps. Its move should be
    SYSTEMATIC, not sampling noise, and the contiguous year is the fix.
  * `ncsu1` — cubic in gradients, so the heaviest tail of the 21 and the
    largest sampling error on any percentile of it.

A large move in anything else is worth understanding before it is used.

USAGE
-----
    pixi run python ada/compare_thresholds.py                 # two newest
    pixi run python ada/compare_thresholds.py OLD.json NEW.json
    pixi run python ada/compare_thresholds.py --csv out.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
sys.path.insert(0, str(REPO))

import calibration                                            # noqa: E402

# Williams (2017) Table 1, adopted verbatim by Prosser (2023). Exceedance
# probability implied by each rung: (100 - percentile)/100.
LADDER = {"light": 97.0, "light_to_moderate": 99.1, "moderate": 99.6,
          "moderate_to_severe": 99.8, "severe": 99.9}
ORDER = list(LADDER)
LABEL = {"light": "LOG", "light_to_moderate": "LMOG", "moderate": "MOG",
         "moderate_to_severe": "MSOG", "severe": "SOG"}
EXCEED = {k: (100.0 - v) / 100.0 for k, v in LADDER.items()}


def local_dP_dthr(thr: dict[str, float], sev: str) -> float | None:
    """dP/dthreshold at rung `sev`, by central difference on its neighbours.

    Raw scale, not log: three of the 21 diagnostics carry negative thresholds
    and a log-log elasticity is undefined for them.
    """
    i = ORDER.index(sev)
    lo, hi = ORDER[max(0, i - 1)], ORDER[min(len(ORDER) - 1, i + 1)]
    dthr = thr[hi] - thr[lo]
    if dthr == 0 or not np.isfinite(dthr):
        return None
    return (EXCEED[hi] - EXCEED[lo]) / dthr


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    ap.add_argument("paths", nargs="*", help="OLD.json NEW.json; default = two newest")
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--csv", default=None, help="also write the full table here")
    ap.add_argument("--flag-above", type=float, default=5.0,
                    help="call out any implied frequency change above this %%")
    args = ap.parse_args()

    base = Path(args.base)
    if args.paths:
        if len(args.paths) != 2:
            ap.error("give exactly two threshold files, or none")
        old_p, new_p = (Path(p) if Path(p).is_absolute() else base / p
                        for p in args.paths)
    else:
        found = sorted((base / "calibration").glob("thresholds_*.json"))
        if len(found) < 2:
            raise SystemExit(f"need two threshold files, found {len(found)} in "
                             f"{base / 'calibration'}")
        old_p, new_p = found[-2], found[-1]

    old, _, old_prov = calibration.load_thresholds(old_p)
    new, _, new_prov = calibration.load_thresholds(new_p)

    print(f"OLD  {old_p.name}   period={old_prov.get('period')}")
    print(f"NEW  {new_p.name}   period={new_prov.get('period')}")
    print()

    rows = []
    for name in sorted(set(old) & set(new)):
        for sev in ORDER:
            if sev not in old[name] or sev not in new[name]:
                continue
            o, n = float(old[name][sev]), float(new[name][sev])
            d_thr = n - o
            rel_thr = 100.0 * d_thr / abs(o) if o != 0 else np.nan
            slope = local_dP_dthr(new[name], sev)
            if slope is None:
                rel_freq = np.nan
            else:
                # dP is negative when the threshold rises: fewer exceedances.
                rel_freq = 100.0 * (slope * d_thr) / EXCEED[sev]
            rows.append((name, sev, o, n, rel_thr, rel_freq))

    missing = (set(old) ^ set(new))
    if missing:
        print(f"!! present in only one file: {sorted(missing)}\n")

    hdr = f"{'diagnostic':24s}{'sev':6s}{'old':>14s}{'new':>14s}{'d thr %':>10s}{'d freq %':>11s}"
    print(hdr)
    print("-" * len(hdr))
    for name, sev, o, n, rt, rf in rows:
        flag = ""
        if np.isfinite(rf) and abs(rf) > args.flag_above:
            flag = "  <--"
        print(f"{name:24s}{LABEL[sev]:6s}{o:14.6g}{n:14.6g}{rt:10.2f}{rf:11.2f}{flag}")

    thr_moves = np.array([abs(r[4]) for r in rows if np.isfinite(r[4])])
    frq_moves = np.array([abs(r[5]) for r in rows if np.isfinite(r[5])])
    print("\n" + "=" * len(hdr))
    print(f"threshold move   median {np.median(thr_moves):6.2f}%   "
          f"p90 {np.percentile(thr_moves, 90):6.2f}%   max {thr_moves.max():6.2f}%")
    print(f"implied freq     median {np.median(frq_moves):6.2f}%   "
          f"p90 {np.percentile(frq_moves, 90):6.2f}%   max {frq_moves.max():6.2f}%")

    # The 21 sampling errors are independent draws and aggregate.py averages
    # the 21 exceedance FIELDS, so the ensemble error falls by sqrt(21).
    n_diag = len(set(r[0] for r in rows))
    print(f"ensemble of {n_diag} (independent draws, /sqrt(n)):"
          f"   median {np.median(frq_moves)/np.sqrt(n_diag):5.2f}%"
          f"   p90 {np.percentile(frq_moves, 90)/np.sqrt(n_diag):5.2f}%")
    print()
    print("NOTE ON THE TREND: a threshold shift scales exceedance by nearly the")
    print("same factor in 1979 and 2020, and P2020/P1979 is EXACTLY invariant to")
    print("a common factor. Only the difference in local slope between the two")
    print("years survives. Expect the LEVEL to move by the figures above and the")
    print("TREND to be near-unchanged -- and confirm it by re-running jobs/16,")
    print("which is the only thing that has to be recomputed. No diagnostic")
    print("field is affected: thresholds are applied at AGGREGATION time.")

    worst = sorted((r for r in rows if np.isfinite(r[5])),
                   key=lambda r: -abs(r[5]))[:5]
    print("\nlargest implied frequency moves:")
    for name, sev, o, n, rt, rf in worst:
        note = ""
        if name == "f2d":
            note = "  (expected: 25% of its 48-day sample was on a broken time stencil)"
        elif name == "ncsu1":
            note = "  (expected: cubic in gradients, heaviest tail of the 21)"
        print(f"  {name:24s}{LABEL[sev]:6s}{rf:+8.2f}%{note}")

    if args.csv:
        out = Path(args.csv)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w") as fh:
            fh.write("diagnostic,severity,old,new,threshold_pct,implied_freq_pct\n")
            for name, sev, o, n, rt, rf in rows:
                fh.write(f"{name},{LABEL[sev]},{o!r},{n!r},{rt:.6f},{rf:.6f}\n")
        print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
