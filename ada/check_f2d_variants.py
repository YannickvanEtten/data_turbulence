#!/usr/bin/env python
"""
ada/check_f2d_variants.py
=========================
Which reading of Sharman A9 is #20? FORMULA_AUDIT.md §4.5.

    pixi run python ada/check_f2d_variants.py <raw.grib> [<raw.grib> ...] \\
        [--box williams|prosser|none]

THE QUESTION
------------
`frontogenesis_isentropic` computes

    +0.5 D/Dt[ (du/dtheta)^2 + (dv/dtheta)^2 ]                      (variant A)

which is NOT A9 as printed. A9 carries a |dv/dtheta|^-1 normalisation and a
leading minus sign; Sharman's Table B1 units support dropping the
normalisation and say nothing about the sign. FORMULA_AUDIT.md §4 sets out the
argument in full. Four readings survive:

    A   +0.5 D/Dt[Q]        current behaviour
    B   -0.5 D/Dt[Q]        A9's leading minus applied
    C   |0.5 D/Dt[Q]|       magnitude
    D   -D/Dt[sqrt(Q)]      literal A9, normalisation included

WHAT DECIDES IT
---------------
Not the paper -- the paper is internally inconsistent on this equation, which
is why the question is open at all. Two measurable things decide it:

1. DISTRIBUTION SHAPE. Williams (2017) Table 2 and Williams & Joshi (2013)
   Table 1 give p97 = 770 and median = 56.6 for frontogenesis, from the same
   model, box, season and daily-mean sampling: a ratio of 13.6. A SIGNED
   material derivative cannot produce that. D/Dt of a positive quantity in a
   statistically stationary atmosphere is centred on zero, so its median
   collapses and the ratio diverges. Variants A and B are signed; C is not;
   D is signed. If A's measured ratio comes back in the hundreds while C's
   lands near 13.6, that is decisive and it does not depend on reading the
   equation correctly.

2. THE SIGN OF THE TREND, against Prosser's Figure S5, whose frontogenesis
   panel is near-zero to slightly negative over the North Atlantic where
   variant A gives the second-largest positive trend of the 21. That test
   needs a full re-run (diagnostics -> calibration -> trend) and is step 2;
   this script is step 1 and narrows the field first.

WHY A AND B CANNOT BE TOLD APART BY SHAPE ALONE
-----------------------------------------------
B = -A exactly, so their distributions are mirror images and every symmetric
statistic agrees. What differs is WHICH TAIL a percentile threshold selects,
and therefore which grid cells are flagged -- disjoint sets. The tail-overlap
matrix below makes that concrete: if A and B overlap near 0 %, then choosing
between them is choosing between two completely different sets of turbulent
cells, and it cannot be left to a default.

RUN IT ON WHAT
--------------
A global calibration month with --box williams reproduces the published
sampling most closely. An NA month with --box prosser is cheaper and is the
box the replication actually reports on. Both are informative; the shape
verdict should not depend on which.

MEMORY -- WHY --chunk-days EXISTS HERE TOO
------------------------------------------
Same reason as ada/diagnostics_global.py. A global file is 721 x 1440 x 32 =
3.3e7 points, the diagnostics touch three pressure levels, and this script
computes two variants rather than one. Unchunked that runs well past the 24 GB
that the full 21-diagnostic global run needed -- and STATUS.md §11.9 measured a
120 GB request scheduled 22 hours out while a 24 GB one started in 13 seconds,
so an over-large request costs a day rather than some memory.

Chunking is exact here for the same reason it is there: each chunk is computed
with one extra timestep on each side and then trimmed, so every retained step
saw the neighbours it would have had in an unchunked run. Only the file's true
first and last step are dropped, and those are dropped deliberately -- d/dt is
one-sided there, and a one-sided slope on a diurnally varying field is a
different estimator from the centred one used everywhere else. Letting it into
a tail statistic would put an artefact exactly where the comparison lives.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _checkutil import (PROSSER_BOX, WILLIAMS_BOX, cos_phi_weights,  # noqa: E402
                        load_module, peak_rss_gb, spearman, subset_box)
from calib_weighted_percentile import weighted_percentile  # noqa: E402

PUBLISHED_RATIO = 770.0 / 56.6      # 13.6, Williams (2017) / Williams & Joshi (2013)
MOG_PERCENTILE = 99.6               # Prosser's MOG rung, for the tail-overlap sets
STEPS_PER_DAY = 8                   # 3-hourly ERA5
OVERLAP = 1                         # buffer steps each side; d/dt needs exactly 1


def time_dim_of(ds) -> str | None:
    for candidate in ("time", "valid_time"):
        if candidate in ds.dims:
            return candidate
    return None


def describe(values: np.ndarray, weights: np.ndarray) -> dict:
    ok = np.isfinite(values)
    v, w = values[ok], weights[ok]
    if v.size < 100:
        return dict(n=int(v.size))
    q = weighted_percentile(v, w, np.array([50.0, 97.0, MOG_PERCENTILE]))
    p50, p97, pmog = float(q[0]), float(q[1]), float(q[2])
    return dict(
        n=int(v.size),
        p50=p50, p97=p97, pmog=pmog,
        ratio=(p97 / p50) if p50 != 0 else float("inf"),
        frac_positive=float(np.average(v > 0, weights=w)),
        mean=float(np.average(v, weights=w)),
    )


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("grib", nargs="+", help="raw ERA5 GRIB file(s)")
    ap.add_argument("--box", default="williams",
                    choices=["williams", "prosser", "none"])
    ap.add_argument("--level", type=int, default=200)
    ap.add_argument("--chunk-days", type=int, default=1,
                    help="process this many days at a time, with a "
                         "1-timestep overlap buffer so the retained values are "
                         "identical to an unchunked run. Default 1, matching "
                         "ada/diagnostics_global.py, which keeps peak memory "
                         "near a North Atlantic month's ~15 GB. 0 = whole file.")
    args = ap.parse_args()

    box = {"williams": WILLIAMS_BOX, "prosser": PROSSER_BOX, "none": None}[args.box]

    diag = load_module("diagnostics", "2_diagnostics.py")

    per_file: dict[str, list[np.ndarray]] = {k: [] for k in "ABCD"}
    weights: list[np.ndarray] = []

    for path in args.grib:
        p = Path(path)
        if not p.exists():
            print(f"!! not found: {p}")
            return 1
        print(f">>> {p.name}  ({p.stat().st_size / 1e9:.2f} GB)")
        ds_raw = diag.load_era5(p)
        tname = time_dim_of(ds_raw)
        if tname is None:
            print(f"!! no recognised time dimension in {list(ds_raw.dims)}")
            return 1
        n_steps = ds_raw.sizes[tname]
        if n_steps < 3:
            print("!! fewer than 3 timesteps: d/dt is one-sided everywhere and")
            print("   the tail statistic would be an artefact. Refusing.")
            return 1

        step = args.chunk_days * STEPS_PER_DAY if args.chunk_days > 0 else n_steps
        bounds = [(s, min(s + step, n_steps)) for s in range(0, n_steps, step)]
        print(f"    {n_steps} timesteps -> {len(bounds)} chunk(s), overlap {OVERLAP}")

        for i, (b0, b1) in enumerate(bounds):
            lo, hi = max(0, b0 - OVERLAP), min(n_steps, b1 + OVERLAP)
            sub = ds_raw.isel({tname: slice(lo, hi)})
            prepared = diag.prepare_for_rojak(sub)
            ds = prepared._dataset

            # A and D are genuinely different computations; B and C are exact
            # functions of A, so computing those separately would add cost and
            # a chance of the three drifting apart.
            a = diag.frontogenesis_isentropic(ds, target_level=args.level, variant="A")
            d = diag.frontogenesis_isentropic(ds, target_level=args.level, variant="D")

            # Keep only steps that had a neighbour on BOTH sides in the
            # original file: the buffer steps this chunk borrowed, and the
            # file's own first and last step, all go.
            keep0 = max(b0, 1) - lo
            keep1 = min(b1, n_steps - 1) - lo
            if keep1 <= keep0:
                continue
            tdim = time_dim_of(a)
            a = a.isel({tdim: slice(keep0, keep1)})
            d = d.isel({tdim: slice(keep0, keep1)})

            if box is not None:
                a = subset_box(a.to_dataset(name="f2d"), **box)["f2d"]
                d = subset_box(d.to_dataset(name="f2d"), **box)["f2d"]

            av = np.asarray(a.values, dtype=np.float64).ravel()
            dv = np.asarray(d.values, dtype=np.float64).ravel()
            per_file["A"].append(av)
            per_file["B"].append(-av)
            per_file["C"].append(np.abs(av))
            per_file["D"].append(dv)
            weights.append(cos_phi_weights(a))
            print(f"    chunk {i + 1}/{len(bounds)}  steps {b0}:{b1}  "
                  f"{av.size:,} cells   [rss {peak_rss_gb():.1f} GB]")
            del sub, prepared, ds, a, d

        del ds_raw

    if not weights:
        print("!! no usable chunks — every chunk was buffer-only. Check the "
              "input file's timestep count.")
        return 1

    vals = {k: np.concatenate(v) for k, v in per_file.items()}
    w = np.concatenate(weights)
    stats = {k: describe(v, w) for k, v in vals.items()}

    # ------------------------------------------------------------------ shape
    print("\n" + "=" * 84)
    print("DISTRIBUTION SHAPE PER VARIANT")
    print("=" * 84)
    print(f"   box {args.box}, {vals['A'].size:,} cells, "
          f"published p97/median = {PUBLISHED_RATIO:.1f}\n")
    print(f"   {'':<3}{'formula':<42}{'median':>12}{'p97':>12}{'p97/med':>10}{'frac>0':>9}")
    for k in "ABCD":
        s = stats[k]
        if "p50" not in s:
            print(f"   {k:<3}{diag.F2D_VARIANTS[k][:40]:<42}{'too few finite values':>43}")
            continue
        print(f"   {k:<3}{diag.F2D_VARIANTS[k][:40]:<42}"
              f"{s['p50']:>12.3e}{s['p97']:>12.3e}{s['ratio']:>10.1f}"
              f"{s['frac_positive']:>9.1%}")

    print("\n   A signed, centred field shows frac>0 near 50 % and a p97/median")
    print("   ratio that is large and unstable — the median sits wherever the")
    print("   zero crossing happens to fall. A one-sided field shows frac>0 at")
    print("   100 % and a stable, modest ratio. The published pair describes")
    print("   the second kind.")

    # ---------------------------------------------------------- tail overlap
    print("\n" + "=" * 84)
    print(f"TAIL OVERLAP — which cells each variant flags at p{MOG_PERCENTILE} (MOG)")
    print("=" * 84)
    masks = {}
    for k in "ABCD":
        s = stats[k]
        if "pmog" not in s:
            continue
        masks[k] = np.isfinite(vals[k]) & (vals[k] >= s["pmog"])
    keys = list(masks)
    print("   " + " " * 6 + "".join(f"{k:>10}" for k in keys))
    for a_k in keys:
        row = []
        for b_k in keys:
            inter = np.count_nonzero(masks[a_k] & masks[b_k])
            union = np.count_nonzero(masks[a_k] | masks[b_k])
            row.append(inter / union if union else float("nan"))
        print(f"   {a_k:<6}" + "".join(f"{v:>10.3f}" for v in row))
    print("\n   Jaccard overlap of the flagged sets. A vs B near 0.000 is the")
    print("   expected result and is the whole point: the sign choice selects")
    print("   two disjoint populations of grid cells, so it cannot be left to a")
    print("   default. A vs C and B vs C near 0.5 each is also expected, since")
    print("   C is the union of both tails.")

    # ------------------------------------------------------------- ranks A/D
    rho = spearman(vals["A"], vals["D"])
    print("\n" + "-" * 84)
    print(f"   Spearman rho(A, D) = {rho:.4f}")
    print("   A and D differ by the |dv/dtheta|^-1 normalisation. Near 1.0 would")
    print("   mean the normalisation barely reorders anything and the question")
    print("   is only about sign; well below 1.0 means it is a real rank change")
    print("   and FORMULA_AUDIT.md §4.1's category-1 concern is live.")

    # ------------------------------------------------------------- verdict
    print("\n" + "=" * 84)
    scored = [(k, abs(np.log(stats[k]["ratio"] / PUBLISHED_RATIO)))
              for k in "ABCD"
              if "ratio" in stats[k] and np.isfinite(stats[k]["ratio"])
              and stats[k]["ratio"] > 0]
    if scored:
        scored.sort(key=lambda t: t[1])
        best = scored[0][0]
        print(f"   CLOSEST TO THE PUBLISHED SHAPE: variant {best}  "
              f"({diag.F2D_VARIANTS[best]})")
        print(f"   ranking: " + ", ".join(f"{k} ({stats[k]['ratio']:.1f})"
                                          for k, _ in scored))
        print()
        print("   This is one line of evidence, not a decision. It says which")
        print("   variant is distributed like the published one. Confirm with")
        print("   the sign of the trend against Prosser Figure S5 before")
        print("   changing --f2d-variant in the production run: re-run the")
        print("   NA months and the calibration under the candidate, then")
        print("   ada/trend_check.py --per-diagnostic.")
    else:
        print("   No variant produced a usable ratio — check the input file.")
    print(f"\n   PEAK RSS {peak_rss_gb():.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
