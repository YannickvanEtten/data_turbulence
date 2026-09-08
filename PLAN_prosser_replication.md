# PLAN — full replication of Prosser (2023), figure by figure

Written 2026-09-08, from the paper and the Supporting Information directly
(`Turbulence project/Articles/`). Companion to `STATUS.md` §15–16,
`RESULT_full_trend_42yr.md` and `CALIBRATION_REFERENCE.md` §4.4.

The question this answers: *what would it take to reproduce **everything**
Prosser publishes — every table, figure and appendix panel — so that a
disagreement can be traced to a specific diagnostic rather than absorbed into
an ensemble mean?*

**Short answer.** The paper is two objects. The North Atlantic box analysis is
already complete and its figures cost nothing but plotting. The global map
analysis is a separate 42-year global run: not blocked by storage, blocked by
about two weeks of CDS queue and compute.

---

## 1. The complete inventory

Everything Prosser publishes, with what it needs and where we stand.

| Item | What it is | Domain | Status |
|---|---|---|---|
| **Table 1** | Fitted 1979→2020 hours, 5 seasons × 5 severities | NA box | ✅ **replicated** — 25/25 significant, 25/25 containing his value (`RESULT_full_trend_42yr.md`) |
| **Figure 3a** | Annual diagnostic-mean MOG probability, 42 points + fit | NA box | ✅ data on disk (`per_diagnostic_annual_moderate_series.csv`, `ensemble` column) |
| **Figure 4** | The same fit for each of the 21 diagnostics | NA box | ✅ data on disk (same file, 21 diagnostic columns) |
| **Figure S3-LOG/SOG** | Figure 3 at LOG and SOG | NA + USA | ◑ NA-SOG on disk; **NA-LOG needs one `jobs/21` run**; USA panel needs a download |
| **Figure S4-LOG/SOG** | Figure 4 at LOG and SOG | NA box | ◑ SOG on disk; **LOG needs one `jobs/21` run** |
| **Figure 3b** | Figure 3 for the **USA box** (30–55°N, 124–60°W) | USA box | ❌ not downloaded |
| **Figure 1a–d** | **Global maps** of annual-mean MOG probability: 1979, 2020, and both inferred from the regression | global, 42 y | ❌ |
| **Figure 2a,b** | **Global maps** of absolute and relative change, Wald-test stippling at p=0.05 | global, 42 y | ❌ |
| **Figure S1-LOG/SOG** | Figure 1 at LOG and SOG | global, 42 y | ❌ |
| **Figure S2-LOG/SOG** | Figure 2 at LOG and SOG | global, 42 y | ❌ |
| **Figure S5** | Figure 2a broken down by all 21 diagnostics | global, 42 y | ❌ |

### 1.1 A precision caveat that applies to every figure row

**Table 1 is the only thing Prosser tabulates.** His Supporting Information
contains Figures S1–S5 and no tables at all (`CALIBRATION_REFERENCE.md` §3).
So every figure comparison is a value read off an image: a sign-and-magnitude
check, not a precision one. The scorecard records this per row rather than
letting exact-looking numbers imply exact comparisons.

Two of his figure-level claims *are* stated numerically in the main text and
can be compared exactly:

> "…significant (p < 0.05) upward trends, with relative changes of up to
> 75.6%. The remaining 4 diagnostics show no significant trend, and none of the
> diagnostics shows a significant downward trend."

That is **17 of 21 significant, maximum +75.6 %, zero significant decreases** —
three hard numbers for Figure 4, which is exactly the per-diagnostic check this
project wants. `ada/prosser_scorecard.py` compares against them directly.

---

## 2. Why the global half is expensive, and why it is not a storage problem

The instinct is that global × 42 years is impossible on a 2.5 TB share. Run the
numbers and the constraint moves.

**Storage is not the obstacle.** The maps need only an *annual* exceedance
field per gridpoint. At 721 × 1440 float32 that is **4 MB per year per
severity** — even keeping all 21 diagnostics separately at all 5 severities is
~17 GB for the whole 42 years. The 3-hourly derived field is an intermediate
that can be deleted as soon as the year is aggregated. Transient footprint is
~324 GB per year in flight (128 GB raw + 196 GB derived, both measured in
§16.2), which fits several times over in the 1.5 TB free.

**Time is the obstacle.**

| Stage | Rate (measured, §15.7) | 42 years |
|---|---|---|
| CDS download, 3 blocks per global month, `%4` | 5–21 min per block | **~1,500 requests, 3–4 days** |
| Global diagnostics, `%8` (QOS cap) | 3:18–5:01 per month | **~504 months, ~10 days** |
| Annual aggregation + per-gridpoint OLS | minutes per year | negligible |

Roughly **two weeks of mostly unattended wall clock.** This is the one scenario
that would make the `unlimited` QOS worth asking for after all (§16.3 dropped
it precisely because no such job was planned) — it would roughly halve the
compute half. The download half is bounded by Copernicus and would not move.

**Two shortcuts that look tempting and must not be taken.** Coarsening the grid
or dropping to 6-hourly would cut this substantially, and both change the
gradients the thresholds were calibrated on (`PLAN_full_year_calibration.md`
§9.1). The thresholds must be applied to data sampled exactly as the
calibration was.

---

## 3. Tiers, in the order they earn their cost

**Tier 1 — the North Atlantic figures. No new data. CHOSEN 2026-09-08.**
Figures 3a, 4, S3-a and S4 come from CSVs already on disk, plus one `jobs/21`
run for the LOG severity. Days, not weeks. This is where the "is one of the 21
diagnostics broken" question actually gets answered, because Figure 4 is
per-diagnostic and Prosser states its summary statistics in the text.

**Tier 2 — the USA box.** ~470 GB and ~1.5 days, giving Figure 3b and the USA
halves of S3. Worth more than the missing panels suggest: it is a genuinely
independent second region, so it tests the whole chain rather than adding
another view of the same data.

**Tier 3 — the global maps.** Figures 1, 2, S1, S2, S5. Two weeks, rolling
per-year storage. Do this only if Tier 1 leaves a question that the maps would
settle — and note that **S5 is the only published per-diagnostic map in
existence**, so if a diagnostic turns out to be spatially pathological, this is
the figure that shows it.

---

## 4. The mechanism: reduce on ADA, render anywhere

```
reduce   ADA, hours, jobs/16 and jobs/21   ->  cat_outputs/*.csv   (~100 KB)
render   anywhere, seconds, one script     ->  cat_outputs/figures/*.png
score    anywhere, seconds, one script     ->  cat_outputs/SCORECARD.md
```

**The CSVs are the archived result; the figures are disposable.** A plot can
then be iterated on in seconds on Windows instead of in a queue, and the whole
figure set regenerates identically months later from one command.

A Jupyter session on ADA (via Open OnDemand — the `ood` QOS in §11.9 implies it
exists) is the right tool for *looking* at a field when a diagnostic
misbehaves, because that needs the 3-hourly data that never leaves the cluster.
But a finding from a notebook gets written into a script or into `STATUS.md`;
the notebook is never the record. Cluster notebooks are hard to version, easy
to leave half-run, and impossible to re-execute reliably a month later.

### 4.1 New files

```
ada/prosser_figures.py     renders Figures 3a, 4, S3-a, S4 from the CSVs.
                           Recomputes each fit from the series and CROSS-CHECKS
                           it against the archived fit CSV, failing loudly if
                           the picture and the table disagree by >0.5 pp.
ada/prosser_scorecard.py   one row per published item: ours, Prosser's, ratio,
                           the precision of that particular comparison, and a
                           verdict. Parses the 25 Table 1 cells out of the
                           jobs/16 log rather than re-running a 6 h 20 m job.
```

Both need only numpy (+ matplotlib for the figures) and neither imports the
pipeline, so they run on ADA or on Windows.

### 4.2 Running it

```bash
# 1. on ADA — the one missing reduction (LOG at annual, for S3/S4-LOG)
sbatch --export=ALL,SEASON=annual,SEVERITY=light jobs/21_per_diagnostic_trend.sbatch

# 2. on ADA — make the small CSVs travel
git add cat_outputs/*.csv && git commit -m "per-diagnostic reductions" && git push

# 3. on Windows — pull, render, score
git pull
python ada/prosser_figures.py
python ada/prosser_scorecard.py --log logs/full-trend-1101414.out
```

Step 3's `--log` is the jobs/16 log that produced
`RESULT_full_trend_42yr.md`; copy it into the repo's `logs/` first, or point
`--log` at it on ADA and run the scorecard there.

**`matplotlib` is not declared in `pixi.toml`.** Rendering on ADA needs
`matplotlib-base = "*"` added to it, which forces a pixi re-solve. Rendering on
Windows needs nothing — which is the intended workflow anyway, and the reason
the split above exists.

---

## 5. What Tier 1 will and will not settle

**Will settle.** Whether our 21 diagnostics individually behave like Prosser's:
his stated 17-of-21 significant, +75.6 % maximum, and zero significant
decreases are three exact numbers to check against. `horizontal_divergence` is
already the standing suspect (§15.6a — no trend in 0 of 7 fits, negative in two
seasons), and Figure 4's panel for it is the published comparison that either
convicts or clears it.

**Will not settle.** Anything spatial. If a diagnostic is wrong in a
latitude-structured way — the known failure mode of `ubf` (§4g), and the
reason §9 of `PLAN_full_year_calibration.md` argued for keeping the full global
derived field — no North Atlantic box average will show it. That is Figure S5,
and it is Tier 3.

---

## 6. Record-keeping

This file is the plan. When Tier 1 produces its figures and scorecard, the
findings go into `STATUS.md` as a new section — not into this file, and not
into a notebook. `SCORECARD.md` is regenerated, never hand-edited.
