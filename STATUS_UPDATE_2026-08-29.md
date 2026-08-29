# STATUS update — 2026-08-28/29

Paste §11 below into `STATUS.md` after §10. Then apply the corrections in
"Edits to existing sections" at the end — several earlier statements are now
wrong, and leaving them is worse than not writing this at all.

Detail on the literature, the reference tables and the full result is in
`CALIBRATION_REFERENCE.md`. This file is the summary and the decisions.

---

## 11. Stage 2 executed — the calibration check, 2026-08-28/29

**The headline: the chain reproduces Prosser (2023) Table 1 to within 24 %,
and the identity check on real data is exact.** The dataset does not exist
yet, but the method that will produce it has now been validated end to end
against a published result.

### 11.1 What was run

| Step | Job | Cost |
|---|---|---|
| Global year 2000, days 1/9/17/25 of each month | 1092393 (12 tasks) | ~35 min wall, 16 GB |
| North Atlantic Feb + Dec 1979 | 1092394 | ~10 min, 0.7 GB |
| Diagnostics, 12 global months | 1092536 | 15–36 min each, 27 GB out |
| Diagnostics, NA Jan/Feb/Dec 1979 | 1092537 | 5–8 min each |
| Calibration check | 1092630 | 117 min |

**14 CDS requests and about four hours of compute**, against 504 requests and
days for the production run. All 15 input files passed integrity checks
(`ada/check_calib_downloads.py`); all 15 diagnostics runs produced 21 of 21
with no failures and no all-NaN placeholders.

### 11.2 The result

North Atlantic DJF 1979, in **Prosser's own box (36–60°N / 55–10°W)**:

| | ours | Prosser | ratio |
|---|---|---|---|
| LOG | 5.037 % | 5.968 % | 0.84 |
| LMOG | 1.684 % | 2.111 % | 0.80 |
| **MOG** | **0.785 %** | **1.032 %** | **0.76** |
| MSOG | 0.399 % | 0.560 % | 0.71 |
| SOG | 0.200 % | 0.296 % | 0.68 |

Order of magnitude **PASS**. Ladder shape, worst deviation 11 % — **PASS**.
Absolute agreement 24 % — **PASS**, and this one was not expected to be met.

**The identity check came back exact**: applying the year-2000 thresholds back
to the year-2000 data returns 3.0000 / 0.9000 / 0.4000 / 0.2000 / 0.1000 %,
0.00 % relative error on all five. That is what licenses reading the table
above; without it, agreement with Prosser would be coincidence.

### 11.3 The one systematic pattern, and its cause

The ratios decline monotonically (0.84 → 0.68): our tail is thinner than
Prosser's, more so the further out you go.

**The vertical stencil predicts exactly this.** We differentiate across
175→225 hPa (50 hPa); Prosser uses 188→206 (18 hPa). A 2.8× wider stencil
smooths across sharp shear layers and under-resolves shear *maxima* — barely
visible in the bulk, progressively damaging in the extreme tail.

This is §5.4's documented level substitution, finally measured. And per §11.6
it cannot be closed, so it should be quoted as a known property of the dataset
rather than treated as an open defect.

### 11.4 Method facts now pinned

- **The percentile ladder is LOCKED**: LOG p97.0, LMOG p99.1, MOG p99.6,
  MSOG p99.8, SOG p99.9. Williams (2017) Table 1, adopted verbatim by Prosser.
  Origin: a log-normal EDR distribution constrained by an observed 3.0 % LOG
  and 0.4 % MOG probability.
- **Prosser recomputes thresholds from ERA5; he does not reuse Williams'
  numbers.** This was genuinely open and is now settled from the paper's own
  text. Latitude weighting is required by the method, not an embellishment.
- **Prosser's North Atlantic box is 36–60°N / 55–10°W** (his Figure 2 caption).
  It is a strict subset of our 30–60°N / 75–0°W download, so Table 1 can be
  compared exactly. **This coordinate pair was recorded nowhere before.**
  Comparing on our larger box gives systematically different frequencies and
  looks like an error.
- **Prosser publishes no threshold values.** Confirmed by reading his
  Supporting Information directly: Figures S1–S5, no tables. The only published
  per-diagnostic threshold table anywhere is Williams (2017) Table 2, and it is
  GFDL-CM2.1 — magnitudes are not comparable, as Williams states himself.

### 11.5 Bugs found and fixed

1. **Colson–Panofsky was compared in the wrong units.**
   `3_pipeline.SCALE_TO_TABLE["colson_panofsky"]` was `1.0` with a comment
   promising a "see note below" that was never written, so a native value in
   m² s⁻² was compared against a table printed in **10³ kt²**. That single line
   is the entire reason CP has read `1-2 ORDERS OFF` since the comparison table
   was built. Correct factor 1/(1000 × 0.2646526) = **3.778538 × 10⁻³**.
   CP is the only one of the 21 not published in SI.

2. **`calibration.compute_thresholds` sorts each array five times.** It asks
   `weighted_percentile` for one severity at a time, and that function sorts
   its whole input on every call — 105 full sorts of 4 × 10⁸ values where 21
   would do. The first attempt at the check would not have finished inside any
   sane walltime. `ada/calibration_check.py` memoises the five severities onto
   one sort via the injected-function seam; measured 5.67× on a benchmark,
   bit-identical output. **The library itself is still doing 5× the work for
   any other caller** — worth fixing in `calibration.py` before phase 4.

3. **A unit error in the new comparison code itself** (mine): Williams prints
   NVA/RVA in 10⁻⁹ s⁻² where `REFERENCE_TABLE` uses 10⁻¹⁰; converting requires
   multiplying by ten and the first version divided, inflating those two ratios
   by a hundred and making them look like wild outliers.

### 11.6 The 188/197/206 levels are MODEL levels — the branch is closed

Checked against ECMWF's L137 definitions: **188.29, 197.37 and 206.81 hPa are
model levels 73, 74 and 75.** Prosser used ERA5 model-level data described by
nominal pressures, not the pressure-level product.

- They **cannot** be requested from `reanalysis-era5-pressure-levels`, whose
  fixed 37-level set offers 150, 175, 200, 225, 250 around 200 hPa and nothing
  between.
- Model-level ERA5 is on the **MARS tape archive** (`reanalysis-era5-complete`)
  — different request format, slow access. Same obstacle as ERA5.1.
- **Therefore 175/200/225 is not a compromise: it is the finest stencil
  available at 200 hPa on CDS.** §6 phase 3's "accept 175/200/225, or pursue
  model levels" now has a price attached, and the answer is accept.

**What can be done instead** (§11.1 of `CALIBRATION_REFERENCE.md`): request
150/175/200/225/250 in one download and compute at 200 hPa with both a 50 hPa
and a 100 hPa stencil. Two points on the stencil-width curve let the 18 hPa
case be *extrapolated* rather than measured. ~1 day. Worth doing before
phase 4, because "why are your levels different from Prosser's" is a certain
referee question and this converts the answer from an apology into a number.

### 11.7 ERA5 vs ERA5.1 — also MARS-only

Prosser used **ERA5.1** for 2000–2006, which corrects a lower-stratospheric
cold bias — and year 2000 is our calibration year. ERA5.1 is **not a CDS disk
dataset**; it lives in MARS. So this is not a one-constant change either.
ECMWF's own note is that behaviour "in most of the troposphere is similar to
that in ERA5.1". Proceeding with ERA5; revisit at publication, not before.

### 11.8 The sub-sampled calibration design

Days 1/9/17/25 of each month, all 8 three-hourly steps: **48 days,
3.99 × 10⁸ points per diagnostic, 1.59 × 10⁶ above the MOG threshold.**

- Split-half (days 1+17 vs 9+25): median disagreement **0.81 %**, 90th
  percentile **5.90 %**. **48 days is enough.** The script's "marginal" verdict
  came from a 92 % worst case on `colson_panofsky @ LMOG`, which is an artifact
  — CP's LMOG threshold sits at −1.04, essentially zero, and a relative
  difference across a sign change is meaningless. The metric should divide by
  the diagnostic's spread, not its value.
- It also **made stage 2 possible at all**: a full global year is 3.04 × 10⁹
  points per diagnostic, and `compute_thresholds` ravels the whole field —
  over 100 GB peak. It does not fit any ADA partition.
- And it sidestepped the unconfirmed CDS per-request field limit: 672 fields
  per request instead of 5,208.

### 11.9 ADA facts learned the hard way

- **`defq` is heterogeneous.** 59 GB (003, 004, 220–222), 124 GB (007, 008,
  016, 017), 252 GB (009–011, 013–015), 1031 GB (001, 002). The fat and GPU
  nodes are *also* in `defq`.
- **Memory footprint is the biggest lever on time-to-start.** A 120 GB request
  was scheduled **22 hours out** while 24 GB started in 13 seconds. Do not
  round requests up "to be safe" — it costs a day.
- **Do not submit to `defq-thin` or `defq-fat` directly.** Their nodes are
  shared with `defq`, which has higher priority, so a job addressed to the
  narrower partition queues behind everything. Ask `defq` for the memory you
  need and let the scheduler pick.
- **Job-step accounting is DISABLED.** `sacct` and `sstat` both return empty
  MaxRSS. Every driver must measure and print its own peak RSS
  (`resource.getrusage`). **This means §10.12's 14.8 GB figure cannot be
  reproduced by the method that produced it** — though the global run's 17.7 GB
  on a comparable point count corroborates it.
- **`export PYTHONUNBUFFERED=1` in every job script.** Python block-buffers
  stdout to a file, so a working job is indistinguishable from a hung one.
  This cost an hour of watching a blank log.
- **`mkdir -p logs` before `sbatch`, not inside the script.** SLURM opens the
  output file before the script runs.
- **Short, honest walltimes get backfilled.** 6 h waits; 1 h runs.
- **HARD CAP: 8 CONCURRENTLY RUNNING JOBS PER USER.** The default `normal`
  QOS sets `MaxJobsPU = 8`; confirmed both in `sacctmgr show qos` and in the
  association manager (`yen230 ... MaxJobsPU=8`). Array tasks count as jobs,
  so a 9th sits as `QOSMaxJobsPerUserLimit`.

  **This makes any `%` throttle above 8 a no-op**, and §10.12's plan to run
  the 504-month diagnostics array at `%16` unreachable as written.

  **There is an `unlimited` QOS on this cluster** with no `MaxJobsPU`, and an
  `ood` QOS capped at 2. Whether `yen230` may use `unlimited` is a question
  for `itvo.ucit@vu.nl`, and it is the single highest-leverage thing to ask
  before phase 4 — it is the difference between 8-way and node-limited
  concurrency across 504 months. Check what is currently granted with:

      sacctmgr show assoc user=yen230 format=User,Account,Partition,QOS,DefQOS

### 11.10 Measured costs, for sizing phase 4

| Quantity | Measured |
|---|---|
| Global month, 4 days sub-sampled, download | 78 s – 23 min (CDS cache vs tape) |
| Global month, 21 diagnostics, `--chunk-days 1` | 30 min, **17.7 GB peak**, 2.28 GB zarr |
| NA month, 21 diagnostics, whole file | 5–8 min, 12.9–14.3 GB peak |
| Calibration check, 12 global months | 117 min, 32 GB |
| Diagnostics under 6-way concurrency | 15–36 min (vs 30 min alone) |

**§10.12's phase-4 timing is wrong, for two independent reasons.** It plans
the 504-month diagnostics array at `%16` and projects 4.3 h, reasoning purely
from memory per node. But:

1. **The QOS caps concurrency at 8** (§11.9), so `%16` cannot happen.
2. **Shared-filesystem I/O contends.** Six concurrent tasks already stretched
   a 30-minute job to 36, mostly in GRIB decode. Mild — ~20 % on the worst
   task — but real, and it grows with concurrency.

Corrected estimates at the 8-job ceiling:

| Stage | 504 months at 8-way |
|---|---|
| Download (measured 7–23 min each, CDS-dependent) | **~8–16 h** |
| Diagnostics (measured 5–8 min each) | **~7 h** |

So roughly **a day, possibly spread over two** — comfortably inside the 7-day
walltime cap, but not the 4.3 h §10.12 implies. If the `unlimited` QOS can be
granted, the diagnostics stage shrinks toward the I/O ceiling instead; the
download stage will not, because it is bounded by Copernicus, not by ADA.

### 11.11 New files

```
ada/calibration_check.py        stage 2->3: calibrate, verify, compare to Prosser
ada/diagnostics_global.py       21 diagnostics, overlap-buffered time chunking
ada/check_calib_downloads.py    verify all 15 inputs in one pass
ada/progress.sh                 one line per job; reads logs, liveness from squeue
jobs/02b_download_global_calib.sbatch
jobs/04_diagnostics_global.sbatch
jobs/05_diagnostics_na_djf.sbatch
jobs/06_calibration_check.sbatch
tests/test_calibration_roundtrip.py   31 tests, the identity, no ERA5 needed
RUNBOOK_calibration_check.md
calibration/thresholds_2026-08-29.json
```

`--days` added to `1_download_hpc.py`, with a distinct filename
(`era5_glob_2000-01_d01-09-17-25.grib`) so a partial month can never be
mistaken for a complete one.

**`ada/diagnostics_global.py`'s chunking is bit-identical to whole-file
processing** — verified across 1-, 2-, 4-day and whole-file chunk sizes. Each
chunk is computed with one extra timestep on each side and trimmed back, so
F2D loses only the 2 steps at a file's true ends. Naive chunking would cost
it 8. Same device as `chunk_stitch.py`, one level down.

---

## Edits to existing sections

- **§5.4** — add: the level substitution's damping is now measured, and it is
  *not* constant across the distribution. It is ~16 % at LOG and ~32 % at SOG
  (§11.3). The percentile argument holds for the bulk; the tail is where it
  leaks. This is the proviso §5.4 itself flagged, now with a number.
- **§6 phase 3, "Pressure levels"** — resolved. Model levels are MARS-only
  (§11.6); 175/200/225 is the finest available on CDS. Not a decision, a
  constraint.
- **§6 phase 3, "Calibration domain"** — resolved in favour of global, and
  demonstrated feasible via sub-sampling (§11.8).
- **§7** — the units half is closed: all 21 `REFERENCE_TABLE` units and
  `wj_median` values were checked line by line against both papers and all 21
  match. Colson–Panofsky's *scale factor* was wrong (§11.5) but its recorded
  unit was right. The negative-Ri branch coverage gap in §4a is also closed
  (§10.5 of `CALIBRATION_REFERENCE.md`). What remains open is the formula
  derivations against Sharman (2006) Appendix A and Ellrod & Knapp (1992) —
  and `brown2` is now a specific suspect, five orders from the published value
  in a way resolution does not explain.
- **§10.3** — add the per-node memory table (§11.9). "17 nodes" is true but
  useless without knowing that five of them cannot hold a 60 GB job.
- **§10.12** — flag that the 14.8 GB MaxRSS figure is not currently
  reproducible, since step accounting is off (§11.9), and that 17.7 GB from
  the global run corroborates it independently.
- **§10.9** — the ERA5 Terms of Use item is closed; the 1979-01 download on
  2026-08-27 could not have succeeded otherwise.
- **§3 "Not done"** — still true that no month of the *production* series has
  been through the full stack, but the stack itself has now been exercised end
  to end on 15 months and validated against Prosser. Worth saying, because the
  sentence as written understates where the project is.
