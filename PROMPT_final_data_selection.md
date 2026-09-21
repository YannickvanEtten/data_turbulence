# TASK — finalise the CAT indicator dataset

Close out the data-gathering phase. The outcome I need is **one fixed dataset**
I can build econometric research on with confidence, plus a written
justification a referee in atmospheric science would accept. Not a new
comparison, not another variant, not another week of runs.

---

## READ FIRST, IN THIS ORDER

These are project docs, not files on disk. Read them with `project_read`.

1. **`claude/DECISION_RECORD.md`** — the decision as it currently stands.
2. **`claude/FORMULAS_AND_DECISIONS.md`** — the 21 formulas with equation and
   page, the decision register D1–D11, the measured sensitivities, the units
   caveats, the primary-source list.
3. **`claude/STATUS.md` §18** — the audit and the re-derive. §18.12 through
   §18.17 are the current state; §1–17 are history and their forward-looking
   lines are superseded.
4. **`claude/RESULT_full_trend_42yr.md`** — the previous-conventions result.
5. `claude/AUDIT_diagnostics_vs_literature.md` — only if you need a specific
   diagnostic's individual audit.

**Treat these as the record.** Do not re-derive what they establish. Your job
is to audit the conclusion and close the remaining gaps, not to repeat three
weeks of work.

---

## WHAT EXISTS

On the ADA share `/scistor/SBE-EDS-ClimateKoopman/yen230/`:

```
raw/north_atlantic/          504 GRIB months, 1979-2020, 3-hourly, 7 variables, 3 levels
raw/global/                  the full contiguous year 2000 (+ 10-day day-blocks)
raw/global_sub48/            the superseded 48-day sub-sample, kept as a comparison set

derived/north_atlantic_audit/   THE CANDIDATE. 504 zarr, F1+F12+F6 conventions
derived/north_atlantic/         the previous-conventions series, 504 zarr
derived/global_audit/           reference year under the audit conventions
derived/global/                 reference year under the previous conventions

calibration/thresholds_2026-09-11.json   pairs with the *_audit series
calibration/thresholds_2026-09-07.json   pairs with the baseline series
calibration/tails_2026-09-11/            41 GB tail archive, GLOBAL REFERENCE YEAR ONLY

cat_outputs/                 scorecard_audit.csv, SCORECARD_audit.md,
                             per_diagnostic_annual_moderate_audit.csv, threshold_shift_*.csv
```

Repo: `data_turbulence/` on ADA and on Windows (OneDrive + GitHub).
Share is at ~1.5 TB of 2.5 TB. Storage is not a constraint.

---

## HARD CONSTRAINTS — these are not negotiable

1. **No new CDS or MARS downloads.** None. The raw data is complete.
2. **No new full derived series.** The two that exist are the two that will
   ever exist. Do not propose a third set of conventions.
3. **Total new cluster time ≤ 3 hours**, submitted as one batch, unattended.
   If a check cannot fit in that, it does not happen — say so and move on.
4. **Every check must end in a decision.** If a proposed check could come back
   ambiguous and leave a new open question, do not propose it. I have had
   enough of measurements that spawn more measurements.
5. **Prefer desk work to cluster work.** The most valuable findings in this
   project came from reading a table, a source file, and a documentation
   section — each about ten minutes — after weeks of jobs. Check what is
   already on disk before you queue anything.

---

## WHAT I NEED YOU TO PRODUCE

### 1. `FINAL_DATASET.md` — the selection

State, in one page, exactly which dataset I use and why, at the level of:
paths, conventions, calibration file, pairing rule. If you disagree with
`DECISION_RECORD.md`'s choice, say so with evidence and say what would change
your mind. If you agree, say why the evidence is sufficient rather than
restating it.

### 2. The honest limitations section

For every place this dataset does **not** reach Prosser or the literature:
what it is, **why** it could not be closed, and **what the measured impact
is**. I want the impact quantified, not described. Known items to cover, plus
anything you find:

- the vertical stencil (175/200/225 vs his model levels 73–75)
- ERA5 rather than ERA5.1 for 2000–2006
- `f2d`'s diurnal sampling damping
- `ncsu1`'s level ratio against his Figure 4
- `endlich`'s level ratio, and the unresolved reading of its vertical derivative
- `brown2`'s dimensional gap
- the four diagnostics that do not reach significance
- the global and USA-box figures that were never computed

For each: could it ever be closed, at what cost, and does it affect the trend,
the level, or only a magnitude comparison? **Separate "we chose not to" from
"it cannot be done" — they are different sentences in a paper.**

### 3. A bounded list of additional checks

Propose the checks that are worth running **given what is already on disk**.
For each, state: what decision it changes, what it costs, and what the two
possible outcomes are. Then run the ones that pass your own bar, in one batch.

Candidates I am aware of — evaluate them, do not assume they are all worth it:

- the per-diagnostic 1979 level ratios on the audit series (a join of
  `cat_outputs/per_diagnostic_annual_moderate_audit.csv` with
  `ada/prosser_published.py`; minutes, no cluster) — this rebuilds STATUS
  §17.5's investigate list on the series I will actually use, and it is the
  table that decides my per-indicator caveat list
- a provenance sweep of all 504 audit stores (`jobs/15` currently hardcodes the
  baseline paths; count, size and one store's attributes are verified, all 504
  are not)
- a silent-failure sweep of the 504 job logs for all-NaN placeholders
- per-diagnostic fits at `light` and `severe` (the scorecard lists Figures
  S4-LOG and S4-SOG as NOT RUN)
- per-diagnostic fits for DJF, to check whether any indicator's verdict is
  season-dependent
- regenerating the baseline `per_diagnostic_annual_moderate.csv`, which a
  mismatched run overwrote on 2026-09-15

### 4. The per-indicator verdict

A table: for each of the 21, can I use it without caveat, with a stated
caveat, or not at all. This is the thing I actually need for the econometrics.

---

## RULES OF EVIDENCE — learned expensively, please follow them

- **A magnitude comparison cannot settle a formula question.** Thresholds are
  percentiles of the data itself, so any constant multiplicative factor cancels
  exactly. Use exceedance-set flip rates, not ratios against published medians.
- **A series must be scored against thresholds built under the same
  conventions.** `full_trend_check.py` and `per_diagnostic_trend.py` enforce
  this and exit 2 on a mismatch. A mismatched run on 2026-09-15 produced a
  plausible-looking table in which one diagnostic silently had zero exceedance.
- **Verify the file on ADA, not the commit on Windows.** Files written into the
  Windows repo do not reach ADA until committed, pushed and pulled, and this
  has cost two jobs. `grep` the thing you changed on ADA before you `sbatch`.
- **Re-read the directory, not your note about the directory.** A paper was
  reported missing for a week while sitting in `Articles/`.
- **Check the output at 60 seconds, not at 90 minutes.** Every driver prints
  its configuration on the first lines.
- **Do not automate a verdict.** Report the statistic and an interval; a
  threshold-based pass/fail on a small sample has flipped on rounding twice in
  this project.

---

## SPECIFIC QUESTIONS TO ANSWER

1. Is `derived/north_atlantic_audit/` + `thresholds_2026-09-11.json` the right
   final choice? What is the strongest argument against it?
2. Which of the 21 indicators would you not hand to an econometrician without
   a flag, and exactly what is the flag?
3. `ncsu1` sat at a 6.07 level ratio against Prosser on the previous
   conventions and three candidate explanations were tested and failed. Where
   does it stand now, and is it usable?
4. Is there anything in the pipeline that is still **unverified** rather than
   verified-and-documented? I want to know what has been assumed.
5. What is the one thing most likely to be challenged by a referee, and what is
   the answer?

---

## WHAT I DO NOT WANT

A third convention set. A new download. A check whose outcome is "interesting,
needs another run". A restatement of the documents I have already read. An
optimistic summary — if something is weak, say it is weak and say how weak.

The deliverable is a dataset I stop thinking about, and a page that explains it.
