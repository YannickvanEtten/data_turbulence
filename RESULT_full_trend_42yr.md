# RESULT — the 42-year, all-season trend

`ada/full_trend_check.py` via `jobs/16_full_trend_check.sbatch`, run twice:
job `1098646` (2026-09-04, 48-day thresholds, 4 h 06 m) and job **`1101414`**
(2026-09-07, **full-year thresholds**, 6 h 20 m, exit 0), whose figures are the
ones reported here. Logs: `logs/full-trend-1098646.out`,
`logs/full-trend-1101414.out`. §5 compares the two.

> **PROVENANCE — RESOLVED 2026-09-07.**
> The figures below are computed on **`calibration/thresholds_2026-09-07.json`**,
> the FULL contiguous reference year (job `1101414`, 6 h 20 m, exit 0). An
> earlier version of this document reported the same run on the 48-day
> sub-sampled thresholds and carried a warning here. **That warning is
> discharged, and the difference between the two is now measured rather than
> assumed — see §5.** Short version: recalibrating on the full reference year
> moves the fitted trend by at most ONE percentage point in any of the 25
> season x severity cells, and changes no verdict.

This supersedes STATUS.md §12, whose numbers are a 9-season, DJF-only
subsample. **Two sections of STATUS.md are now wrong and are corrected in §3
below.** The first attempt (job `1098191`) hit a 4 h walltime at 40 of 42 years
and wrote nothing; `jobs/16` is now at 12 h.

---

## 1. The headline

**25 of 25 fitted trends are significant at 5 %. 25 of 25 contain Prosser's
published value inside their 95 % interval.**

| season | severity | ours | Prosser | ratio | t |
|---|---|---|---|---|---|
| **Annual** | LOG | 16 % | 17 % | **0.96** | 4.08 |
| | LMOG | 27 % | 28 % | **0.96** | 4.86 |
| | MOG | 36 % | 37 % | **0.98** | 5.36 |
| | MSOG | 45 % | 46 % | **0.98** | 5.67 |
| | SOG | 55 % | 55 % | **0.99** | 5.84 |

The annual relative change is reproduced to **within 0–4 %** at every severity.
The seasonal pattern reproduces too — MOG trend by season, ours against
Prosser: DJF 34/37, MAM 60/57, JJA 31/31, SON 30/31, all within ±3
percentage points, with MAM correctly the strongest season in both.

### What changed relative to the subsample, and why it matters

| | n = 9, DJF only (§12) | n = 42, all seasons |
|---|---|---|
| significant at 5 % | **2 of 5** | **25 of 25** |
| Prosser inside 95 % CI | 5 of 5 | 25 of 25 |
| MOG interval | [−3 %, +89 %] | [+23 %, +50 %] |

**STATUS.md §12.3's caveat is now discharged.** It said, correctly at the time,
that *"what has been established is consistency with Prosser, not independent
confirmation of a trend"*, because the intervals were wide enough to contain
both Prosser's value and zero. At n = 42 the MOG interval is [+23 %, +50 %] and
excludes zero comfortably. This is now an **independent confirmation**, and it
should be written that way.

---

## 2. The level deficit is real, stable, and does not drift

Annual, ours ÷ Prosser, in hours per year at an average point in his box:

| | level 1979 | level 2020 | drift | trend ratio |
|---|---|---|---|---|
| LOG | 0.881 | 0.876 | −0.005 | 0.96 |
| LMOG | 0.857 | 0.851 | −0.006 | 0.96 |
| MOG | 0.841 | 0.837 | −0.005 | 0.98 |
| MSOG | 0.820 | 0.817 | −0.003 | 0.98 |
| SOG | 0.797 | 0.796 | −0.001 | 0.99 |

Two things to read here.

**The deficit still deepens with severity** — 0.88 at LOG to 0.80 at SOG. That
is the signature §11.3 attributes to the 50 hPa vertical stencil against
Prosser's 18 hPa: a wider stencil under-resolves sharp shear maxima, and the
damage grows toward the extreme tail. **That explanation survives.**

**The deficit does not move across the record.** Drift is −0.005 to +0.006 over
41 years — flat to three decimal places. So the bias is a stable offset on the
*level*, and it cancels almost exactly out of the *change*. That is why the
trend ratios sit at 0.96–0.99 while the level ratios sit at 0.80–0.88.

**This is the best possible shape for the result.** A bias that is stable in
time is a bias that does not contaminate a trend estimate, which is the only
quantity the project actually claims.

---

## 3. Two corrections to STATUS.md

### 3.1 §12.5 — "level deficit and trend excess are one fact" — DOES NOT HOLD

§12.5 reported that the trend ratio rose monotonically **above** 1 as the level
ratio fell below it, products landing at 0.89–0.98, and concluded that this was
"one parameter seen twice": sitting lower on the tail makes the same physical
intensification produce a *larger* relative change, because exceedance
probability is convex in the shift.

On the full record the trend excess is not there. DJF, like for like:

| | trend ratio n = 9 | trend ratio n = 42 | shift |
|---|---|---|---|
| LOG | 1.03 | 0.84 | −0.19 |
| LMOG | 1.08 | 0.90 | −0.18 |
| MOG | 1.16 | 0.92 | −0.24 |
| MSOG | 1.29 | 0.95 | −0.34 |
| SOG | 1.43 | 0.98 | −0.45 |

Every ratio falls, and the shift is largest exactly where §12.5's argument was
strongest. The monotone *increase with severity* survives (0.84 → 0.98), but it
now runs up toward 1 from below rather than away from 1 from above, and there
is no excess to pair against the deficit. The annual products are 0.85, 0.81,
0.81, 0.78, 0.76 — declining, not the near-constant 0.89–0.98 that made the
"one parameter" reading attractive.

**§12.5 was a subsample artifact, and §12.4 already predicted it.** §12.4
recorded that the nine-season fit leaned on two strongly positive-NAO winters
(2015 and 2020, both ~1.29 % MOG against a 0.785 % start), and that across nine
winters the MOG frequency spanned a factor of 1.65 against a fitted change of
43 %. A fit anchored on two high winters at the end of the sample overstates
the trend, and overstates it most at the severities where the tail is thinnest.
That is precisely the pattern in the table above. **§12.4 was right and §12.5
is its casualty** — which is a good outcome for the project's own reasoning,
because the mechanism that invalidated §12.5 was already written down.

Also to correct in §12.5: *"our agreement with Prosser improves across the
record — the MOG level ratio goes from 0.79 in 1979 to 0.83 in 2020"*. On the
full record the MOG level ratio goes from **0.823 to 0.820**. It does not
improve; it is flat.

### 3.2 §5 point 5 — the amendment should itself be amended

§5 point 5 originally argued that a constant magnitude offset cancels out of a
percentile-calibrated exceedance field. On 2026-08-29 it was amended:
*"the proviso has been measured, and it does not hold … the percentile argument
holds for the bulk and leaks in the tail."*

That amendment rested on §12.5. With the trend ratios at 0.96–1.00 annually and
the level ratio drifting by less than 0.006 over 41 years, **the offset does
behave as a stable multiplicative bias for the purpose of a trend**. It does
*not* cancel out of the level — the deficit is real, monotone in severity, and
must still be reported. The correct statement is narrower than either version:

> The stencil bias is stable in time. It biases the exceedance **level** by
> 12–20 %, deepening with severity, and it cancels out of the **relative
> change** to within 0–4 %.

---

## 4. What this settles, and what it does not

**Settled.** The pipeline reproduces Prosser (2023) Table 1 on the quantity the
paper is about — the change — across all four seasons and the annual figure,
at every severity, with independently significant trends. Every gate the
project set for itself now reads clean.

**Not settled, and not addressed by this run.**

- The **level deficit** (§2). It is characterised, stable and explained, but it
  is still a 12–20 % deficit and cannot be closed on CDS (§11.6).
- ~~**`ubf`** (§7, §12.6)~~ — **CHECKED 2026-09-07 and resolved.**
  `ada/per_diagnostic_trend.py` refit all 21 over 42 years with intervals.
  20 of 21 fitted trends came in BELOW their §12.6 endpoint numbers and `ubf`
  is the only one that ROSE (+11 % → +18 % DJF, +20 % Annual). §12.6's two
  claims are both false: `ubf` is rank 3/21 (DJF) and 2/21 (Annual), not last,
  and the siblings span +5 % to +83 %, not +38 % to +226 %. The consistent low
  outlier is **`horizontal_divergence`** (+5 % annual, +10 % DJF, not
  significant in either), not `ubf`. See STATUS §15.6. What remains open:
  the per-diagnostic breakdown may move just as much, and `ubf` may no longer
  be anomalous at all. **This needs no new data**, only the production series
  that already exists.
- The **21 sign entries** and **Brown (1973)** (§7). Unaffected by any of this
  and still the highest-value open items.

---

## 5. The full-year recalibration, MEASURED

This section previously carried a *prediction*. Job `1101414` (2026-09-07,
6 h 20 m, exit 0) re-ran the whole fit on `thresholds_2026-09-07.json`, so the
prediction and the outcome can be compared. Both are kept, because a forecast
that was written down before the run is worth more than one reconstructed after.

**The method.** `ada/compare_thresholds.py` diffs two threshold sets and
converts the difference into implied exceedance-frequency change, using each
diagnostic's own severity ladder to estimate the local elasticity (central
difference on the raw scale, so it also works for the three diagnostics whose
thresholds are negative). Measured moves: threshold **0.39 % median / 2.27 %
p90**, per-diagnostic frequency **1.94 % / 8.36 %**, ensemble **0.42 % / 1.83 %**.

**Predicted vs actual — annual 1979 exceedance level:**

| | predicted | actual |
|---|---|---|
| LOG | −0.00 % | **−0.06 %** |
| LMOG | +1.24 % | **+1.01 %** |
| MOG | +1.68 % | **+2.13 %** |
| MSOG | +2.17 % | **+3.43 %** |
| SOG | +4.60 % | **+4.58 %** |

**The level deficit shrank at high severity, as predicted:**

| annual level ratio | 48-day | predicted | full-year |
|---|---|---|---|
| LOG | 0.881 | 0.881 | **0.881** |
| MOG | 0.823 | 0.837 | **0.841** |
| SOG | 0.757 | 0.792 | **0.797** |

**And the trend did not move.** Across all 25 season × severity cells:
**16 identical, 9 changed by exactly one percentage point, none by more.**
Annual trend ratios went 0.96 / 0.96 / 0.98 / 0.99 / 1.00 → 0.96 / 0.96 /
0.98 / 0.98 / 0.99 — a maximum shift of 0.01. Significance stayed 25/25 and
Prosser-inside-CI stayed 25/25. **No verdict anywhere changed.**

**This is the level/trend separation of §2, confirmed quantitatively rather
than argued.** A threshold shift scales exceedance by nearly the same factor in
1979 and in 2020, and `P₂₀₂₀/P₁₉₇₉` is exactly invariant to a common factor;
only the difference in local tail slope between the two years survives, and it
is worth ≤ 1 percentage point.

**The sentence this buys for a write-up:** *recalibrating the reference year
from a 48-day sub-sample to the full contiguous year changes the fitted trend
by at most one percentage point at any season or severity, and changes no
significance verdict.* That is a measured sensitivity, not an assurance, and it
is the direct answer to the obvious referee question about the sub-sample.

**One caveat on the level, which did move.** The sub-sample was not merely
noisier — it was **systematically different**, and the difference grows with
severity (ensemble mean Δ frequency: −0.0 % at LOG rising to +4.6 % at SOG,
with 16 of 21 diagnostics moving the same direction). The two that move
opposite, `colson_panofsky` and `negative_richardson`, are both functions of
Richardson number, and `vertical_wind_shear`'s threshold moved −0.09 %, i.e.
nothing — so with Ri = N²/Sv² and the shear term static, the difference is in
**N², the static stability**. Forty-eight days spread eight apart is ~50
independent synoptic states, and the effective sample for an extreme quantile
of a spatially correlated field is set by the number of weather realisations,
not the number of gridpoints. Reported levels should come from the full-year
calibration for that reason, independent of the precision argument.

---

## 6. Run status

| | |
|---|---|
| `jobs/16` full trend, 48-day thresholds | job `1098646`, 4 h 06 m, exit 0 |
| `jobs/16` full trend, **full-year thresholds** | job `1101414`, **6 h 20 m, exit 0** — the figures in this document |
| `jobs/17`–`19` full-year global | complete; `jobs/15` **ALL CLEAR**, 12/12 stores |
| `jobs/20` full-year calibration | 117 min → `thresholds_2026-09-07.json` + 41 GB tail archive |
| `jobs/21` per-diagnostic trend | Annual 72 min, DJF 20 min |
| share | 928 GB of 2.5 TB, 38 % |

`1101414` ran 54 % slower than `1098646` on identical data and code. Cause:
three `jobs/21` tasks were reading **the same** `derived/north_atlantic` stores
concurrently. Once they cleared, the rate returned to ~6 min/year. Worth
remembering when sizing a walltime — **what else is running against the same
data matters as much as the job itself**, and at the original 4 h limit this
run would have died at 33 of 42 years.
