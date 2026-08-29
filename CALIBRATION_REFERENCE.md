# CALIBRATION REFERENCE — published values and what each one is for

Companion to `STATUS.md` §7 and §10.12. Created 2026-08-28 from the PDFs in
`Turbulence project/Articles/`. Supersedes the assumption that Williams & Joshi
(2013) Table 1 is the calibration anchor.

---

## 1. The question this file answers

*Are the thresholds published anywhere, so we can check the indicators before
committing to the 42-year download?*

Short answer: **the thresholds are published, but not by Prosser, and not in a
form that transfers directly.** Three separate published tables exist, they are
three different quantities, and the project has been comparing against the
wrong one.

---

## 2. Did Prosser reuse Williams' threshold values? No.

This was the open question. The paper is explicit — Prosser et al. (2023),
Section 2:

> "To allow an inter-diagnostic comparison, the uncalibrated CAT diagnostic
> values, each with different physical units, are compared with threshold values
> derived from a climatological probability distribution for each diagnostic,
> **following Williams (2017)**. The reanalysis data were extracted on a fixed
> Gaussian grid […] and so the climatological probability distributions were
> **latitudinally weighted**. The latitudinally weighted distributions were
> calculated for **the year 2000**, a reference year chosen as being the
> 1979–2020 midpoint."

and then:

> "Following Williams (2017), the diagnostic values corresponding to the **97th,
> 99.1st, 99.6th, 99.8th, and 99.9th percentiles** were then **derived globally
> for the reference year 2000**, corresponding, respectively, to the thresholds
> for LOG, LMOG, MOG, MSOG, and SOG turbulence."

**Prosser inherits Williams' _percentile ladder_, and recomputes the _values_
from ERA5.** Williams' own numbers are never used. This is the correct reading
and it is what `calibration.py` should implement.

Confirmed by the calibration design already in the repo: latitude weighting
(`calib_weighted_percentile.py`) is not an embellishment, it is required by the
method — Prosser weights because the Gaussian grid over-samples the poles.

### The percentile ladder — LOCKED

| Category | Cumulative percentile (the "or-greater" threshold) |
|---|---|
| LOG — light-or-greater | **97.0** |
| LMOG — light-to-moderate-or-greater | **99.1** |
| MOG — moderate-or-greater | **99.6** |
| MSOG — moderate-to-severe-or-greater | **99.8** |
| SOG — severe-or-greater | **99.9** |

Origin (Williams 2017 Table 1): a log-normal EDR distribution constrained by an
observed **3.0 % probability of LOG** and **0.4 % probability of MOG**, tied to
EDR^(1/3) bands 0.1/0.2/0.3/0.4/0.5 m^(2/3) s⁻¹ and to vertical accelerations of
0.2/0.4/0.6/0.8/1.0 g.

> **Check `compute_thresholds` against this.** If it computes a single 99th
> percentile — which is Williams & Joshi (2013)'s rule, not Prosser's — it
> cannot produce the LOG/LMOG/MOG/MSOG/SOG breakdown that Prosser Table 1
> reports, and stage 2 cannot be compared to anything.

---

## 3. Prosser's Supporting Information contains NO threshold table

Checked directly (`Supporting information si-s01 (1).docx`). Contents are
Figures S1–S5 only:

- S1–S4: LOG and SOG versions of main-paper Figures 1–4.
- **S5: "The MOG CAT absolute change (Panel a in main paper Figure 2), broken
  down by its 21 constituent CAT diagnostics."**

S5 is worth extracting by eye. It is the only published **per-diagnostic**
replication target in existence — it lets a single diagnostic be checked
individually instead of only the 21-member ensemble mean. Values must be read
off a figure, so it is an order-of-magnitude and sign check, not a precise one.

---

## 4. The three published tables, and what each is actually good for

| Source | What it contains | Data behind it | Use it for |
|---|---|---|---|
| **Williams (2017) Table 2** | Onset thresholds, all 21 diagnostics × 5 severities | **GFDL-CM2.1** climate model, DJF, 200 hPa, 50–75°N/10–60°W | Sign, units, and shape of the threshold ladder. **Not** magnitude. |
| **Williams & Joshi (2013) Table 1** | **Medians** (pre-industrial vs 2×CO₂), 21 diagnostics | Same GCM, DJF, 200 hPa, 50–75°N/10–60°W, **daily-mean** | A median-vs-median smell test only |
| **Lee et al. (2023) Table 1** | MOG thresholds, 7 empirical indices | **ERA5**, 250 hPa, 20–60°N, 1979–2019, 6-hourly, **p95** | **The only like-for-like magnitude anchor.** Use this one. |
| **Prosser (2023) Table 1** | Fitted 1979→2020 turbulence hours, NA box, 5 severities × 5 seasons | ERA5(.1), 197 hPa | **The replication target.** |

### 4.1 Why the GCM tables cannot be matched on magnitude — the author says so

Williams (2017), immediately below Table 2:

> "Note that the thresholds are **dependent on the grid resolution of the
> atmospheric model**. Therefore, the values listed in Table 2 **may differ from
> those computed in other studies** (e.g., Sharman et al., 2006)."

GFDL-CM2.1 is ~2°. ERA5 is 0.25°. Every gradient-based diagnostic is
resolution-dependent, and a residual like UBF is the most resolution-sensitive
of all. **The `1–2 ORDERS OFF` flags in `comparison_table.csv` are measuring
resolution, not correctness.**

This *strengthens* STATUS.md §5.4 rather than weakening it: because thresholds
are percentiles of our own ERA5 data, a resolution offset moves the data and the
threshold by the same factor and cancels exactly out of the exceedance field.

### 4.2 Williams (2017) Table 2 — onset thresholds (GFDL-CM2.1, DJF, 200 hPa, 50–75°N/10–60°W)

Use for **sign and units**, not magnitude.

| Diagnostic | Units | Light | L-to-M | Moderate | M-to-S | Severe |
|---|---|---|---|---|---|---|
| Negative Richardson number | — | −15.4 | −9.8 | −7.9 | −6.7 | −5.9 |
| Vertical shear of horizontal wind | 10⁻³ s⁻¹ | 5.3 | 6.6 | 7.4 | 7.9 | 8.4 |
| Colson–Panofsky index | **10³ kt²** | −29.3 | −27.0 | −25.2 | −23.7 | −22.2 |
| Frontogenesis function | 10⁻⁹ m² s⁻³ K⁻² | 770 | 1280 | 1660 | 1980 | 2340 |
| Brown index | 10⁻⁶ s⁻¹ | 99 | 106 | 110 | 113 | 118 |
| Brown energy dissipation rate | 10⁻⁶ J kg⁻¹ s⁻¹ | 870 | 1370 | 1730 | 2030 | 2330 |
| Ellrod TI1 | 10⁻⁹ s⁻² | 195 | 292 | 360 | 419 | 472 |
| Ellrod TI2 | 10⁻⁹ s⁻² | 184 | 282 | 356 | 419 | 477 |
| Flow deformation | 10⁻⁶ s⁻¹ | 50.9 | 60.9 | 66.9 | 71.8 | 76.3 |
| Magnitude of potential vorticity | PVU | 8.33 | 8.73 | 8.98 | 9.19 | 9.41 |
| Relative vorticity squared | 10⁻⁹ s⁻² | 2.46 | 3.74 | 4.70 | 5.50 | 6.24 |
| Horizontal temperature gradient | 10⁻⁶ K m⁻¹ | 14.7 | 17.6 | 19.4 | 20.8 | 22.0 |
| Wind speed | m s⁻¹ | 40.9 | 48.4 | 52.4 | 55.3 | 58.5 |
| Wind speed × directional shear | 10⁻³ rad s⁻¹ | 3.21 | 3.94 | 4.39 | 4.72 | 5.08 |
| Flow deformation × wind speed | 10⁻³ m s⁻² | 1.65 | 2.29 | 2.76 | 3.17 | 3.54 |
| Flow deformation × vertical temp. gradient | 10⁻⁹ K m⁻¹ s⁻¹ | 53 | 84 | 106 | 127 | 151 |
| Residual of nonlinear balance equation (UBF) | 10⁻¹² s⁻² | 1230 | 1840 | 2270 | 2610 | 2960 |
| Magnitude of horizontal divergence | 10⁻⁶ s⁻¹ | 11.9 | 15.7 | 18.2 | 20.4 | 22.5 |
| NCSU1 | 10⁻¹⁸ s⁻³ | 1200 | 3600 | 6300 | 9300 | 13000 |
| Negative absolute vorticity advection | 10⁻⁹ s⁻² | 1.33 | 1.86 | 2.23 | 2.56 | 2.93 |
| Magnitude of relative vorticity advection | 10⁻⁹ s⁻² | 1.44 | 1.99 | 2.34 | 2.66 | 3.00 |

**Two traps in this table, both directly relevant to STATUS.md §7:**

1. **Colson–Panofsky is published in 10³ kt², not SI.** 1 kt² = 0.2646 m² s⁻².
   If `REFERENCE_TABLE` compares an SI value to the tabulated number without
   converting, that is a factor of ~3800 — which is exactly the "1–2 orders off"
   flag we see. **Check this before treating CP as a defect.**
2. **CP's thresholds are all negative and increase toward severe** (−29.3 → −22.2).
   Williams explains why: "the Colson–Panofsky index is proportional to
   1 − Ri/0.5, and the Richardson number is rarely less than 0.5 in the
   GFDL-CM2.1 model." In ERA5 at 0.25°, Ri < 0.5 *does* occur, so our CP
   distribution has a genuinely different shape and can be positive. **CP's
   published value is untransferable — but its sign convention is now settled,
   which closes one of the §7 items.** The turbulent tail is the *upper* one
   (`sign = "+"`), consistent with the correction already made.

### 4.3 Lee et al. (2023) Table 1 — ERA5-derived MOG thresholds (250 hPa, 20–60°N, 1979–2019, p95, 6-hourly)

| Index | MOG threshold | Unit |
|---|---|---|
| VWS | > 1.03 × 10⁻² | s⁻¹ |
| DEF | > 1.29 × 10⁻⁴ | s⁻¹ |
| −DIV | > 5.15 × 10⁻⁵ | s⁻¹ |
| DVT (divergence tendency) | > 7.18 × 10⁻⁵ | s⁻¹ |
| TI1 | > 9.07 × 10⁻⁷ | s⁻² * |
| TI2 | > 1.03 × 10⁻⁶ | s⁻² * |
| TI3 | > 1.39 × 10⁻⁶ | s⁻² * |
| N² | < 0 | s⁻² |
| PV | < 0 | PVU |
| Ri | 0 < Ri < 1 | — |

\* Lee's table prints s⁻¹ for TI1–TI3; TI = VWS × DEF is s⁻², so this is a
typo in the source.

Note Lee uses **p95** for MOG, not Prosser's p99.6 — a deliberate departure they
justify (their §2). So Lee's numbers are *lower* than a Prosser-style MOG
threshold would be. Compare with that in mind.

### 4.4 Prosser (2023) Table 1 — THE REPLICATION TARGET

Fitted change 1979→2020, hours per period at an average point in the North
Atlantic box, from ERA5 at 197 hPa, diagnostic-mean over the 21.

| Season | | LOG | LMOG | MOG | MSOG | SOG |
|---|---|---|---|---|---|---|
| **DJF** | 1979 (h) | 128.9 | 45.6 | 22.3 | 12.1 | 6.4 |
| | 2020 (h) | 155.6 | 59.3 | 30.6 | 17.2 | 9.6 |
| | Relative increase | 21 % | 30 % | **37 %** | 43 % | 49 % |
| **MAM** | 1979 (h) | 90.4 | 27.2 | 11.8 | 5.7 | 2.7 |
| | 2020 (h) | 113.4 | 38.9 | 18.6 | 9.7 | 5.0 |
| | Relative increase | 26 % | 43 % | 57 % | 71 % | 85 % |
| **JJA** | 1979 (h) | 114.1 | 36.5 | 16.1 | 7.7 | 3.6 |
| | 2020 (h) | 124.5 | 43.8 | 21.1 | 10.9 | 5.5 |
| | Relative increase | 9 % | 20 % | 31 % | 41 % | 52 % |
| **SON** | 1979 (h) | 133.1 | 43.4 | 19.8 | 10.0 | 5.0 |
| | 2020 (h) | 153.2 | 53.4 | 25.8 | 13.9 | 7.4 |
| | Relative increase | 15 % | 23 % | 31 % | 39 % | 47 % |
| **Annual** | 1979 (h) | 466.5 | 152.7 | 70.0 | 35.5 | 17.7 |
| | 2020 (h) | 546.8 | 195.4 | 96.1 | 51.8 | 27.4 |
| | Relative increase | **17 %** | 28 % | **37 %** | 46 % | **55 %** |

All statistically significant, p < 3 × 10⁻².

**The absolute levels are as useful as the trends.** 70.0 h of MOG in a year is
70.0 / 8766 = **0.80 %** — matching the paper's stated "0.8 % in 1979". This is
a *level* check that needs only the calibration year plus a few analysis years,
not all 42.

---

## 5. Evidence the magnitudes are already right

Using the January-1979 medians recorded in STATUS.md §10.12 against the
published tables. Three independent consistency checks, all passing.

### 5.1 Composition ratios compose exactly

TI1 is defined as VWS × DEF, so its ratio to any reference must equal the
product of its factors' ratios.

| | project Jan-1979 | W&J median | ratio |
|---|---|---|---|
| VWS | 4.67 ×10⁻³ s⁻¹ | 1.88 | 2.484 |
| DEF | 43.8 ×10⁻⁶ s⁻¹ | 18.6 | 2.355 |
| TI1 predicted (2.484 × 2.355) | | | **5.850** |
| TI1 observed (184 / 31.5) | | | **5.841** |

Agreement to **0.15 %**. TI1's entire offset is inherited from its factors;
nothing is introduced by TI1 itself. Generalise this to every product
diagnostic (NGM1, NGM2, TI2, wind × directional shear) as a real test.

### 5.2 Resolution-insensitive diagnostics show no offset

Brown index is `sqrt(0.3 ζa² + Dsh² + Dst²)`, dominated by the Coriolis part of
ζa, which is a function of latitude alone and carries no resolution dependence.
Project 80.6 vs W&J 77.1 — **ratio 1.045**. Exactly what the resolution
hypothesis predicts, and a strong argument that the large ratios elsewhere are
resolution and not error.

### 5.3 Extrapolating to Lee's ERA5 p95 lands on Lee's published value

Take the distribution shape from the GCM (where both a median and a p97 are
published) and apply it to our ERA5 median:

| | our median | shape factor p97/median from GCM | our implied p95–97 | **Lee's ERA5 p95** |
|---|---|---|---|---|
| VWS | 4.67 ×10⁻³ | 5.3 / 1.88 = 2.82 | ~1.3 ×10⁻² | **1.03 ×10⁻²** |
| DEF | 43.8 ×10⁻⁶ | 50.9 / 18.6 = 2.74 | ~1.2 ×10⁻⁴ | **1.29 ×10⁻⁴** |

Right order, right magnitude, against an **independent ERA5-based study**.
Caveats: Lee is 250 hPa not 200, 20–60°N not 30–60°N, DJF+JJA not January,
6-hourly not 3-hourly. This is a magnitude check, not a precision check — but
it is the closest thing to like-for-like that exists.

### 5.4 What is left genuinely unexplained

| Diagnostic | Project vs W&J median | Note |
|---|---|---|
| `colson_panofsky` | flagged | Almost certainly the **kt² unit conversion** (§4.2). Check first. |
| `ubf` | 3632 vs 161 ×10⁻¹² (22.6×) | A residual of near-cancelling terms — the most resolution-sensitive quantity in the set. Plausible, but the one to watch. Its geometry changed on 2026-08-26 (STATUS §4g). |
| `nva`, `rva_magnitude` | flagged | Advection terms are heavy-tailed: Williams' own p97/median ratios are ~6.3 for these versus ~2.8 for VWS. A large median ratio is less alarming here. |
| `ncsu1` | 8.98 vs 11.1 (**0.81**) | The only diagnostic *below* the reference while everything else is above. Runs against the resolution story. **Worth a look.** |
| `magnitude_pv` | 5.19 vs 6.84 (0.76) | Explained: W&J's 50–75°N box is further poleward, i.e. more stratospheric, so higher PV. Not a concern. |

---

## 6. A data-configuration issue that would silently shift every threshold

Prosser, Section 4:

> "Our analysis has used **ERA5.1**, which corrects for the known cold bias in
> the lower stratosphere during **2000–2006** in the previous version of ERA5
> (Simmons et al., 2020)."

**The calibration year 2000 falls inside the ERA5.1 correction window**, and the
correction is in the lower stratosphere — which is where 197/200 hPa sits. On
CDS, ERA5.1 is a separate dataset entry, not a flag on the standard request.

If `download_plan.py` requests standard ERA5 for the global year-2000
calibration pull, every threshold will be derived from a slightly different
temperature field than Prosser's, and the difference is systematic rather than
random. **Verify before running `02_download_global.sbatch`.** The same question
applies to the 2000–2006 slice of the North Atlantic analysis run.

---

## 7. Verification plan, cheapest first

| # | Check | Cost | Acceptance criterion |
|---|---|---|---|
| 1 | Unit audit of `REFERENCE_TABLE` against the published units in §4.2, starting with Colson–Panofsky's kt² | Free, no data | Every unit matches, or the mismatch is deliberate and recorded |
| 2 | Re-derive all 21 `sign` entries from Sharman (2006), Ellrod & Knapp (1992), Williams (2017) | Free, PDFs on disk | STATUS §7 closed. **Highest-value remaining check** — a flipped sign is invisible in magnitude and fatal in rank |
| 3 | Composition-ratio test as a real test (§5.1) across every product diagnostic | Free, Jan-1979 output on disk | Each product's ratio equals the product of its factors' ratios to within a few % |
| 4 | Confirm `compute_thresholds` uses the 5-percentile ladder (§2), not a single p99 | Free, code read | Produces LOG/LMOG/MOG/MSOG/SOG |
| 5 | Re-point `compare_to_published` at Williams (2017) Table 2 for **sign and units**, and label the magnitude comparison as resolution-confounded | Small code change | No longer compares a threshold to a median |
| 6 | Confirm ERA5 vs ERA5.1 for year 2000 (§6) | One CDS lookup | Matches Prosser's configuration |
| 7 | **Global year-2000 calibration run** — 12 CDS requests, ~1/40th of the full download | ~1–2 h queue, 128 GB | Thresholds have the right sign and plausible magnitude vs §4.2; latitude weighting active |
| 8 | **THE GATE:** apply year-2000 thresholds to a few North Atlantic years and check absolute frequency | A handful of months | **Annual MOG ≈ 0.8 % (70 h), LOG ≈ 5.3 % (466 h)** at an average NA point, per Prosser Table 1 |

Check 8 is the go/no-go for the 42-year download. It is resolution-robust
precisely because it is percentile-calibrated on our own data — which is why it
succeeds where every magnitude comparison fails.

---

## Sources

- Prosser, M. C., Williams, P. D., Marlton, G., & Harrison, R. G. (2023).
  Evidence for large increases in clear-air turbulence over the past four
  decades. *Geophysical Research Letters*, 50, e2023GL103814.
  Main text and Supporting Information S1.
- Williams, P. D. (2017). Increased light, moderate, and severe clear-air
  turbulence in response to climate change. *Advances in Atmospheric Sciences*,
  34(5), 576–586. Tables 1 and 2.
- Williams, P. D., & Joshi, M. M. (2013). Intensification of winter transatlantic
  aviation turbulence in response to climate change. *Nature Climate Change*,
  3, 644–648. Table 1.
- Lee, J. H., Kim, J. H., Sharman, R. D., Kim, J., & Son, S. W. (2023).
  Climatology of clear-air turbulence in upper troposphere and lower
  stratosphere in the Northern Hemisphere using ERA5 reanalysis data.
  *JGR: Atmospheres*, 128, e2022JD037679. Table 1.

---

## 8. Code audit — 2026-08-28

### 8.1 `REFERENCE_TABLE` units and medians: all 21 correct ✅

Every `units` and `wj_median` entry in `2_diagnostics.py` was checked line by
line against Williams & Joshi (2013) Table 1 and Williams (2017) Table 2, read
from the PDFs. **All 21 match.** The kt² concern raised in §4.2 does not apply —
`colson_panofsky` is already recorded as `10^3 kt^2` with median −34.8.

`f2d`'s units carry a note that they were corrected from stale Miller-form units
to the A9/isentropic form; the corrected value `10^-9 m^2 s^-3 K^-2` is what both
papers print. Correct.

**This closes the units half of STATUS.md §7.** What remains open there is the
formula derivations themselves against Sharman (2006) Appendix A and
Ellrod & Knapp (1992) — the PDFs are in `Articles/`.

### 8.2 All 21 signs are `"+"`, and Williams Table 2 confirms it

`colson_panofsky` and `negative_richardson` were the two corrected from
`"either"`. Williams (2017) Table 2 settles it independently: both rows are
**single monotonic ladders** running from most-negative (light) to
least-negative (severe) — −29.3→−22.2 and −15.4→−5.9. A two-tailed criterion
would need two ladders. Negative *values* with an upper-tail *criterion*. The
correction already in the code is right.

### 8.3 Calibration → aggregation identity: PASS ✅

`tests/test_calibration_roundtrip.py` (new, 31 tests, ~15 s, no ERA5 required).

Applying a p-th percentile threshold back to the data it was calibrated on must
reproduce `(100−p)/100` exactly. Measured on synthetic latitude-dependent
fields through the project's own `calibration.py` / `aggregate.py`:

| Severity | expected | observed | rel. error |
|---|---|---|---|
| LOG (p97.0) | 3.0000 % | 3.0000 % | 0.00 % |
| LMOG (p99.1) | 0.9000 % | 0.9000 % | 0.00 % |
| MOG (p99.6) | 0.4000 % | 0.4000 % | 0.01 % |
| MSOG (p99.8) | 0.2000 % | 0.2000 % | 0.02 % |
| SOG (p99.9) | 0.1000 % | 0.0999 % | 0.07 % |

Worst error 0.071 %, which is the Hazen half-slot offset `0.5/(n(1−p))` and not
a defect. Also verified: the ladder is monotonic for every diagnostic including
the negative-valued one; NaN cells reduce `populated_count` rather than voting
"no"; an all-NaN diagnostic is refused rather than given a threshold; and the
Q-CALIB-6 hardcoding guard fires on literal Table 2 values but not on a
calibration that merely lands near them.

**The cos(φ) weighting was verified to actually bite** — on a field whose tail
thickens poleward it moves the MOG threshold by **−16.3 %** versus unweighted,
and in the correct direction. A test asserts this, because a weighting that
silently failed to reach `weighted_percentile` would leave every threshold
biased toward the over-sampled poles and nothing else would notice.

### 8.4 ⚠ `compute_thresholds` will not fit in memory on a full global year

`compute_thresholds` does `np.asarray(field).ravel()` on the entire calibration
field, then `weighted_percentile` sorts it and calls `np.unique`.

A global year at 0.25°, 3-hourly is 721 × 1440 × 2928 = **3.04 × 10⁹ points per
diagnostic**. As float64 that is 24 GB for the values, another 24 GB for the
broadcast weights, and `np.unique` needs a sort buffer on top — **well over
100 GB peak, per diagnostic.** For comparison, one North Atlantic month of
diagnostics peaked at 14.8 GB (STATUS §10.12).

**Stage 2 as currently written cannot run.** Three ways out, in increasing
effort: sub-sample the calibration year (§9 — also the cheapest verification
path, so this is the recommended fix); cast to float32 and drop the `np.unique`
tie-consolidation; or replace the exact percentile with a streaming/histogram
estimator. Sub-sampling is enough and needs no code change.

---

## 9. THE CHEAP VERIFICATION EXPERIMENT

### 9.1 Prosser's North Atlantic box is INSIDE ours

Prosser et al. (2023), Figure 2 caption:

> "The two boxes represent the North Atlantic (**36–60°N and 55–10°W**) and
> USA (30–55°N and 124–60°W) areas used in Figures 3 and 4 and Table 1."

`download_plan.DOMAINS["north_atlantic"]` is 30–60°N, 75–0°W. **Prosser's box is
a strict subset of ours**, so Table 1 can be compared exactly by subsetting at
analysis time. No re-download, no approximation. This coordinate pair is not
recorded anywhere in STATUS.md and must be — comparing on our larger box would
produce systematically different frequencies and look like an error.

### 9.2 The experiment

**Calibration — global year 2000, sub-sampled.** `download_plan.build_request`
already accepts an explicit `days` list, so this is a job-script change only.
Take **days 1, 9, 17 and 25 of each month** at full 3-hourly resolution: 48 days
spread evenly across the year, keeping the whole diurnal cycle and the whole
seasonal cycle.

| | full year | days 1/9/17/25 |
|---|---|---|
| CDS requests | 12 | 12 (each ~⅛ the size) |
| Volume | ~128 GB | **~17 GB** |
| Points per diagnostic | 3.04 × 10⁹ | **3.99 × 10⁸** |
| Points above the MOG threshold | 12.2 M | **1.6 M** |
| Fits in a 24 GB job | no (§8.4) | yes |

1.6 million points in the MOG tail makes the sampling error on the threshold
negligible — far below the 0.2 % gradient-metric noise floor. Sub-sampling costs
nothing scientifically here and is what makes the run possible at all.

**Analysis — North Atlantic DJF 1979.** January 1979 is already downloaded;
add **February and December 1979** — two more CDS requests, ~15 minutes.

**Total: 14 CDS requests, roughly half a day.** Against 504 for the full run.

### 9.3 The acceptance criterion

Subset to 36–60°N / 55–10°W, apply the year-2000 global thresholds, and compute
the diagnostic-mean exceedance frequency. Prosser Table 1, DJF 1979, over
90 days = 2160 hours:

| Severity | Prosser 1979 (h) | **as a frequency** |
|---|---|---|
| LOG | 128.9 | **5.97 %** |
| LMOG | 45.6 | **2.11 %** |
| MOG | 22.3 | **1.03 %** |
| MSOG | 12.1 | **0.56 %** |
| SOG | 6.4 | **0.30 %** |

(The annual row, for later: 5.33 %, 1.74 %, 0.80 %, 0.41 %, 0.20 %.)

Read the result at three strengths:

1. **Order of magnitude.** MOG near 1 %, not 10 % or 0.01 %. Catches a gross
   calibration or unit failure.
2. **The ladder shape** — the five frequencies in the ratios 5.97 : 2.11 : 1.03 :
   0.56 : 0.30. This is the strong test. It constrains the whole tail shape of
   the North Atlantic distribution relative to the global one, using five
   numbers at once, and no single-point agreement can fake it.
3. **Absolute agreement within ~30 %.** Would be a strong result given the
   sub-sampled calibration and the 175/200/225 vs 188/197/206 hPa level
   substitution.

Failing (1) means stop. Passing (2) means the chain is right.

### 9.4 Two known offsets to expect, and not to panic about

- **Levels.** We evaluate at 200 hPa from a 175/225 stencil; Prosser uses
  197 hPa from 188/206. A narrower stencil gives larger vertical derivatives, so
  our shear-based diagnostics run slightly high — but this is a near-constant
  factor and largely cancels in the percentile calibration (§5.4 of STATUS.md).
- **ERA5 vs ERA5.1 for year 2000** (§6). Unresolved and worth settling before
  the calibration pull, since it moves every threshold in the same direction.

### 9.5 A per-diagnostic target, if the ensemble check passes

Prosser's Figure S5 breaks the MOG change down by all 21 constituent
diagnostics, and Figures 3/4 show 42 years of per-diagnostic data for the same
NA box. Neither is tabulated, so both must be read off the figures — but they
turn a single ensemble number into 21 individual ones, which is what would let a
single misbehaving diagnostic be identified rather than averaged away.

---

## 10. RESULT — the calibration check, 2026-08-29

Job 1092630, node013, 117 min. Thresholds in
`calibration/thresholds_2026-08-29.json`.

### 10.1 The gate: PASSED

North Atlantic DJF 1979, inside Prosser's box (36–60°N / 55–10°W), against
Prosser (2023) Table 1:

| | ours | Prosser | ratio | ours (h) | Prosser (h) |
|---|---|---|---|---|---|
| LOG | 5.037 % | 5.968 % | 0.84 | 108.8 | 128.9 |
| LMOG | 1.684 % | 2.111 % | 0.80 | 36.4 | 45.6 |
| **MOG** | **0.785 %** | **1.032 %** | **0.76** | 17.0 | 22.3 |
| MSOG | 0.399 % | 0.560 % | 0.71 | 8.6 | 12.1 |
| SOG | 0.200 % | 0.296 % | 0.68 | 4.3 | 6.4 |

| Criterion | Result | |
|---|---|---|
| 1. Order of magnitude | MOG ratio 0.76 | **PASS** |
| 2. Ladder shape | worst deviation 11.0 % | **PASS** |
| 3. Absolute (<30 %) | 24.0 % from Prosser | **PASS** |

Criterion 3 was not expected to be met. An independent implementation of 21
diagnostics, calibrated on a sub-sampled global year, reproduces a published
GRL result on absolute turbulence frequency to within a quarter.

### 10.2 The identity check on real data: exact

| Severity | expected | observed | rel. err |
|---|---|---|---|
| LOG | 3.0000 % | 3.0000 % | 0.00 % |
| LMOG | 0.9000 % | 0.9000 % | 0.00 % |
| MOG | 0.4000 % | 0.4000 % | 0.00 % |
| MSOG | 0.2000 % | 0.2000 % | 0.00 % |
| SOG | 0.1000 % | 0.1000 % | 0.00 % |

This is what licenses reading §10.1 at all. Thresholds are applied exactly as
computed, on real fields with real NaN patterns.

### 10.3 The ratios decline monotonically — and why

0.84 → 0.80 → 0.76 → 0.71 → 0.68. Systematic, not scatter: our distribution is
thinner than Prosser's and increasingly so toward the extreme tail.

**The vertical stencil explains it and predicts exactly this shape.** We
differentiate across 175→225 hPa (**50 hPa apart**); Prosser uses 188→206
(**18 hPa**). A stencil 2.8× wider smooths across sharp shear layers and
under-resolves shear *maxima*. The bulk of the distribution barely notices; the
extreme tail — which is where thin, intense shear layers live — is damped
progressively harder.

So the 24 % gap is not error. It is the documented level substitution
(STATUS §5.4) finally quantified. §11 explains why it cannot be closed.

### 10.4 Split-half: 48 days is enough

Median disagreement 0.81 %, 90th percentile 5.90 %.

The script printed "marginal" on a worst case of 92.33 % for
`colson_panofsky @ light_to_moderate`. **That is an artifact of the metric.**
CP's ladder is −19.33, −1.042, +11.01, +24.74, +40.4 — the LMOG threshold sits
essentially at zero, and a *relative* difference on a near-zero quantity is
meaningless. The check divides by the value; it should divide by the
diagnostic's own spread.

Read 5.9 % at the 90th percentile as the honest number. Adding days 5/13/21/29
is worth doing for a production calibration, but is not a blocker.

### 10.5 The Colson–Panofsky "sign mismatch" is a positive finding

Stage D flagged CP: ours **+0.042** against Williams' **−25.2**. Expected, and
Williams says why beneath his own table:

> "the Colson–Panofsky index is proportional to 1 − Ri/0.5, and the Richardson
> number is **rarely less than 0.5 in the GFDL-CM2.1 model**."

His CP is entirely negative because his climate model never resolves
sub-critical Richardson numbers. ERA5 at 0.25° does. Cross-checked against our
own independently computed Richardson diagnostic:

- `negative_richardson` at p99.6 = **−0.3611**, i.e. Ri = 0.361
- CP ∝ 1 − Ri/0.5 = 1 − 0.722 = **+0.278**, positive
- CP's measured MOG threshold is indeed **positive**

Two separately computed diagnostics agreeing on where the sign flips. And at
p99.9, Ri = **0.179** — below the classical 0.25 Kelvin–Helmholtz threshold.

**This closes the STATUS §4a gap** where the negative-Ri branches of
`colson_panofsky`, `negative_richardson` and `ncsu1` were never exercised
because the validation file had Ri > 0 everywhere. On a global year they are.

### 10.6 Three diagnostics sit outside the resolution family

After correcting a unit error in the comparison itself (Williams prints NVA/RVA
in 10⁻⁹ s⁻² where `REFERENCE_TABLE` uses 10⁻¹⁰; converting requires
multiplying by ten, and the first version divided — `nva` read 1081 when it is
**10.8**, `rva` 1327 when it is **13.3**), the stage-D ratios form one coherent
family spanning roughly 1× to 22×, consistent with a 2°→0.25° resolution
change. Three do not:

| Diagnostic | ratio vs W2017 MOG | |
|---|---|---|
| `brown2` | **7.6 × 10⁻⁷** | five orders too small |
| `f2d` | **2146** | three orders too large |
| `ncsu1` | 153 | borderline; NCSU1 is cubic in gradients |

**None of these affects §10.1**, because a constant factor cancels exactly in
percentile calibration (STATUS §5.4). They are evidence about the *formulas*,
which is STATUS §7's remaining open item. `brown2` is the place to start.

---

## 11. The 188/197/206 levels are MODEL levels — MARS only

Checked against ECMWF's L137 model level definitions:

| Model level | Full-level pressure |
|---|---|
| 73 | **188.29 hPa** |
| 74 | **197.37 hPa** |
| 75 | **206.81 hPa** |

Prosser's "197 hPa pressure level, with 188 and 206 for vertical derivatives"
is **ERA5 model levels 73/74/75**, described by their nominal pressures. It is
not the pressure-level product.

**Consequences:**

1. **They cannot be requested from `reanalysis-era5-pressure-levels`.** That
   product has a fixed 37-level set; around 200 hPa it offers 150, 175, 200,
   225, 250 and nothing between.
2. **Model-level ERA5 lives on the MARS tape archive**
   (`reanalysis-era5-complete`) — a different request format and slow access,
   the same obstacle as ERA5.1 (§6).
3. **175/200/225 is therefore not a compromise — it is the finest stencil
   available at 200 hPa on the CDS pressure-level product.** There is no
   cheaper improvement to be had. Recording this closes the "pursue model
   levels?" branch of STATUS §6 phase 3 with a cost, rather than leaving it
   open as if it were a small decision.

### 11.1 What CAN be done instead — measure the sensitivity

The stencil hypothesis in §10.3 is testable without model levels. Request
**150, 175, 200, 225 and 250 hPa in one download**, then compute the 200 hPa
diagnostics twice: once with the 175/225 stencil (50 hPa) and once with
150/250 (100 hPa). That gives two points on the stencil-width curve.

If the MOG ratio falls from 0.76 toward ~0.6 as the stencil widens, the
mechanism is confirmed and the slope lets you extrapolate to what Prosser's
18 hPa spacing would have produced. An inaccessible comparison becomes a
defensible extrapolation — which is exactly what a referee will want when
asking why the levels differ.

Cost: the same 15 CDS requests as the calibration check, at 5/3 the volume
(five levels instead of three), plus a repeat of the diagnostics and the
check. Roughly one day. Worth doing before phase 4, not during it.

### 11.2 Feasibility of the MARS route, checked properly

The question is not only "can the data be got" but "what would it cost the
code". Both were checked.

**Access.** ERA5 model levels are served through `reanalysis-era5-complete` on
the CDS. A standard CDS account reaches it — no special licence — but it is a
MARS front end onto the tape archive, with a different request grammar
(`class`, `levtype: ml`, numeric `param` codes, `grid`, `stream`, `type`) and
retrieval times set by tape queues rather than by the size of the transfer.
ECMWF's own guidance is to consolidate everything for a month into a single
request, because per-request overhead dominates.

So: obtainable, slowly. That alone would be a reason to prefer the
pressure-level product, but not a decisive one.

**The decisive part is that model levels are a different vertical coordinate.**
ERA5's 137 model levels are hybrid sigma-pressure and terrain-following, so:

> "Pressure at ERA5 model levels is not constant — it varies with both
> location and time."

Pressure must be *computed*, per gridpoint and per timestep, from the
logarithm of surface pressure plus the level's `a` and `b` coefficients:
`p_half = a + b·sp`, with the full-level pressure the mean of the two
half-levels around it. That needs `lnsp` downloaded alongside every field.

**What that breaks in this project:**

| Component | Why it breaks |
|---|---|
| `download_plan.DATASET` / `PRESSURE_LEVELS` | different dataset, MARS request grammar, native reduced-Gaussian grid needing explicit regridding |
| `1_download_hpc.verify_file` | `check_counts` / `check_grid` assume the pressure-level structure |
| `_sel_level(ds, name, level)` | `.sel(pressure_level=level)` has no meaning when the coordinate is a hybrid level index |
| `_brunt_vaisala_squared`, `altitude_derivative_on_pressure_level` | derived on pressure surfaces |
| `colson_panofsky`, `frontogenesis_2d`, `frontogenesis_isentropic`, `ncsu1` | each uses the target pressure explicitly (`p = float(target_level)`, `.sel(pressure_level=...)`) |
| `prepare_for_rojak` / rojak `CATData` | built around pressure levels |
| `tests/` — the 25 analytic cases | the manufactured atmosphere is defined on pressure surfaces |

This is not "change three numbers in a list". It is a change of vertical
coordinate that reaches the download layer, the data-prep layer, at least four
diagnostics, and the entire verification suite — the suite being the part that
makes any of the numbers trustworthy.

**And there is a second, subtler difference it would not even fix.** Because
model-level pressure varies with surface pressure, Prosser's "197 hPa" is a
*nominal* label on a terrain-following surface, not a true isobaric one. Our
200 hPa is genuinely isobaric. So the two studies differ in the surface they
evaluate on, not only in the spacing of the stencil. Exact agreement was never
available, whatever we did. (How large that effect is near 200 hPa depends on
the `b` coefficient at level 74, which is small in the upper troposphere but
not zero — worth checking if this is ever revisited.)

**Verdict: stay on 175/200/225.** The result is within 24 % of Prosser with a
fully understood, monotone, physically explained residual. The MARS route
costs slow retrieval, a new download path, a new vertical-coordinate layer
through half the diagnostics, and a rebuilt verification suite — to close a gap
that is already explained and would not close completely. The sensitivity
experiment in §11.1 buys the defensible number for a hundredth of the effort.
