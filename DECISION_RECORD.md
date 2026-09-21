# DECISION RECORD — the CAT indicator dataset

**Decided 2026-09-16.** One page. What the dataset is, why these choices, what
the evidence is, and what must be disclosed. Detail lives in
`FORMULAS_AND_DECISIONS.md`; the run history is `STATUS.md` §18.

---

## 1. THE DATASET

```
/scistor/SBE-EDS-ClimateKoopman/yen230/derived/north_atlantic_audit/
    504 zarr stores, 1979-01 .. 2020-12, 3-hourly, 200 hPa, 21 diagnostics, float32
    store attributes: audit_fixes = meridional_metric+four_variable_provenance
                      f2d_variant = A
                      deformation_convention = DEF (un-squared, Sharman A17)

/scistor/SBE-EDS-ClimateKoopman/yen230/calibration/thresholds_2026-09-11.json
    cos(phi)-weighted global percentiles of the full contiguous year 2000,
    computed under the same conventions
```

**This is the dataset. Use it.** `derived/north_atlantic/` +
`thresholds_2026-09-07.json` is the same 42 years under the previously used
conventions; it is retained as the appendix sensitivity and is not otherwise
needed.

**Pairing rule, enforced in code.** A series must be scored against thresholds
built from a reference year computed under the same conventions.
`full_trend_check.py` and `per_diagnostic_trend.py` refuse a mismatch and exit
2. This is not pedantry: on 2026-09-15 a mismatched run produced a
plausible-looking table in which `magnitude_pv` silently had zero exceedance.

---

## 2. THE CHOICES

| | choice | in line with |
|---|---|---|
| Fields | ζ, δ, PV **computed** from u, v, T, z | Prosser (2023) p. 2 — he downloads those four and nothing else |
| Horizontal metric | MetPy spherical operator; meridional spacing as plate-carrée map distance | MetPy convention, with a documented departure on the spacing — 2e-7 vs 4.2e-3 against exact ellipsoidal geometry |
| Frontogenesis | signed material derivative (variant A) | Williams & Storer (2022) Eq. (3) p. 1427; Williams (2017) Table 2 units |
| Deformation | DEF, not DEF² | Sharman (2006) A17; Williams (2017) Table 2 units |
| Vertical stencil | 175/200/225 hPa at 200 hPa | **a constraint** — the finest available in the CDS pressure-level product |
| Thresholds | 97.0/99.1/99.6/99.8/99.9, cos φ-weighted, full reference year | Williams (2017) Table 1 |
| Comparison box | 36–60°N, 55–10°W | Prosser (2023) Fig. 2 caption |
| Signs | all 21 one-tailed upper | Williams (2017) Table 2 — all 21 ladders monotonically increasing |
| Endlich ∂ψ/∂z | angle secant (rojak) | the literal reading of Sharman A25 + "centred second-order finite differences"; the closed form is an appendix sensitivity |
| CP length scale λ | ICAO Δz — constant in space and time | Sharman A4 "the local value of vertical grid increment Δz"; the alternative reading is documented |

**On the field choice, stated honestly because a referee will ask.** ERA5's
archived ζ and δ are the IFS's own prognostic spectral variables and u, v are
derived *from* them (IFS CY49R1 Part III §2.2.7), so computing them back from
u, v is a lossy round trip rather than a refinement. The reasons to do it
anyway are that (i) it is the method of the paper being replicated, (ii) it
makes the 21-member ensemble internally consistent at one effective
resolution, which matters because the percentile calibration is
resolution-specific (Williams 2017), and (iii) **it is the only version with
positive external evidence** — see §3. There is no observational ground truth
for 200 hPa vorticity, so "which ζ is more correct" has no empirical answer;
Prosser is the only anchor that exists.

---

## 3. THE EVIDENCE

**Against every published quantity that can be compared: 28 inside our 95 %
interval, 0 outside** (`cat_outputs/SCORECARD_audit.md`).

| | audit series | Prosser |
|---|---|---|
| Table 1 cells inside 95 % CI | **25/25** | — |
| Table 1 cells significant at 5 % | **25/25** | — |
| diagnostics with significant upward trends | **17/21** | **17/21** |
| largest relative change | **+74 %** | +75.6 % |
| significant downward trends | **0** | 0 |
| annual trend ratio, LOG → SOG | 0.97 / 0.98 / 1.01 / 1.02 / 1.04 | — |
| annual 1979 level ratio, LOG → SOG | 0.904 / 0.888 / 0.880 / 0.865 / 0.842 | — |

**The trends are invariant to the convention choices.** Baseline → audit, the
25 trend ratios move by at most 0.05 and no verdict changes — across three
choices that alter 16 of the 21 derived fields and shift calibration
thresholds by up to 8 %.

**The levels improved at every severity.** Deficit 11.9–20.3 % → 9.6–15.8 %,
roughly a fifth to a quarter of the gap closed, monotonically. This was not
guaranteed: it was on record beforehand that A18 might make `magnitude_pv`
verify worse.

**The 17/21 match is the single strongest result.** On the previous
conventions `magnitude_pv` came in at t = 1.97 against a 2.021 critical value —
the only reason the count read 16/21. Under the four-field method it is +42 %,
t = 4.65, significant, and the count matches Prosser exactly.

**Verification of the code itself:** 90 tests pass; 16 diagnostics checked
against closed-form manufactured solutions; 5 cross-checked against an
independent implementation at ρ = 1.0000; the horizontal gradient operator
checked against exact WGS84 geometry.

---

## 4. PER-INDICATOR STATUS

**17 of 21 carry a significant positive trend and need no caveat.**

**Four do not reach significance at 5 %** — report them as such, not as
absent trends:

| | change | t | note |
|---|---|---|---|
| `horizontal_divergence` | +5 % | 1.01 | **Prosser gets the same flat line** (+6.7 %, p = 0.2). Not a defect — a property of the field. See below. |
| `colson_panofsky` | +62 % | 1.85 | sparse tail (0.039 % exceedance), wide interval |
| `f2d` | +25 % | 0.91 | sparse, and systematically damped ~36 % by 3-hourly sampling |
| `negative_richardson` | +41 % | 1.00 | sparsest of the 21 (0.019 %) |

**`horizontal_divergence` — the pre-registered question is CLOSED.** The
hypothesis was that its flatness came from ERA5's *archived* divergence being
weakly observation-constrained at 200 hPa, i.e. §5 risk 3, observing-system
changes across the record. Under computed divergence it is **still flat**
(+5 %, t = 1.01, z = −3.57 against the ensemble). The hypothesis loses its best
remaining candidate. Report it as a diagnostic without a trend, in both
studies.

**`ubf` is NOT anomalous.** §12.6 called it last on trend at +11 %; on the full
record it is rank 2 of 21 at +20 %, t = 2.33, significant. Third independent
confirmation. §12.6 was a nine-season endpoint artifact.

**`brown2` is rank-only.** Sharman A14 yields s⁻³ while the published tables
are m² s⁻³ — the dimensional gap is in the published equation, not the
pipeline, and the implied L² is constant so it cannot affect any exceedance
decision. Never compare its magnitude to a published J kg⁻¹ s⁻¹ figure.

**`magnitude_pv` carries two constant factors** (×100 from the hPa pressure
coordinate, and PVU = 10⁻⁶ SI). Both cancel in the calibration; neither is
applied. Do not quote its magnitude against a published PVU number.

---

## 5. WHAT MUST BE DISCLOSED

1. **The vertical stencil.** 175/200/225 hPa against Prosser's model levels
   73–75 (≈188/197/206 hPa, MARS-only). It biases the exceedance **level** by
   10–16 %, deepening with severity in winter and spring and **flat in summer**
   (JJA LOG−SOG spread −0.018 against DJF +0.097), and changes the fitted
   relative trend by at most one percentage point. Four independent lines of
   evidence, `FORMULAS_AND_DECISIONS.md` D6. Model levels via MARS are a
   costed, declined option.
2. **`f2d` is damped ~36 %** for sub-diurnal signals: a centred difference of a
   diurnal signal at 3-hourly sampling returns sin(Wh)/(Wh) = 0.637 of the true
   derivative. A property of the sampling; no code change repairs it.
3. **Three implementation choices are not specified anywhere in the CAT
   literature** — discretisation, field provenance, vertical stencil. That
   absence is itself reportable, and this dataset measures what each is worth.
4. **ERA5, not ERA5.1** for 2000–2006 (Prosser used ERA5.1; it is MARS-only).

---

## 6. NOT DONE, AND WHY IT DOES NOT BLOCK

- **Per-diagnostic 1979 level ratios on the audit series.** §17.5's investigate
  list was built on the baseline. Needs one join of
  `cat_outputs/per_diagnostic_annual_moderate_audit.csv` with
  `ada/prosser_published.py`. Minutes, no cluster. Affects the caveat list,
  not the dataset.
- **Provenance sweep of all 504 audit stores.** `jobs/15` hardcodes the
  baseline paths. Count, size and one store's attributes are verified; all 504
  are not.
- **Baseline `per_diagnostic_annual_moderate.csv`** was overwritten by the
  2026-09-15 mismatched run. Regenerate with
  `THRESHOLDS=.../thresholds_2026-09-07.json`. Appendix only.
- **Global and USA-box figures** (Prosser Figs 1, 2, 3b, S1, S2, S5) — a
  separate download, costed in `PLAN_prosser_replication.md`. Not needed for a
  North Atlantic result.
- **Sharman & Pearson (2017) Part I** and **Brown (1973)** — provenance, not
  correctness. Every formula is already traced to an equation this project has.

---

## 7. BEFORE THE ECONOMETRICS

**One thing genuinely gates it.** The peaks-over-threshold archive for the
North Atlantic does not exist — the tail archives are the global reference
year only. EVT or GPD work on 1979–2020 needs per-gridpoint excesses across
42 years (~23 GB, one pass over the existing series). Build it once the unit
of observation is settled, because the aggregation shape depends on that.

**And one risk that comes due at the same moment.** The ADA share has no
backup and no published retention policy (§16.3). While the derived data is
reproducible in principle that costs compute, not information — but the moment
a published result is fitted on the tail archive, "we could rebuild it in a
week" stops being an acceptable answer. Worth an email to ITvO before fitting
anything intended for publication.
