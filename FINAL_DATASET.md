# FINAL DATASET — CAT indicators, North Atlantic, 1979–2020

**Written 2026-09-16, updated 2026-09-21 with the batch results.** This closes
the data-gathering phase. It does not restate `DECISION_RECORD.md`. It checks
that document's conclusion, corrects eleven statements in the record that are
wrong or overstated (§2), puts a number on every limitation (§3), gives each of
the 21 indicators a class (§5), and records what to keep (§8).

**Evidence status.** Job **1124312** (2026-09-21, 32 min, exit 0) ran the batch
of §4: the store sweep, the zonal tail profile, per-diagnostic fits at all five
severities, and the scorecard. Every number in §5 is now measured on
`derived/north_atlantic_audit` + `thresholds_2026-09-11.json`; nothing is
pending. The archived CSVs are `cat_outputs/final/`; the run is
`logs/final-1124312.out`; the figures come from
`ada/prosser_figures.py` and `ada/final_figures.py`, and
`notebooks/final_dataset_report.ipynb` renders them all in one pass.

---

## 0. The answer

1. **Use `derived/north_atlantic_audit/` with `thresholds_2026-09-11.json`.**
   The sweep certifies it: 504 of 504 stores under the audit fix set, 12 of 12
   for the reference year, 21 threshold ladders finite and increasing, every
   value finite in every year, 0 items to repair.
2. **The replication is a tail problem, not a replication problem.**
   **59 of 63** published per-diagnostic cells lie inside our 95 % interval —
   **21/21 at LOG, 20/21 at MOG, 18/21 at SOG** — and the ensemble reproduces
   Table 1 at every severity (+17/+28/+37/+47/+57 % against his
   +17/+28/+37/+46/+55 %). The level ratios are near 1 in the bulk and degrade
   into the tail: the stencil group's median runs 0.93 → 0.85 → 0.83 and the
   single-level group's 1.00 → 0.96 → 0.96. **The econometrics wants the tail,
   which is where the agreement is weakest** (§3 L15).
3. **13 indicators need no caveat, 2 need a stated caveat, and 6 are ensemble
   members only:** `negative_richardson`, `colson_panofsky`, `f2d`, `ubf`,
   `ncsu1` and — moved by the batch — `endlich`. Dropping all six moves the
   ensemble trend by **−0.3 pp (LOG), +0.1 pp (MOG), +0.4 pp (SOG)**.
4. **`ncsu1` is explained.** The prediction registered before the run held: at
   LOG it agrees with Prosser (+20 % against +21 %, level ratio 1.23, same
   significance call) and departs only above it, where his tail is dominated by
   the `1/max(Ri,1e-5)` floor that a finer vertical stencil triggers and ours
   does not (§6 Q3). This closes a question open since 2026-09-09.
5. **`ubf` is now the one real defect, and it is ours.** Its reference-year tail
   is polar — density 4.7 and 3.7 in the two caps at SOG, 56 % of the tail mass
   poleward of 60° — while its 30–60°N band drains from 2.0 to 1.0. Every
   jet-aligned diagnostic does the opposite. Ensemble member only (§3 L10).
6. **The zonal profiles reproduce the only spatial statements the literature
   makes.** Williams & Storer (2022) describe −Ri and Colson–Panofsky as having
   a single tropical peak, TI1 two midlatitude peaks, and frontogenesis no clear
   midlatitude peaks. Ours: −Ri 1.9/2.0 in the tropics and 0.01 at the poles,
   CP 1.9/1.9, `f2d` 1.7/2.0, TI1 1.74/1.79 in the two midlatitude bands. That
   is an independent check of spatial structure that no North Atlantic box
   average can give.
7. **The referee question to prepare for is Prosser's own argument:** finer
   vertical differences are "more accurate" (Prosser p. 5, against Lee et al.
   2023). The answer is in §6 Q5. One part of it cannot be closed: there is
   no 18-hPa series, so the trend bound covers all the differences from him
   at once, not the stencil on its own.

---

## 1. The selection

```
DATASET      /scistor/SBE-EDS-ClimateKoopman/yen230/derived/north_atlantic_audit/
             504 zarr stores  diagnostics_na_YYYY-MM.zarr, 1979-01 .. 2020-12
             21 variables, float32, 200 hPa, 0.25 deg, 3-hourly
             domain 30-60N, 75W-0; analysis box 36-60N, 55-10W (Prosser Fig. 2)
             .zattrs: audit_fixes            = meridional_metric+four_variable_provenance
                      f2d_variant            = A
                      deformation_convention = DEF (un-squared, Sharman A17)
                      target_level_hPa       = 200

THRESHOLDS   /scistor/SBE-EDS-ClimateKoopman/yen230/calibration/thresholds_2026-09-11.json
             provenance.period = "derived/global_audit (12 stores)"
             cos(phi)-weighted percentiles 97.0/99.1/99.6/99.8/99.9 of the full
             contiguous global year 2000 under the same conventions; all 21 signs "+"

REFERENCE    derived/global_audit/  (12 stores)   the year the thresholds came from
TAILS        calibration/tails_2026-09-11/        reference year only, value + latitude row

PAIRING      Score audit stores only with thresholds whose provenance period names
             derived/global_audit. full_trend_check.py, per_diagnostic_trend.py and
             per_diagnostic_all_severities.py refuse anything else (exit 2).
             The new driver also checks the first store's own audit_fixes attribute.

SENSITIVITY  derived/north_atlantic/ + thresholds_2026-09-07.json (baseline).
             For the 10 diagnostics that four-variable provenance changes, the
             baseline columns ARE the archived-field (ERA5 vo/d/pv) version.
             This is the ready answer to "your vorticity is worse physics" (below).

NOT FOR USE  derived/global_sub48, raw/global_sub48 (the 48-day comparison set),
             thresholds_2026-08-29/-09-03 (superseded), any mixed pairing.
```

**Why the evidence is enough.** You chose this dataset for three reasons.
Each one needs a different kind of support.

- **The choice affects levels, not trends, and trends are what the
  econometrics will estimate.** Moving from the baseline to the audit
  conventions changes 16 of the 21 fields and moves thresholds by up to 8 %.
  Every one of the 25 Table-1 trend ratios moves by 0.05 or less, and no
  significance call changes. On the level side, Prosser is the only external
  anchor, and the audit series follows his four-field method (p. 2). The
  choice therefore rests on method fidelity, and the levels are simply
  reported.
- **Two of the three changes are corrections.** F1 (the metric) is four
  orders of magnitude closer to exact ellipsoidal geometry. F6 (the `f2d`
  variant) follows the only published equation (W&S 2022 Eq. 3). Nothing
  argues for keeping the baseline on either of them.
- **The known cost of the third change can be measured without a third
  series.** Four-variable provenance (F12) replaces ERA5's spectral ζ, δ and
  PV with lower-fidelity finite differences. The baseline series holds the
  archived-field version of exactly those 10 diagnostics. It differs only
  through F1, whose flip rate is below 0.8 %. So the sensitivity is already
  on disk.

**The strongest argument against.** The econometric object is turbulence in
the atmosphere, not Prosser's indicator. On that view the audit series
knowingly uses the worse ζ, δ and PV. The IFS carries ζ and D as prognostic
spectral variables, and archived PV comes from L137. So for the vorticity
family, the baseline columns are arguably the better measurement.

A second, weaker point: the claim that "levels improved at every severity"
probably rests on one diagnostic. The baseline MOG ensemble deficit
attributable to `magnitude_pv` is −0.037. The audit improvement is +0.039.
Prediction P3 (§4.3) tests this. If it holds, the correct sentence becomes:
"switching PV to Prosser's A18 form moved one diagnostic's level onto his;
the rest barely moved".

**What would change my mind:**
1. The batch's regression check fails, meaning the new MOG table ≠ the
   job-1108148 table. Then something differs between two runs on the same
   stores. Stop, and use no audit per-diagnostic number until that is
   explained.
2. The sweep finds stores that were not written under the audit fix set.
   That calls for a repair, not a new choice:
   `sbatch --array=<idx>%8 jobs/25_diagnostics_na_audit.sbatch` recomputes
   only those months, at about 8 min each.
3. The econometric design treats vorticity-family magnitudes as physical
   quantities (EVT scale parameters in s⁻¹, for example). Then make the
   baseline columns primary for those 10, and report the audit columns as
   the sensitivity. That is a per-column swap, not a new dataset.

---

## 2. Corrections to the record

| # | The record says | What is true | Consequence |
|---|---|---|---|
| **C1** | `f2d` is damped ~36 % (sin(Wh)/(Wh) = 0.637) | For a diurnal signal with h = 3 h, Wh = π/4, so the factor is **0.900 (10 %)**. 0.637 is the **semi-diurnal** figure (or diurnal with h = 6 h). It applies only to the ∂Q/∂t term, not to the advective terms. Prosser also samples 3-hourly (p. 2: "every three hours"). | This is a magnitude-versus-truth caveat only. It is **not** a difference from Prosser and does not affect trends. Fix it in DECISION_RECORD §5.2, FORMULAS D3, STATUS §4d and essentials/DECISIONS §7. |
| **C2** | (not recorded) | Each month is computed from its own GRIB file, and `chunk_stitch.py` was never used (essentials/DECISIONS §9). So ∂/∂t at the first and last step of every month is either one-sided or NaN, which is **24 of 2920 steps a year (0.82 %)**. The record disagrees with itself on which: STATUS §15.3 says finite, while the `tail_thresholds.py` docstring says `f2d` "loses the two timesteps at each file's true ends". A one-sided difference passes a diurnal signal at 0.974 instead of 0.900. | **Settled by the census (P6): the edge steps are FINITE, i.e. one-sided.** STATUS §15.3 was right and the `tail_thresholds.py` docstring is wrong. ≤ 0.82 % of steps carry a first-order difference, the same count every year, so no trend effect. |
| **C3** | 17/21 significant = 17/21, "the single strongest result" | That is a **count** match. The two sets differ on two diagnostics: `colson_panofsky` (ours n.s., his p = 0.002) and `ncsu1` (ours t ≈ 4, his p = 0.7). The significance calls agree on **19 of 21**. | Report it as 19/21 agreement, with the two named. |
| **C4** | "The residual 10–16 % level deficit is the vertical stencil" | **True on the final dataset, and NOT on the baseline — the split is a property of the conventions.** Level-weighted ensemble gap, stencil group (11 diagnostics under A18) versus single-level (10): **LOG −0.059 / −0.039 (60 % stencil), MOG −0.083 / −0.039 (68 %), SOG −0.110 / −0.045 (71 %)**. On the *baseline* series the same split at MOG is −0.073 / −0.088, i.e. 45 %, because there `magnitude_pv` is archived Ertel PV on a single level and carries −0.037 of the gap on its own. Median level ratios, audit: 0.93/0.85/0.83 (stencil) against 1.00/0.96/0.96 (single-level). | Quote the split **with its severity and its series**. The earlier even split was measured on the baseline and does not describe the dataset. What the stencil still does not explain is `ubf`, which is single-level and is the largest single contributor at every severity (−0.030 at SOG). |
| **C5** | W&S (2022) Table 1 shows `endlich` is robust to vertical coarsening, so the stencil cannot explain its 0.51 | Table 1 is the **inter-dataset standard deviation** of LOG probability across HadGEM2-ES, ERA-Interim on the model grid and native ERA-Interim, averaged globally. It is not a test of vertical stencil width. FORMULAS D10 also says `endlich` is "not in the stencil-sensitive group", but it differentiates across 175/225 and `prosser_published.py` puts it in that group. | `endlich` at 0.51 sits inside its group's spread (−Ri 0.49, `f2d` 0.43, `ngm2` 0.61). There is no evidence that its level gap is anything other than the stencil (§3 L6). |
| **C6** | F5: "the Ri-floor candidate is DEAD" | F5 measured how often the floor binds on **our** stencil (0.0356 % of cells). That cannot rule out the floor's role on **his** stencil. | The `ncsu1` explanation is back on the table (§6 Q3). |
| **C7** | Scorecard: Figure S4-SOG is "figure only" | S4-LOG and S4-SOG print `rel`, `abs` and `p` in every panel. There are 42 exact numbers, now transcribed in `ada/prosser_published_s4.py`. The read-off levels average to within 0.2–0.4 % of Table 1. | Per-diagnostic comparison is now exact at three severities. |
| **C8** | `.gitignore` excludes `.pixi/` | The line is `.pixi/          # pixi builds...`, and git does not strip trailing comments, so the pattern never matched (tested today). On ADA, `git add -A` would stage the 1.2 GB environment. | Fixed in the repo: the comment now sits on its own line. |
| **C9** | FORMULAS D10 / §7: "`endlich` OPEN, needs a decision and one A/B run" | DECISION_RECORD §2 made the decision (angle secant, the literal reading). The A/B run on the full series does not fit in this phase (§3 L6). | Change FORMULAS D10 to "decided; sensitivity not measured on the series". |
| **C10** | ERA5.1: "behaviour in most of the troposphere is similar" | ECMWF says ERA5.1 gives "better global-mean temperatures in the **stratosphere and uppermost troposphere**". Poleward of the jet, 200 hPa lies in that layer, and 2000 is the calibration year. | This is not a clearance for 200 hPa. It is bounded instead (§3 L2). |
| **C11** | The stencil "changes the fitted relative trend by at most one percentage point" | That bound is the **recalibration** test (48-day vs full year, §15.11). No 18-hPa series exists, so the stencil's effect on the trend has never been isolated. What is measured is ours vs Prosser: annual ratios 0.97–1.04, which bounds **all** differences together (stencil, ERA5.1, grid, code). | Use: "our annual trends are within 4 % of Prosser's at every severity, which bounds the combined effect of every difference between the two pipelines, the stencil included". |

---

## 3. Limitations, with numbers

*"Cannot"* means no action inside this project's data access could close the
item. *"Chose not to"* means it could be closed at the stated cost, and that
cost was declined.

| | What | Why it stays open | Measured impact | Affects | Could it be closed? |
|---|---|---|---|---|---|
| **L1** | **Vertical stencil** 175/200/225 hPa vs Prosser's model levels 73–75 (≈188/197/206) | The CDS pressure-level product has nothing between 175 and 225 | Ensemble level: −9.6 / −11.2 / −12.0 / −13.5 / −15.8 % (LOG→SOG, annual, audit). The gap deepens with severity in DJF/MAM/SON and is flat in JJA. Stencil-group median level ratio **0.93 / 0.85 / 0.83** (LOG/MOG/SOG, audit) against **1.00 / 0.96 / 0.96** for the 10 single-level diagnostics; a quarter of the SOG level gap sits in diagnostics the stencil cannot touch. Trend: **combined** bound only, ≤ 4 % relative (C11) | Level, strongly. Trend within the combined bound. Per-diagnostic magnitudes | **Chose not to.** MARS "ERA5 complete" model levels, same volume as held; weeks of tape queue. Violates constraint 1 now |
| **L2** | **ERA5 instead of ERA5.1** for 2000–2006, including the calibration year | ERA5.1 is on MARS only | Block test (2000–06 dummy, baseline): ensemble anomaly **−5.5 % (t −1.23) MOG, −7.1 % (t −1.24) SOG**. Dummying the block out moves the ensemble change **+36.2 → +37.1 % (MOG)** and **+54.6 → +55.6 % (SOG)**. By leverage, a 2000–06 bias of δ moves the fitted change by 0.163 δ. Diagnostics with \|t\| ≥ 2: `magnitude_pv` (−17 %, −22 %) and `wind_speed` (−32 %, −46 %). The wind-speed anomaly shows that variability produces blocks of this size, so nothing here can be attributed to ERA5 | Trend ≤ 1 pp (ensemble). Level of N²-dependent diagnostics via the calibration year: not quantified | **Cannot** on CDS. Chose not to via MARS |
| **L3** | **`f2d` diurnal damping** | Property of 3-hourly sampling | Factor **0.900** on the ∂Q/∂t term for a diurnal cycle, 0.637 for a semi-diurnal one (C1). Identical in Prosser | Magnitude vs truth only | **Cannot** without hourly data (a new download) |
| **L4** | **`f2d` at month boundaries** (C2) | Months are processed per file | **Measured: the edges are one-sided, not NaN** (finite fraction 1.000 at the first and last step of every month). 24 of ~2920 steps a year, 0.82 %, carry a first-order difference, which passes a diurnal signal at 0.974 instead of 0.900 | Level ≤ ~1 %; trend none (constant count) | Chose not to. A re-derive with `chunk_stitch.py` is a new series, which constraint 2 rules out |
| **L5** | **60°N edge row.** The box's northern row is the edge of our download; Prosser downloaded the whole globe | Box and domain share their top edge | The row carries **0.78 %** of the box's cos-weight. Its horizontal derivatives use a one-sided stencil (assumed from MetPy; the census shows NaN if not) | Level < 1 %; trend none | Chose not to. It needs a download past 60.25°N (constraint 1) |
| **L6** | **`endlich`**: level ratio 0.51 (MOG) / 0.46 (SOG); ∂ψ/∂z reading | The literature says only "centred second-order differences". The two readings flip **41–66 %** of the exceedance set | Final values (convention-invariant): MOG +38 % [+23, +52] vs his +27.7 % (inside); **SOG +58 % [+37, +78] vs his +36.3 % (outside)**. Dropping `endlich` moves the ensemble trend by −0.1 pp | One diagnostic's level and SOG trend | **Chose not to.** A closed-form `endlich` series means 504 NA + 12 global months with a u/v-only reader (new code) plus recalibration: about 8–10 CPU-h, ~2 h wall. Over budget |
| **L7** | **`ncsu1`** level ratio 1.23 (LOG) / 6.38 (MOG) / 2.57 (SOG) | **Explained** — §6 Q3, prediction P5 confirmed | **At LOG the two agree**: +20 % [+10, +30] against his +21 %, same significance call. Above LOG: MOG +39 % [+20, +57] vs +4 %, SOG +47 % [+21, +73] vs +10 %, both outside. NA share of the global tail, his 0.92 → 0.18 → 0.36 against ours 1.13 → 1.12 → 0.93 | One diagnostic, above LOG only | Only with model levels (L1). The mechanism is his stencil, not our code |
| **L8** | **`brown2` dimensional gap** | Sharman A14 gives s⁻³; the tables print m² s⁻³. No length scale is published | **Exactly zero** effect on any exceedance decision, because a constant factor cancels in a percentile. Magnitudes are off by an unknown constant L² | Magnitude comparison and GPD scale σ only; ξ unaffected | **Cannot.** Sharman & Pearson (2017) Part I might name L; that is provenance, not a fix |
| **L9** | **The four non-significant diagnostics** (audit MOG): `horizontal_divergence`, `colson_panofsky`, `f2d`, `negative_richardson` | Three are properties of the field; one is statistical power | `horizontal_divergence` (+5 %), `f2d`, `negative_richardson`: **Prosser is also non-significant** and his rel lies inside our CI. `colson_panofsky`: ours +62 % [−6, +130] vs his +45 %, inside, but our NA exceedance is **7× sparser** (NA share 0.10 vs 0.71), so t = 1.85 is a power shortfall, not a disagreement | Significance calls only | Cannot for the first three. `colson_panofsky` is tied to L1 |
| **L10** | **`ubf`** level ratio 0.39 / 0.19 / **0.08**; SOG trend | **Structural, and ours** (P8). Its reference-year tail is polar: density **4.67 (90–60°S) and 3.70 (60–90°N) at SOG**, 56 % of the tail mass poleward of 60°, while the 30–60°N band drains 2.01 → 1.75 → 1.01 and NH/SH falls 1.59 → 0.99. Every jet diagnostic moves the other way (`vertical_wind_shear` 30–60°N 1.60 → 2.51, poles ≤ 0.41). W&S 2022 describe the published UBF as *stronger in the Northern Hemisphere* | LOG +13 % vs +18 % (inside), MOG +20 % vs +32 % (inside), **SOG +18 % [−6, +42] vs +44 %: outside**, significance differs. Largest single contributor to the ensemble level gap at every severity. Ensemble effect on the trend ≤ 0.1 pp | One diagnostic: level everywhere, SOG trend, and its magnitudes | **Cannot** without Prosser's code, but the next step is ours: UBF is a residual of near-cancelling terms and our spherical geometry was verified at mid-latitudes only (STATUS §4g). Repeat the Stokes and divergence-theorem checks on a high-latitude circuit. Note `magnitude_pv` and `temperature_gradient` are also polar-weighted for formula reasons (ζ+f, the polar front), so a polar tail is not by itself a defect — the draining midlatitude band is |
| **L11** | **Global and USA-box figures** (Prosser Figs 1, 2, 3b, S1, S2, S5) | They need new downloads | None on the NA dataset. **No claim outside 36–60N, 55–10W is supported.** Step 1b adds a zonal profile of each diagnostic's tail in the reference year, the only spatial check that fits the constraints | Generalisation only | **Chose not to.** Global 42 yr: ~1,500 CDS requests (3–4 days) + ~10 days of compute. USA box: ~470 GB, ~1.5 days download + ~7 h. Both violate constraint 1 |
| **L12** | **ERA5 observing-system changes** (STATUS §5 risk 3) | No independent reference at 200 hPa | Not measured. `horizontal_divergence` is flat in both studies. L2's block test covers 2000–06 only | Trend, in principle | Cannot within this dataset |
| **L13** | **No NA peaks-over-threshold archive**; tails exist for the reference year only | Its shape depends on the unit of observation | — | EVT readiness, not the indicators | Chose not to *yet*: one pass, ~23 GB, once the unit is fixed |
| **L14** | **No backup** on the share (§16.3) | Not asked | Loss costs compute (~7 h NA + ~14 h global + ~2 h calibration), not information | — | One email to ITvO before anything is published |

| **L15** | **The agreement is severity-dependent** | A property of the comparison, now measured at three severities | His value inside our CI: **21/21 (LOG), 20/21 (MOG), 18/21 (SOG)**. Ensemble level ratio 0.90 → 0.88 → 0.84. Stencil-group median ratio 0.93 → 0.85 → 0.83 | **Everything the econometrics does in the tail.** EVT and peaks-over-threshold sit at the severities where our agreement with the literature is weakest, and where four of the six class-C diagnostics do their damage | Cannot: it is L1 seen through the severity ladder. It is a disclosure, and it argues for reporting LOG and MOG alongside any SOG result |

Three limitations are new since the first version:
- **L5**, the edge row.
- **L10**, `ubf`'s polar tail — the batch's most consequential finding.
- **L15**, the severity dependence, which is the one that bears directly on
  the econometric design.

---

## 4. Additional checks

### 4.1 The list, and what each one decides

| | Check | Decision it changes | Cost | Outcome A → | Outcome B → | Run? |
|---|---|---|---|---|---|---|
| D1 | **Restore the baseline per-diagnostic table** from `data_prosser/` (intact copy from 2026-09-08, before the overwrite) | Whether the appendix sensitivity table can be cited | Seconds, login node | Step 3's 5-row identity passes → the copy is the 09-07 baseline; restore it | Identity fails → do not restore; regenerate (72 min, next phase) | **Yes** (desk) |
| D2 | **Read the F10 latitude-tilt table** already in the job-1102936 log | The wording of the `endlich` caveat | Seconds | F10 ran on the **global** file and the closed form *raises* 30–60°N exceedance → the alternative reading would move `endlich` toward Prosser; caveat reads "level depends on the reading" | Tilt ≤ 0 there, or F10 ran on the NA file (no global contrast) → caveat reads "level gap consistent with the stencil" | **Yes** (desk) |
| D3 | S4 transcription, baseline scorecard at MOG/SOG, ensemble leave-out, ERA5.1 block test, level-gap decomposition | §2 C3–C7, §3, §5 | Done today | — | — | **Done** |
| B1 | **Sweep of 504 + 12 stores, thresholds and logs** (`ada/final_sweep.py`) | Certify the dataset, or produce a repair list | ~5–15 min | 0 items → certified | N items → recompute those months with jobs/25 `--skip-if-matching` | **Yes** |
| B1b | **Zonal profile of each tail** in the reference year (`ada/final_tail_latitude.py`) | `ubf` caveat: regional vs structural. Consistency with W&S 2022's four spatial statements | ~5–10 min | `ubf` tail NH-weighted (NH/SH > 1, 30–60N density > 1) → the NA deficit is regional | NH/SH ≤ 1 or tropics densest → ours differs structurally from the Williams-lineage UBF; say so in the paper | **Yes** |
| B2 | **Per-diagnostic fits at all 5 severities in one pass, plus a finite census** (`ada/per_diagnostic_all_severities.py`) | §5 LOG/SOG columns for the 16 changed diagnostics; `ncsu1` P5; silent-failure sweep | ~40–75 min (est.), checkpointed | See P1–P8 | See P1–P8 | **Yes** |
| B3 | **Scorecard** (`ada/final_scorecard.py`): regression, identity, exact comparison at 3 severities, ensemble sensitivity, block test | Everything in §5 marked "pending" | Seconds | — | — | **Yes** |
| — | Per-diagnostic **DJF** fits | Nothing. Prosser publishes no seasonal per-diagnostic numbers, so a DJF run can only report significance calls on n = 42, which have flipped on rounding twice | 20 min | — | — | **No** |
| — | Separate `jobs/21` runs at light and severe | Superseded by B2, which is ~2× cheaper and adds a regression test | 144 min | — | — | **No** |
| — | Regenerate the baseline per-diagnostic CSV | Replaced by D1 | 72 min | — | — | **No** |
| — | `essentials/jobs/0` (bit-identity of the re-implementation on ADA) | Whether `essentials/` can be cited as the release code. It says nothing about the dataset | ~1 h (3 h requested) | — | — | **No.** Run it when the notebook work starts |
| — | F10 `endlich` series; F7 stencil-width experiment; USA / global | — | Over budget / needs downloads | — | — | **No** (L6, L1, L11) |

**Run 2026-09-21 as job 1124312: 32 min 34 s, exit 0**, all four steps clean,
all six file hashes matching. Sweep 8 min (0 items), tail profile 2 min, the
42-year fit 21 min at ~30 s a year, scorecard seconds. Well inside the 2 h 45
cap and the 3 h budget. D1 and D2 are the two desk checks; D1 is resolved by
P2 below.

### 4.2 How to run it

```bash
# Windows
cd data_turbulence && git status          # expect the new files listed below, nothing else
git add FINAL_DATASET.md .gitignore jobs/27_final_checks.sbatch ada/final_sweep.py \
        ada/final_tail_latitude.py ada/per_diagnostic_all_severities.py \
        ada/final_scorecard.py ada/prosser_published_s4.py \
        data_prosser/SCORECARD_baseline_S4.*
git commit -m "Final dataset selection and the one closing batch" && git push

# ADA, login node
cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence && git pull
grep -l FINAL-CHECKS-2026-09-16 jobs/27_final_checks.sbatch ada/final_sweep.py \
     ada/final_tail_latitude.py ada/per_diagnostic_all_severities.py \
     ada/final_scorecard.py ada/prosser_published_s4.py        # must list all 6
sha256sum jobs/27_final_checks.sbatch ada/final_*.py ada/per_diagnostic_all_severities.py ada/prosser_published_s4.py
#   ae1973789166…  jobs/27_final_checks.sbatch
#   7556777bef3f…  ada/final_sweep.py
#   670c81b6e544…  ada/final_tail_latitude.py
#   ce1504c438ab…  ada/final_scorecard.py
#   5ad7ffd71ecb…  ada/per_diagnostic_all_severities.py
#   60d20690e097…  ada/prosser_published_s4.py

# D1: restore the baseline table (the batch's step 3 then verifies it)
git ls-files --error-unmatch data_prosser/per_diagnostic_annual_moderate.csv \
  && cp data_prosser/per_diagnostic_annual_moderate.csv \
        data_prosser/per_diagnostic_annual_moderate_series.csv cat_outputs/
# D2: the F10 tilt
grep -l endlich_component_shear logs/*.out          # then, in that file:
grep -n -E "^input|LATITUDE TILT" -A12 <that log> | grep -E "input|band|endlich"

squeue -u yen230                                     # nothing else reading the audit stores
mkdir -p logs && sbatch jobs/27_final_checks.sbatch
# at 60 s: head -40 logs/final-<jobid>.out  -> git rev, 6 sha256, sweep started
# afterwards:
git add cat_outputs/final/*.csv cat_outputs/final/*.md logs/final-<jobid>.out && git commit -m "final checks" && git push
```

All four new scripts were run locally, before delivery, on synthetic zarr
stores that include NaN holes, month-edge NaNs and a whole-NaN month:
- **B2** was checked against `per_diagnostic_trend.py` at three severities,
  including the ensemble column. The outputs are **byte-identical**.
- **Resuming** from a checkpoint gives identical output.
- **B2 refuses a mismatched pairing**, and refuses a checkpoint written under
  a different configuration.
- **B1** flags the synthetic defects:
  - a failed log
  - an all-NaN month, which also shows up as missing chunk files (zarr 3
    does not write fill-value chunks)
  - scattered NaNs
  - a wrong grid
  - an unreadable store, which it reports instead of crashing
- **B1b** reproduces the global target fractions exactly.

### 4.3 Pre-registered predictions

Written before the batch runs. None is a pass/fail gate. Each maps its
outcome to an edit of §5.

| | Prediction | Outcome |
|---|---|---|
| **P1** | The new MOG table equals the job-1108148 audit MOG table to printed precision | **HOLDS.** Largest difference 0.00 × half-last-printed-digit. Two independent code paths, one answer |
| **P2** | The 5 convention-invariant rows are identical in the audit and baseline tables | **HOLDS in substance.** Identical at SOG and in 34 of 35 MOG cells; the exception is `endlich`'s CI upper bound at MOG, 0.522526 against 0.522522 — **4 × 10⁻⁶, i.e. 0.0004 percentage points**, which is float summation order between the xarray path and the numpy one, not a difference in the data. The restored baseline table is the 09-07 one, so D1 stands |
| **P3** | `magnitude_pv` MOG level ratio 0.61 → 0.9–1.1 | **PARTLY.** 0.61 → **0.87** at MOG, 0.72 at LOG, 0.96 at SOG. The direction and most of the size, just under the band. Its share of the ensemble level gap went from −0.037 (baseline MOG) to −0.011 |
| **P4** | `magnitude_pv`'s SOG trend moves toward his +90 % and his value falls inside our CI | **HOLDS.** +70 % [+46, +95] contains +90 %, against +43 % [+11, +75] on the baseline |
| **P5** | `ncsu1` at LOG: ratio within 0.75–1.33 and his +20.7 % inside our CI | **HOLDS.** Ratio 1.23, ours +20 % [+10, +30] against his +21 %, same significance call. The Ri-floor reading of §6 Q3 stands |
| **P6** | Finite census: nothing short except possibly `f2d`'s month edges | **HOLDS.** Every diagnostic finite in every cell of every year, and `f2d`'s first and last step of each month are finite — one-sided, not NaN (C2, L4) |
| **P7** | Sweep: 0 items; thresholds 21 increasing ladders, equal n | **HOLDS.** 504/504 and 12/12 under the audit fix set, n = 3,039,966,720 for every diagnostic, all signs `+`, 1,084 store writes in the logs with no failures |
| **P8** | `ubf`'s reference-year tail is not NH-weighted | **HOLDS, and sharper.** NH/SH falls to 0.99 at SOG and the tail is *polar* at both caps. L10 now reads structural |

**One prediction was wrong in its framing rather than its outcome.** P2's
tolerance was half the last printed digit, which is tighter than the arithmetic
noise between two code paths; the check reported DIFFERENT for a 0.0004 pp
disagreement. Widen the tolerance to one printed digit before reusing it.

**Rule for moving a diagnostic between classes after the batch** (applied by
hand, not in code):
- Move it only if a LOG or SOG row shows **both** of these at once:
  - Prosser's rel outside our 95 % CI
  - a level ratio outside 0.75–1.33
- A change in the significance call alone never moves a class.

---

## 5. Per-indicator verdict

**Classes**
- **A** — use without caveat for exceedance-based work.
- **B** — use, with the stated caveat.
- **C** — keep in the 21-member ensemble (needed for comparability with
  Prosser), but do **not** use as a stand-alone indicator.

**Columns**
- **Level ratios** are our fitted 1979 level ÷ Prosser's read-off 1979 level,
  at LOG / MOG / SOG. His levels are read off his figures and are good to a few
  per cent in aggregate, so read 0.95 and 1.05 as the same number.
- **Trend** is our fitted 1979→2020 change against his printed one, same order,
  and `yyy` says at which severities his value lies inside our 95 % interval.
- **■** marks the five diagnostics that are identical under both convention
  sets (P2), so their rows are also the baseline's.
- Every value is measured on the final dataset, job 1124312.

| W&J # | Diagnostic | Class | Level ratio LOG / MOG / SOG | Trend ours vs his (inside our CI?) | Flag |
|---|---|---|---|---|---|
| 3 | `brown1` | **A** | 1.01 / 1.10 / 1.00 | +18 / +44 / +73 vs +18 / +45 / +71 (yyy) | — |
| 4 | `temperature_gradient` | **A** | 0.93 / 0.91 / 0.91 | +7 / +27 / +41 vs +7 / +26 / +40 (yyy) | — |
| 5 | `horizontal_divergence` | **A** | 0.97 / 0.96 / 0.94 | +3 / +5 / +6 vs +3 / +7 / +8 (yyy) | **No trend, in either study** (+5 % vs +7 % at MOG). A property of the field, not a defect |
| 6 | `vertical_wind_shear` ■ | **A** | 0.89 / 0.85 / 0.89 | +13 / +39 / +60 vs +16 / +38 / +55 (yyy) | — |
| 8 | `deformation` | **A** | 1.02 / 1.01 / 0.97 | +17 / +37 / +55 vs +18 / +38 / +54 (yyy) | — |
| 9 | `wind_speed` ■ | **A** | 1.04 / 0.94 / 1.01 | +33 / +74 / +101 vs +32 / +76 / +102 (yyy) | Sparse at SOG (0.05 %), so the interval is wide |
| 12 | `rva_magnitude` | **A** | 1.01 / 0.97 / 0.95 | +20 / +29 / +34 vs +19 / +29 / +35 (yyy) | — |
| 14 | `nva` | **A** | 1.02 / 1.05 / 0.97 | +17 / +29 / +35 vs +16 / +28 / +36 (yyy) | Clipped at 0: a mass of exact zeros, so a POT threshold must sit above it |
| 15 | `brown2` | **A** | 0.94 / 0.92 / 0.88 | +14 / +46 / +76 vs +17 / +46 / +69 (yyy) | Magnitude only: A14 lacks an unstated L², so its units are not the published J kg⁻¹ s⁻¹. Ranks and exceedances are unaffected |
| 16 | `vorticity_squared` | **A** | 0.98 / 1.19 / 1.29 | +18 / +39 / +73 vs +19 / +37 / +69 (yyy) | — |
| 17 | `ti1` | **A** | 0.93 / 0.88 / 0.83 | +15 / +44 / +69 vs +17 / +42 / +58 (yyy) | — |
| 18 | `ngm1` | **A** | 1.00 / 0.94 / 0.88 | +28 / +52 / +79 vs +27 / +51 / +78 (yyy) | — |
| 19 | `ti2` | **A** | 0.94 / 0.87 / 0.84 | +11 / +29 / +41 vs +13 / +27 / +32 (yyy) | — |
| 1 | `magnitude_pv` | **B** | 0.72 / 0.87 / 0.96 | +28 / +42 / +70 vs +22 / +51 / +90 (yyy) | The quantity is Sharman A18 truncated PV on the 50 hPa stencil, in hPa units (×100 vs SI, not PVU) — **not** ERA5 Ertel PV. Never quote a PVU magnitude. Level and trend now agree with Prosser at every severity |
| 10 | `ngm2` | **B** | 0.83 / 0.62 / 0.57 | +18 / +28 / +37 vs +24 / +34 / +42 (yyy) | Level 17–43 % below his (∂T/∂z on the 50 hPa stencil, plus a 2.8 % truncation error). Trends agree at all three severities |
| 2 | `colson_panofsky` ■ | **C** | 0.20 / 0.14 / 0.12 | +37 / +62 / +68 vs +33 / +45 / +55 (yyy) | Tail is **tropical** (1.9–2.1 in the tropics, 0.00 at the poles), so our box holds 5–8× fewer of its events than his. The trend agrees in size and fails significance on power, not on disagreement |
| 7 | `endlich` ■ | **C** | 0.67 / 0.51 / 0.46 | +18 / +38 / +58 vs +14 / +28 / +36 (yy**N**) | Level about half of his at every severity (0.67 / 0.51 / 0.46) and **SOG trend +58 % [+37, +78] excludes his +36 %**. ∂ψ/∂z uses the angle-secant reading; the closed form flips 41–66 % of the set and was never measured on the series. Usable with the level caveat at LOG and MOG |
| 11 | `negative_richardson` ■ | **C** | 1.53 / 0.49 / 0.26 | +44 / +41 / +18 vs +16 / -8 / -14 (yyy) | Tail is **tropical** (density 1.9 in 0–30°N, 0.01 at the poles), so the North Atlantic sees almost none of it: 0.019 % at MOG, 0.003 % at SOG. No trend in either study above LOG. Richardson family: the most stencil-sensitive |
| 13 | `ubf` | **C** | 0.39 / 0.19 / 0.08 | +13 / +20 / +18 vs +18 / +32 / +44 (yy**N**) | **Structural.** Its extreme tail is polar (density 4.7 / 3.7 in the two caps at SOG, 56 % of the tail mass poleward of 60°) while the 30–60°N band drains from 2.0 to 1.0. Every jet diagnostic does the opposite. Ensemble member only; do not use its magnitudes |
| 20 | `f2d` | **C** | 0.97 / 0.62 / 0.38 | +14 / +25 / +27 vs +14 / +0 / +0 (yyy) | Tail is **tropical**. Signed variant A, thresholds cross zero — read it from flip rates, never from relative threshold differences. ∂Q/∂t damped ×0.900 (diurnal); month-edge steps one-sided (×0.974) |
| 21 | `ncsu1` | **C** | 1.23 / 6.38 / 2.57 | +20 / +39 / +47 vs +21 / +4 / +10 (y**N****N**) | Agrees with Prosser at LOG (+20 % vs +21 %, ratio 1.23). Above LOG his `ncsu1` is dominated by the `1/max(Ri,1e-5)` floor on a stencil that resolves thin unstable layers and ours is not, so the two measure different populations. For EVT, note the floor would plant an artificial mass in the extreme tail |

**The econometric consequences.**
1. **Ensemble.** Use all 21, as Prosser does. Dropping the six class-C members
   moves the ensemble trend by **−0.3 pp (LOG), +0.1 pp (MOG), +0.4 pp (SOG)**.
   Leave-one-out spans 15.8–17.6 % (LOG), 36.3–38.5 % (MOG) and 55.1–59.1 %
   (SOG). `temperature_gradient` is at the high end at all three severities; the
   low end is `ngm1` at LOG and `brown2` at MOG and SOG. No single member drives
   the result.
2. **Per-diagnostic panels.** Use classes A and B. Report class-C results
   only as "ensemble members", each with its flag.
3. **Magnitudes (EVT).** Work only with A-class diagnostics and
   `magnitude_pv`/`brown2` under their unit caveats. `ncsu1`,
   `negative_richardson` and `f2d` have tail structure driven by a floor, a
   sign convention or damping.

---

## 6. The five questions

**Q1. Is `north_atlantic_audit` + `thresholds_2026-09-11.json` the right final
choice? What is the strongest argument against it?**

Yes (§1). The strongest argument against is that it deliberately uses
lower-fidelity ζ, δ and PV than ERA5 archives. For a study of the atmosphere
rather than of Prosser's indicator, the archived fields are closer to the
model state.

The answer to that argument is already on disk: for those 10 diagnostics,
the baseline series is the archived-field version. A second, smaller point:
the "levels improved" claim probably reduces to `magnitude_pv` (P3).

**Q2. Which indicators would you not hand to an econometrician without a
flag, and what is the flag?**

Eight, in two groups. The batch moved `endlich` from B to C on its SOG trend,
so the final split is **13 A / 2 B / 6 C** (§5).

- **Class C — ensemble members only, never stand-alone:** `colson_panofsky`,
  `negative_richardson`, `f2d`, `endlich`, `ubf`, `ncsu1`.
- **Class B — usable with the stated caveat:** `magnitude_pv` (the quantity is
  Sharman A18 truncated PV in hPa units, not ERA5 Ertel PV; never quote a PVU
  magnitude) and `ngm2` (1979 level 17–43 % below his; trends agree).

The six C's fail for three different reasons, and only one of them is a defect
in our pipeline:

1. **Wrong population, correctly computed** — `colson_panofsky`,
   `negative_richardson`, `f2d`. All three have **tropical** reference-year
   tails (B3), so a North Atlantic box holds 5–8× fewer of their events than
   the global calibration implies. The flag is "this diagnostic's extremes are
   not where your sample is", not "this number is wrong".
2. **Level disagreement with Prosser** — `endlich` (0.67 / 0.51 / 0.46) and
   `ubf` (0.39 / 0.19 / 0.08), both with his trend outside our CI at SOG.
   `ubf` is the one structural problem in the set: its tail is polar where
   every jet diagnostic's is midlatitude (L10).
3. **A different population above LOG** — `ncsu1` (Q3).

Separately, four class-A diagnostics carry a **magnitude-only** flag and need
no flag at all for exceedance work: `brown2` (A14 lacks an unstated L², so the
units are not the published J kg⁻¹ s⁻¹), `deformation` (DEF, not DEF²), `nva`
(a mass of exact zeros, so a peaks-over-threshold cut must sit above it) and
`wind_speed` (0.05 % of cells at SOG, so its interval is wide).

**Q3. Where does `ncsu1` stand, and is it usable?**

It is correctly computed. It does not reproduce Prosser's `ncsu1` above LOG.
It is usable as an ensemble member and not as a stand-alone indicator.

The explanation below fits all six published `ncsu1` numbers (level and
trend at LOG, MOG and SOG) and nothing on disk contradicts it:

| | LOG | MOG | SOG |
|---|---|---|---|
| Prosser: NA share of the global tail | **0.92** | **0.18** | **0.36** |
| Prosser: trend | **+20.7 %, p = 0.0009** | +4.3 %, p = 0.7 | +9.7 %, p = 0.5 |
| Ours (final dataset): NA share | **1.13** | **1.12** | **0.93** |
| Ours: trend | **+20 % [+10, +30]** | +39 %, sig | +47 %, sig |

Ours is flat in NA share across all three severities — about 1.1 — and his
collapses to 0.18 and then partly recovers. That is the whole puzzle in one
row.

A36 multiplies by 1/max(Ri, 10⁻⁵). Wherever Ri ≤ 10⁻⁵ (static instability),
the diagnostic is inflated by up to 10⁵, so those cells float to the top of
the global distribution. Suppose they make up about 0.3 % of global cells on
Prosser's 18-hPa stencil, almost all outside the NA box. Then:

- **LOG (top 3 %):** ordinary dynamic cells dominate, so his NA share is
  normal (0.92) and his trend looks like everyone else's (+20.7 %).
- **MOG (top 0.4 %):** mostly floor-bound cells, so his NA share collapses
  (0.18) and his trend becomes that of rare instability events, flat like
  his `negative_richardson`.
- **SOG (top 0.1 %):** among the floor cells, those with the largest dynamic
  product rank first. These favour the jets, so the NA share partly recovers
  (0.36).

No other diagnostic in his S4 shows that non-monotone pattern.

On our 50-hPa stencil, thin unstable layers are averaged away. The floor
binds on **0.0356 %** of cells (F5), below even the SOG fraction. So our
`ncsu1` is dynamics-driven at every severity: NA share ~1, trend +39 / +47 %.
The required rate on his stencil (~0.3 %) is about **nine times** ours.

This is a hypothesis about data we do not have, and it was pre-registered as
P5: at LOG, where his `ncsu1` should still be dynamics-driven, ours must agree
with his on both level and trend. **P5 held.** At LOG our level ratio is
**1.23** and our **+20 % [+10, +30]** contains his **+21 %**, with the same
significance call; above LOG the two diverge exactly as the floor argument
requires (MOG ratio 6.38, SOG 2.57; his rel outside our CI at both). This is
the first reading of `ncsu1` that fits all six of his published numbers and
survives a test it could have failed.

It does not make `ncsu1` usable on its own. The class stays **C**, and for any
extreme-value work note that on Prosser's stencil the floor would plant an
artificial mass in the extreme tail — ours does not have it, which is a reason
our `ncsu1` is *better* behaved and still not comparable to his.

**Q4. What is assumed rather than verified?**

Twelve things were on this list before the batch. **Seven are now measured**
and are struck through; five remain assumptions, and four of those cannot be
settled from inside this project.

*Settled by job 1124312:*

1. ~~The `audit_fixes` attribute and completeness of all 504 audit stores and
   all 12 `global_audit` stores.~~ **Verified.** 504/504 and 12/12 carry
   `meridional_metric+four_variable_provenance`, `f2d_variant = A`, the
   un-squared DEF string and `target_level_hPa = 200`; 0 stores with problems;
   1,084 store writes in the job logs, none failed.
2. ~~That the audit NA series is finite everywhere.~~ **Verified.** Every one
   of the 21 diagnostics is finite in every cell of every one of the 42 years.
   Not "no all-NaN variables" — no NaN at all.
3. ~~`f2d` month-edge handling.~~ **Verified one-sided, not NaN**: finite
   fraction 1.0 at the first and last step of every month (C2, L4). The record
   has been corrected; the damping is ×0.974 at those steps.
4. ~~The 60°N edge row.~~ **Verified** by the same census — one-sided and
   finite, as L5 assumed.
5. ~~That the job-1108148 audit MOG table is right.~~ **Verified.** The new
   table, computed through a different code path (numpy on loaded arrays
   rather than xarray reductions), is identical to printed precision — largest
   difference 0.00 × half the last printed digit (P1).
6. ~~That `data_prosser/` was scored on the 09-07 thresholds.~~ **Verified**
   (P2): the five convention-invariant diagnostics are identical at SOG and
   agree at MOG to 4 × 10⁻⁶, which is float summation order between the two
   code paths, not a difference in the data. D1 stands.
7. ~~That `thresholds_2026-09-11.json` has equal n and increasing ladders.~~
   **Verified.** n = 3,039,966,720 for all 21 (721 × 1440 × 2,928 steps, the
   full leap year 2000), all signs `+`, 21 strictly increasing finite ladders,
   and the provenance string names the audit reference year.

*Still assumptions:*

8. **Every per-diagnostic 1979 level of Prosser.** Read off his figures.
   Validated in aggregate against his exact Table 1 — mean of the 21 read-offs
   5.3381 % vs 5.3253 % at LOG (+0.2 %), 0.8009 % vs 0.7991 % at MOG (+0.2 %),
   0.2012 % vs 0.2021 % at SOG (−0.4 %) — which bounds the *systematic* error
   at well under 1 %, not the error in any one panel. This is why §5 reads 0.95
   and 1.05 as the same number, and why no verdict in §5 rests on a level ratio
   inside 0.75–1.33 alone.
9. **That 188/197/206 hPa are L137 levels 73–75.** Only Prosser can settle it,
   and the conclusion (his stencil is finer than the CDS product offers) does
   not depend on the reading.
10. **The size of the ERA5 vs ERA5.1 difference at 200 hPa.** Still bounded
    rather than measured, but the bound is now empirical as well as
    documentary: a 2000–2006 block dummy moves the 42-year change by
    **−14.3 pp to +10.4 pp across the 21 diagnostics, with |t| < 1.8
    everywhere and 18 of 21 within ±6 pp** (LOG; the pattern holds at MOG and
    SOG). No diagnostic has a significant 2000–06 block. That is not proof of
    no effect — it is a demonstration that the effect is smaller than this
    design can see, which is the honest claim (L2).
11. **That `essentials/` reproduces the ADA stores bit for bit.** Shown on
    local test files; `essentials/jobs/0` has never run on ADA. Matters only if
    `essentials/` becomes the published code.
12. **That rojak treats longitude as periodic** on the global grid, affecting
    2 of 1,440 columns of the calibration year. Below the percentile noise.

One item was *added* to this list by the batch: the six `ubf` zonal-density
numbers in L10 come from the reference year only. Whether the polar tail is a
property of the diagnostic or of the year 2000 specifically is not tested. It
would take a second calibration year, which constraint 2 rules out, and it does
not change the verdict — class C either way.

**Q5. What is a referee most likely to challenge, and what is the answer?**

*"Prosser et al. argue (p. 5) that their 188/206-hPa differences are 'much
finer (and therefore more accurate)', and use that to explain their
disagreement with Lee et al. (2023). You use 175/225. By his own argument
your diagnostics are less accurate. Why should your trends be trusted?"*

The answer has four parts:

1. **It is a constraint.** 175/200/225 is the finest stencil at 200 hPa in
   the ERA5 pressure-level product. Model levels are MARS-only and
   terrain-following. Wong et al. (2025) and Williams & Storer (2022) use
   pressure levels too.
2. **The bias it causes is in the level, and it is stable in time.**
   - The 11 stencil-using diagnostics have median level ratios of **0.93 (LOG),
     0.85 (MOG), 0.83 (SOG)** against **1.00 / 0.96 / 0.96** for the 10
     single-level ones. The gap is real, it is in the direction the physics
     predicts, and it grows with severity.
   - It is a **level** effect, not a trend effect: level ratios drift by
     ≤ 0.006 over 41 years, and the ensemble level gap decomposes as
     −0.059 stencil / −0.039 single-level at LOG, −0.083 / −0.039 at MOG and
     −0.110 / −0.045 at SOG. Even at SOG, **a quarter of the level gap is in
     diagnostics the stencil cannot touch** — so the stencil is not the only
     cause of the level deficit, and whatever else is at work is not
     severity-specific.
   - The deficit deepens with severity in DJF/MAM/SON and not in JJA, which is
     where a stencil that misses sharp jet-level shear should matter.
3. **Our trends match his where it counts.** **59 of 63** per-diagnostic
   published changes (21 diagnostics × 3 severities) lie inside our 95 %
   intervals — 21/21 at LOG, 20/21 at MOG, 18/21 at SOG — and we make the same
   significance call on 19/21, 19/21 and 18/21. The three that fall outside at
   SOG are `endlich`, `ubf` and `ncsu1`, all already class C. At the ensemble
   level our 42-year changes are +16.5 % (LOG), +37.2 % (MOG) and +57.3 % (SOG)
   and all 25 Table-1 cells contain his value. That bounds the combined effect
   of every pipeline difference, the stencil included. Prosser's own contrast
   with Lee et al. was about *where* hotspots appear (and Lee used a 100-hPa
   stencil); within his North Atlantic box our changes are his.
4. **Where the stencil changes behaviour, the diagnostic is flagged and the
   result does not depend on it.** That applies to the Richardson family and to
   `ncsu1` through its floor. Dropping all **six** class-C diagnostics moves
   the ensemble change by **−0.3 pp (LOG), +0.1 pp (MOG), +0.4 pp (SOG)**, and
   leave-one-out over all 21 spans 15.8–17.6 %, 36.3–38.5 % and 55.1–59.1 %.
   No single member, flagged or not, carries the result.

What we cannot say is how large the stencil effect on the trend is *on its
own*. Say that in the paper, in those words (C11).

---

## 7. What remains

The batch is run and its outcomes are in §4.3 and §5. Nothing on this list
needs the cluster except item 4, and item 4 is optional.

**Done, and closing the data-gathering phase.**

1. ~~Paste the consistency lines and the changed rows into §5.~~ Done: §5 is
   now generated from `cat_outputs/final/*.csv`, and §4.3 carries the eight
   outcomes.
2. **Corrections still to make in the other documents** (desk work, one
   sitting): `DECISION_RECORD.md` — C1, C4, C5, C9, C11;
   `FORMULAS_AND_DECISIONS.md` — D3 (damping is 0.900 diurnal, not 0.637),
   D10 (`endlich` is decided: class C on its SOG trend), and close the open
   items in its §7; `essentials/DECISIONS.md` §7 the same. Note that
   `ada/prosser_published_s4.py` gained the 42 `abs=` values after the batch
   was submitted, so its sha256 no longer matches the table in §4.2 — the
   version marker is unchanged, so the job guard still passes, but do not
   quote the old hash.
3. **The notebook.** `notebooks/final_dataset_report.ipynb` reproduces every
   figure and table in this document from `cat_outputs/final/`. It
   auto-discovers the final audit outputs, falls back to the baseline with a
   loud banner, and runs end to end. It is the reproducibility artefact; this
   file is the argument.

**Optional, and none of it changes a verdict.**

4. Digitise the 42 blue crosses of Prosser's Figure 3a for a year-by-year
   correlation against our ensemble series. Desk work, one afternoon, and the
   only check left that could still *strengthen* the replication claim rather
   than qualify it. Cost on the cluster: zero.
5. A |lat| ≤ 60 recalibration of `ubf` from the tail archive, to see whether
   its polar tail is the whole story. This would be a **third** set of
   thresholds for one diagnostic, which constraint 2 forbids as a dataset
   change; it could only ever be a footnote, and the verdict is C regardless.
   Recommendation: do not.
6. `essentials/jobs/0` on ADA, if and only if `essentials/` becomes the
   published code (Q4 item 11).

**Then stop.** Next come the unit of observation, then the North Atlantic
peaks-over-threshold archive (L13), then the ITvO email (L14).

---

## 8. What to keep, and what to delete

About **1.5 TB of the 2.5 TB quota** is in use. The dataset this document
certifies is a small part of it. The figures below are the `du` taken while
preparing this batch — **re-run `du -sh` on a directory before deleting it**,
because a delete on the scratch filesystem is not recoverable.

**Keep, permanently. This is the dataset.**

| What | Where | Why |
|---|---|---|
| The audit NA series, 504 stores | `derived/north_atlantic_audit/` | The dataset (§1) |
| The calibration | `calibration/thresholds_2026-09-11.json` | Without it the series cannot be scored |
| The tail archive, ~41 GB | `calibration/tails_2026-09-11/` | The only thing that can answer a new percentile question without a rerun; B3 and every zonal profile in §5 came from it |
| The 12 `global_audit` reference-year stores | `derived/global_audit/` | Provenance for the calibration; small |
| All of `cat_outputs/` and `jobs/` logs | — | The evidence trail; megabytes |

**Delete, in this order.** Stop when there is enough room; each step is
independent.

1. **`raw/global/` day-blocks, ~119 GB.** These were the input to the
   reference-year calibration, which is finished and archived. Reproducible
   from CDS if ever needed. *Delete first.*
2. **`derived/global/` (baseline-convention global year), ~196 GB.** Superseded
   by `derived/global_audit/`. The only thing it still supports is a baseline
   recalibration, which constraint 2 rules out.
3. **The `*_sub48` pair, ~42 GB.** Development subsets. Their job is done.
4. **The baseline NA series, ~294 GB.** *Delete last, and think before you do.*
   This is the robustness arm. Every claim in §6 Q1 — that the
   archived-field version exists and agrees — rests on it, and §4.3 P2/P3/P4
   compare against it. Keeping it costs 294 GB; regenerating it would cost
   about a week of cluster time and would breach constraint 2. **Recommendation:
   keep it until the paper is submitted**, then delete.

Deleting steps 1–3 frees about **357 GB** and touches nothing this document
depends on. That should be enough; step 4 should not be needed.

**Sources.** Prosser et al. (2023) pp. 2, 5, Fig. 4, SI Figs S4-LOG/SOG
(`Articles/`). Williams & Storer (2022) pp. 1428–1432, Tables 1–2. Sharman
et al. (2006) App. A. ECMWF on ERA5.1:
[NCAR GDEX summary](https://gdex.ucar.edu/news/era51-corrections-to-era5-stratospheric-temperature-2000-2006/),
[Simmons et al. 2020, ECMWF TM 859](https://www.ecmwf.int/en/elibrary/81149-global-stratospheric-temperature-bias-and-other-stratospheric-aspects-era5-and).
Everything else is from the project docs and the repository, and from job
1124312's archived outputs in `cat_outputs/final/`.
