# FINAL DATASET — CAT indicators, North Atlantic, 1979–2020

**Written 2026-09-16.** This closes the data-gathering phase. It does not
restate `DECISION_RECORD.md`. It checks that document's conclusion, corrects
eleven statements in the record that are wrong or overstated (§2), puts a
number on every limitation (§3), and gives each of the 21 indicators a class
(§5). The only open work is one batch of about 1–1.7 h, capped at 2 h 45 min
(§4). Its outcomes can change the wording of some caveats. They cannot change
the dataset.

**Evidence status.** Everything marked *baseline* below was computed today
from `data_prosser/`, the copy of the pre-audit per-diagnostic tables that
survived the 2026-09-15 overwrite. For the five convention-invariant
diagnostics, those rows are also the final rows exactly (§5). For the other
sixteen, the audit rows exist only on ADA and come back with the batch.

---

## 0. The answer

1. **Use `derived/north_atlantic_audit/` with `thresholds_2026-09-11.json`.**
   I agree with `DECISION_RECORD.md`. I disagree with three of its supporting
   claims (§2 C3, C4, C11). None of them changes the choice.
2. **13 indicators need no caveat, 3 need a stated caveat, and 5 should stay
   in the 21-member ensemble but should not be used on their own:**
   `negative_richardson`, `colson_panofsky`, `f2d`, `ubf`, `ncsu1`. On the
   baseline tables, dropping all five moves the ensemble trend by
   **+0.0 pp at MOG and +0.3 pp at SOG**.
3. **`ncsu1` now has a specific explanation, and it fits all six of Prosser's
   published `ncsu1` numbers.** The Ri floor was declared dead after it was
   tested on *our* data, where it almost never binds. It was never tested
   against *his* data. His North-Atlantic share of the `ncsu1` tail runs
   0.92 → 0.18 → 0.36 across LOG → MOG → SOG. No other diagnostic shows that
   collapse and partial recovery, and a multiplicative floor produces exactly
   that shape (§6 Q3). The batch runs one test that can refute this.
4. **The record gets the `f2d` damping wrong.** 0.637 is the figure for a
   semi-diurnal signal. For a diurnal signal it is **0.900**, and Prosser
   samples 3-hourly too, so this is not a difference from him (§2 C1).
5. **The referee question to prepare for is Prosser's own argument:** finer
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
| **C2** | (not recorded) | Each month is computed from its own GRIB file, and `chunk_stitch.py` was never used (essentials/DECISIONS §9). So ∂/∂t at the first and last step of every month is either one-sided or NaN, which is **24 of 2920 steps a year (0.82 %)**. The record disagrees with itself on which: STATUS §15.3 says finite, while the `tail_thresholds.py` docstring says `f2d` "loses the two timesteps at each file's true ends". A one-sided difference passes a diurnal signal at 0.974 instead of 0.900. | This affects the `f2d` level by less than 1 % of steps, with the same count every year, so no trend effect. The batch census settles finite or NaN (P6). |
| **C3** | 17/21 significant = 17/21, "the single strongest result" | That is a **count** match. The two sets differ on two diagnostics: `colson_panofsky` (ours n.s., his p = 0.002) and `ncsu1` (ours t ≈ 4, his p = 0.7). The significance calls agree on **19 of 21**. | Report it as 19/21 agreement, with the two named. |
| **C4** | "The residual 10–16 % level deficit is the vertical stencil" | Weighting by level, the baseline MOG ensemble deficit of −0.161 splits into **−0.073 from the 10 stencil diagnostics and −0.088 from the 11 single-level ones**. The largest single contributions are `magnitude_pv` −0.037 and `ubf` −0.035. `ncsu1` (+0.021) hides part of the stencil group's share. SOG: −0.107 / −0.095. The medians (0.73 vs 0.93) are still correct. | The stencil explains the **pattern** within the vertical-derivative group. It does not explain the ensemble deficit **on its own**: `ubf` is single-level and makes up about a fifth of it. Under the audit conventions, A18 puts `magnitude_pv` into the stencil group. |
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
| **L1** | **Vertical stencil** 175/200/225 hPa vs Prosser's model levels 73–75 (≈188/197/206) | The CDS pressure-level product has nothing between 175 and 225 | Ensemble level: −9.6 / −11.2 / −12.0 / −13.5 / −15.8 % (LOG→SOG, annual, audit). The gap deepens with severity in DJF/MAM/SON and is flat in JJA. Stencil-group median level ratio 0.73 vs 0.93 (baseline MOG). Trend: **combined** bound only, ≤ 4 % relative (C11) | Level, strongly. Trend within the combined bound. Per-diagnostic magnitudes | **Chose not to.** MARS "ERA5 complete" model levels, same volume as held; weeks of tape queue. Violates constraint 1 now |
| **L2** | **ERA5 instead of ERA5.1** for 2000–2006, including the calibration year | ERA5.1 is on MARS only | Block test (2000–06 dummy, baseline): ensemble anomaly **−5.5 % (t −1.23) MOG, −7.1 % (t −1.24) SOG**. Dummying the block out moves the ensemble change **+36.2 → +37.1 % (MOG)** and **+54.6 → +55.6 % (SOG)**. By leverage, a 2000–06 bias of δ moves the fitted change by 0.163 δ. Diagnostics with \|t\| ≥ 2: `magnitude_pv` (−17 %, −22 %) and `wind_speed` (−32 %, −46 %). The wind-speed anomaly shows that variability produces blocks of this size, so nothing here can be attributed to ERA5 | Trend ≤ 1 pp (ensemble). Level of N²-dependent diagnostics via the calibration year: not quantified | **Cannot** on CDS. Chose not to via MARS |
| **L3** | **`f2d` diurnal damping** | Property of 3-hourly sampling | Factor **0.900** on the ∂Q/∂t term for a diurnal cycle, 0.637 for a semi-diurnal one (C1). Identical in Prosser | Magnitude vs truth only | **Cannot** without hourly data (a new download) |
| **L4** | **`f2d` at month boundaries** (C2) | Months are processed per file | ≤ 0.82 % of steps are one-sided or NaN. If one-sided, the diurnal factor is 0.974 there | Level ≤ ~1 %; trend none (constant count) | Chose not to. A re-derive with `chunk_stitch.py` is a new series, which constraint 2 rules out |
| **L5** | **60°N edge row.** The box's northern row is the edge of our download; Prosser downloaded the whole globe | Box and domain share their top edge | The row carries **0.78 %** of the box's cos-weight. Its horizontal derivatives use a one-sided stencil (assumed from MetPy; the census shows NaN if not) | Level < 1 %; trend none | Chose not to. It needs a download past 60.25°N (constraint 1) |
| **L6** | **`endlich`**: level ratio 0.51 (MOG) / 0.46 (SOG); ∂ψ/∂z reading | The literature says only "centred second-order differences". The two readings flip **41–66 %** of the exceedance set | Final values (convention-invariant): MOG +38 % [+23, +52] vs his +27.7 % (inside); **SOG +58 % [+37, +78] vs his +36.3 % (outside)**. Dropping `endlich` moves the ensemble trend by −0.1 pp | One diagnostic's level and SOG trend | **Chose not to.** A closed-form `endlich` series means 504 NA + 12 global months with a u/v-only reader (new code) plus recalibration: about 8–10 CPU-h, ~2 h wall. Over budget |
| **L7** | **`ncsu1`** level ratio 6.07 (MOG) / 2.47 (SOG); his trend flat at MOG/SOG | See §6 Q3 | Ours MOG +39 % [+20, +58] vs his +4.3 %; SOG +47 % [+21, +73] vs +9.7 %: **outside at both**. Ensemble effect −0.1 / +0.2 pp | One diagnostic, level and trend | Only with model levels (L1) |
| **L8** | **`brown2` dimensional gap** | Sharman A14 gives s⁻³; the tables print m² s⁻³. No length scale is published | **Exactly zero** effect on any exceedance decision, because a constant factor cancels in a percentile. Magnitudes are off by an unknown constant L² | Magnitude comparison and GPD scale σ only; ξ unaffected | **Cannot.** Sharman & Pearson (2017) Part I might name L; that is provenance, not a fix |
| **L9** | **The four non-significant diagnostics** (audit MOG): `horizontal_divergence`, `colson_panofsky`, `f2d`, `negative_richardson` | Three are properties of the field; one is statistical power | `horizontal_divergence` (+5 %), `f2d`, `negative_richardson`: **Prosser is also non-significant** and his rel lies inside our CI. `colson_panofsky`: ours +62 % [−6, +130] vs his +45 %, inside, but our NA exceedance is **7× sparser** (NA share 0.10 vs 0.71), so t = 1.85 is a power shortfall, not a disagreement | Significance calls only | Cannot for the first three. `colson_panofsky` is tied to L1 |
| **L10** | **`ubf`** level ratio **0.19 (MOG) / 0.08 (SOG)**; SOG trend | Single-level, so not the stencil. W&S 2022 find UBF the **most** dataset-robust of their seven (Table 1: 0.20–0.27 %) and describe it as "stronger in the Northern Hemisphere". Our NA share is 0.35 / 0.10 against his 1.80 / 1.35 | MOG +20 % [+3, +37] vs his +31.6 % (inside). **SOG +18 % [−6, +42] vs +43.7 %: outside**, and the significance calls differ. Ensemble effect +0.1 pp. F12 moves `ubf` by only ~1 %, so the audit rows will be close to these | One diagnostic, level and SOG trend | **Cannot** without Prosser's code. Candidates: the geometry of ∇²Φ and J (the kind of error §4g found in rojak), the Φ convention. Step 1b (§4) decides whether the difference is regional or structural |
| **L11** | **Global and USA-box figures** (Prosser Figs 1, 2, 3b, S1, S2, S5) | They need new downloads | None on the NA dataset. **No claim outside 36–60N, 55–10W is supported.** Step 1b adds a zonal profile of each diagnostic's tail in the reference year, the only spatial check that fits the constraints | Generalisation only | **Chose not to.** Global 42 yr: ~1,500 CDS requests (3–4 days) + ~10 days of compute. USA box: ~470 GB, ~1.5 days download + ~7 h. Both violate constraint 1 |
| **L12** | **ERA5 observing-system changes** (STATUS §5 risk 3) | No independent reference at 200 hPa | Not measured. `horizontal_divergence` is flat in both studies. L2's block test covers 2000–06 only | Trend, in principle | Cannot within this dataset |
| **L13** | **No NA peaks-over-threshold archive**; tails exist for the reference year only | Its shape depends on the unit of observation | — | EVT readiness, not the indicators | Chose not to *yet*: one pass, ~23 GB, once the unit is fixed |
| **L14** | **No backup** on the share (§16.3) | Not asked | Loss costs compute (~7 h NA + ~14 h global + ~2 h calibration), not information | — | One email to ITvO before anything is published |

Two limitations are new today:
- **L5**, the edge row.
- **L10's SOG trend disagreement.** It had never been computed, because S4
  had not been read.

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

**Total new cluster time: one job, 2 h 45 min hard cap, estimated 1–1.7 h.**
If step 2 has not finished by the cap, the budget is spent. The LOG/SOG
columns of §5 then stay as they are (baseline-derived), and nothing is
resubmitted without a new decision.

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

| | Prediction | If it holds | If it fails |
|---|---|---|---|
| **P1** | The new MOG table equals the job-1108148 audit MOG table to printed precision | B2's LOG/SOG tables are usable | Stop. Use no audit per-diagnostic number until the two runs are reconciled |
| **P2** | The 5 convention-invariant rows (`colson_panofsky`, `endlich`, `negative_richardson`, `vertical_wind_shear`, `wind_speed`) are identical in the audit and baseline tables at MOG and SOG | The restored baseline is the 09-07 table (D1) | Do not restore. The claim that these five are "bit-identical" needs rechecking before it goes into the appendix |
| **P3** | The audit `magnitude_pv` MOG level ratio rises from 0.61 to ~0.9–1.1, and the audit level gain concentrates in it | Rewrite "levels improved at every severity" as in §1 | Keep DECISION_RECORD's sentence |
| **P4** | The audit `magnitude_pv` SOG trend moves toward his +90 %, so his value falls inside our CI | `magnitude_pv` stays class B | Add "SOG trend below Prosser's" to its caveat |
| **P5** | **`ncsu1` at LOG:** level ratio within 0.75–1.33 and his +20.7 % inside our CI (§6 Q3) | Caveat: "agrees with Prosser at LOG; differs above LOG through the Ri floor on his finer stencil" | Caveat: "differs at every severity; unexplained" |
| **P6** | Finite census: 20 diagnostics finite everywhere in all 42 years. `f2d` at 1.000 (one-sided edges) or 0.9918 (NaN edges), nothing else | Record which one in L4 | Any other shortfall → locate the month (the census CSV is per year) and repair it via jobs/25 |
| **P7** | Sweep: 0 items. Thresholds: 21 increasing ladders, equal n = 3,039,966,720 for 20 diagnostics, and for `f2d` equal as well or lower by 24 × 721 × 1440 | Dataset certified | Repair list (§1, point 2) |
| **P8** | `ubf`'s reference-year tail is **not** NH-weighted (B1b) | L10 reads "structural" | L10 reads "regional (NA box)" |

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
- **Level ratios** are ours ÷ Prosser's read-off 1979 level.
- **"in"** says whether his rel lies inside our 95 % CI.
- **Values** are baseline. They are *final* for the five diagnostics marked ■
  (identical under both conventions, and verified again by P2).
- **"Pending"** marks cells the batch replaces.

| W&J # | Diagnostic | Class | Lvl MOG / SOG | Trend MOG / SOG vs his (in?) | Flag to hand to the econometrician |
|---|---|---|---|---|---|
| 6 | `vertical_wind_shear` ■ | **A** | 0.85 / 0.89 | +39 / +60 vs +38 / +55 (y/y) | — |
| 3 | `brown1` | **A** | 1.09 / 0.99 | +45 / +73 vs +45 / +71 (y/y) | — (pending: small F1+F12 shift) |
| 15 | `brown2` | **A** | 0.91 / 0.88 | +46 / +76 vs +46 / +69 (y/y) | Magnitude only: units lack an unknown constant L² (s⁻³). Never compare to published J kg⁻¹ s⁻¹. GPD scale σ is not physical, ξ is fine |
| 17 | `ti1` | **A** | 0.88 / 0.83 | +44 / +69 vs +42 / +59 (y/y) | — |
| 19 | `ti2` | **A** | 0.87 / 0.84 | +29 / +40 vs +27 / +32 (y/y) | — (pending) |
| 8 | `deformation` | **A** | 1.00 / 0.96 | +37 / +55 vs +38 / +54 (y/y) | Holds DEF, not DEF² (relevant for EVT) |
| 16 | `vorticity_squared` | **A** | 1.12 / 1.22 | +40 / +73 vs +37 / +69 (y/y) | — (pending; F12 moves 5–8 %) |
| 4 | `temperature_gradient` | **A** | 0.91 / 0.91 | +27 / +41 vs +26 / +40 (y/y) | — |
| 9 | `wind_speed` ■ | **A** | 0.94 / 1.01 | +74 / +101 vs +76 / +102 (y/y) | Sparse (0.26 % / 0.05 %), so the intervals are wide |
| 18 | `ngm1` | **A** | 0.93 / 0.88 | +52 / +79 vs +51 / +78 (y/y) | — |
| 14 | `nva` | **A** | 0.97 / 0.91 | +29 / +34 vs +28 / +36 (y/y) | Clipped at 0, so there is a mass of exact zeros (POT threshold must sit above it) (pending) |
| 12 | `rva_magnitude` | **A** | 0.90 / 0.92 | +29 / +33 vs +29 / +35 (y/y) | — (pending) |
| 5 | `horizontal_divergence` | **A** | 0.87 / 0.85 | +5 / +6 vs +7 / +8 (y/y) | **No trend, in either study.** A property of the field. Do not read it as a defect (pending; audit MOG +5 %, t = 1.01) |
| 1 | `magnitude_pv` | **B** | 0.61 / 0.61 (baseline, archived PV) | Audit MOG +42 %, t = 4.65 vs +50.6 (y). SOG pending | **The quantity is Sharman A18 truncated PV** computed on the 50-hPa stencil, carried in hPa units (×100 vs SI, not PVU). It is **not** ERA5 Ertel PV (64–80 % of its tail set differs). Never quote a PVU magnitude. The archived-PV version (baseline) shows the largest 2000–06 block anomaly (−17/−22 %) |
| 7 | `endlich` ■ | **B** | **0.51 / 0.46** | +38 / **+58** vs +28 / **+36** (y/**N**) | Level about half of Prosser's (stencil). **SOG trend steeper than his, beyond our CI.** ∂ψ/∂z uses the angle-secant reading, and the closed-form reading flips 41–66 % of the set (not measured on the series) |
| 10 | `ngm2` | **B** | **0.61 / 0.56** | +28 / +37 vs +34 / +42 (y/y) | Level ~40 % below Prosser's (∂T/∂z on the 50-hPa stencil; 2.8 % truncation error). Trends agree |
| 11 | `negative_richardson` ■ | **C** | 0.49 / 0.26 | +41 / +18 vs −8 / −14 (y/y), all n.s. | Too sparse to stand alone (0.019 % / 0.003 %: ~0.3 h per point per year at SOG). No trend in either study. Richardson family: the most stencil-sensitive (W&S 2022: thermodynamic diagnostics agree worst across datasets). Upper-tail values are negative numbers, which makes EVT awkward |
| 2 | `colson_panofsky` ■ | **C** | **0.14 / 0.12** | +62 / +68 vs +45 / +55 (y/y); ours n.s., his sig | Our NA exceedance is 7–8× sparser than his (share 0.10 / 0.16 vs 0.71 / 1.30), so this box sees a different subset of events. The trend agrees in size and fails significance on power. λ is constant (ICAO), so scale only |
| 20 | `f2d` | **C** | 0.43 / 0.34 (baseline variant C) | Audit MOG +25 %, t = 0.91 vs +0.3 (y). SOG pending | Sparse. Signed variant A, whose thresholds cross zero (read it from flip rates, never from relative threshold differences). ∂Q/∂t damped ×0.900 (diurnal). Month-edge steps one-sided or NaN (L4). The Williams (2017) Fig. 1 / Table 2 inconsistency about its sign is unresolved in the literature |
| 13 | `ubf` | **C** | **0.19 / 0.08** | +20 / +18 vs +32 / +44 (y/**N**); SOG ours n.s., his sig | Our NA box holds 0.35× / 0.10× an even share of the global tail; his holds 1.80× / 1.35×. It is single-level, so the stencil cannot explain this, and W&S 2022 find UBF the most dataset-robust diagnostic. **An implementation-level difference from the Prosser/Williams UBF that we cannot locate without their code.** Ours follows A30 with verified spherical geometry |
| 21 | `ncsu1` | **C** | **6.07 / 2.47** | +39 / +47 vs +4 / +10 (**N/N**); ours sig, his not | Correct A36 (ρ = 1.0000 cross-check; W&S Eq. 7). **It does not behave like Prosser's `ncsu1` at MOG/SOG.** Leading explanation: the Ri floor (×10⁵ where Ri ≤ 10⁻⁵) binds far more often on his 18-hPa stencil (§6 Q3). For EVT, the floor would plant an artificial mass in the extreme tail. Status at LOG: P5 |

**The econometric consequences.**
1. **Ensemble.** Use all 21, as Prosser does. Dropping the five class-C
   members moves the ensemble trend by +0.0 pp (MOG) and +0.3 pp (SOG).
   Leave-one-out spans 35.0–37.3 % (MOG) and 51.9–56.2 % (SOG), with
   `brown2` and `temperature_gradient` at the ends. No single member drives
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

Eight of them, with the exact flags in §5:
- Class C, as stand-alone indicators: `ncsu1`, `ubf`, `colson_panofsky`,
  `negative_richardson`, `f2d`.
- Class B, with the stated caveat: `magnitude_pv`, `endlich`, `ngm2`.
- Also, for magnitudes only: `brown2` (units) and `deformation` (DEF, not
  DEF²).

**Q3. Where does `ncsu1` stand, and is it usable?**

It is correctly computed. It does not reproduce Prosser's `ncsu1` above LOG.
It is usable as an ensemble member and not as a stand-alone indicator.

The explanation below fits all six published `ncsu1` numbers (level and
trend at LOG, MOG and SOG) and nothing on disk contradicts it:

| | LOG | MOG | SOG |
|---|---|---|---|
| Prosser: NA share of the global tail | **0.92** | **0.18** | **0.36** |
| Prosser: trend | **+20.7 %, p = 0.0009** | +4.3 %, p = 0.7 | +9.7 %, p = 0.5 |
| Ours (baseline): NA share | pending | 1.06 | 0.89 |
| Ours: trend | pending | +39 %, sig | +47 %, sig |

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

This is a hypothesis about data we do not have. Its one testable prediction
on our side is P5: at LOG, where his `ncsu1` is still dynamics-driven, ours
should agree with his on both level and trend. If P5 fails, `ncsu1` goes back
to "unexplained after four candidates". The class stays C either way.

**Q4. What is assumed rather than verified?**

1. The `audit_fixes` attribute and completeness of **all** 504 audit stores
   and all 12 `global_audit` stores. So far: count, size and one store. →
   **B1.**
2. That the audit NA series is finite everywhere. Exit codes exclude
   whole-variable placeholders, but not partial NaN. → **B2 census.**
3. `f2d` month-edge handling: one-sided or NaN. The record contradicts
   itself (C2). → **B2/B1.**
4. That horizontal derivatives on the 60°N edge row are finite and one-sided
   (L5). → **B1/B2.**
5. That the job-1108148 audit per-diagnostic table (produced after two
   script fixes on 09-15/16) is right. → **B2 recomputes it through a
   different code path (P1).**
6. That the `data_prosser/` baseline table was scored on the 09-07
   thresholds. Inferred from its mtime and the t = 1.97 match. → **P2.**
7. That `thresholds_2026-09-11.json` has equal sample sizes for all 21 and
   increasing ladders. The record states n once; the ladders are checked for
   Williams' table, not ours. → **B1.**
8. Every per-diagnostic 1979 level of Prosser. These are read off figures;
   validated in aggregate to 0.2–0.4 %, not per panel.
9. That 188/197/206 hPa are L137 levels 73–75. Prosser's text calls them
   "pressure levels" and says he used "a fixed Gaussian grid" at 0.25°,
   which is not fully consistent with the model-level reading. Only Prosser
   can settle it; the conclusion (a finer stencil than CDS offers) does not
   depend on it.
10. The size of the ERA5 vs ERA5.1 difference at 200 hPa. Bounded (L2), not
    measured.
11. That `essentials/` reproduces the ADA stores bit for bit. Shown locally
    on test files; the ADA check (`essentials/jobs/0`) has never run. This
    matters only if `essentials/` becomes the published code.
12. That rojak's gradients on the global grid treat longitude as periodic,
    which affects 2 of 1440 columns of the calibration year. Not checked;
    the effect is below the percentile noise.

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
2. **The bias it causes is in the level and stable in time.**
   - Level ratios drift by ≤ 0.006 over 41 years.
   - The deficit deepens with severity in DJF/MAM/SON and not in JJA, which
     is where a stencil that misses sharp jet-level shear should matter.
   - The stencil-group median level ratio is 0.73, against 0.93 for the
     single-level diagnostics.
3. **Our trends match his to within 4 %** at every annual severity, and 25/25
   Table-1 cells contain his value. That bounds the combined effect of every
   pipeline difference, the stencil included. Prosser's own contrast with Lee
   et al. was about *where* hotspots appear (and Lee used a 100-hPa stencil);
   within his North Atlantic box our changes are his.
4. **Where the stencil changes behaviour, the diagnostic is flagged and the
   result does not depend on it.** That applies to the Richardson family and
   to `ncsu1` through its floor. Dropping all five flagged diagnostics moves
   the ensemble trend by ≤ 0.3 pp.

What we cannot say is how large the stencil effect on the trend is *on its
own*. Say that in the paper, in those words (C11).

---

## 7. After the batch

1. Paste the three consistency lines (P1, P2, P6/P7) and the LOG/SOG rows of
   the 16 changed diagnostics into §5. Apply the class rule in §4.3.
2. Make the corrections C1, C4, C5, C9 and C11 in `DECISION_RECORD.md`,
   `FORMULAS_AND_DECISIONS.md` and `essentials/DECISIONS.md`.
3. Then stop. Next come the unit of observation, then the NA
   peaks-over-threshold archive (L13), then the ITvO email (L14).

**Sources.** Prosser et al. (2023) pp. 2, 5, Fig. 4, SI Figs S4-LOG/SOG
(`Articles/`). Williams & Storer (2022) pp. 1428–1432, Tables 1–2. Sharman
et al. (2006) App. A. ECMWF on ERA5.1:
[NCAR GDEX summary](https://gdex.ucar.edu/news/era51-corrections-to-era5-stratospheric-temperature-2000-2006/),
[Simmons et al. 2020, ECMWF TM 859](https://www.ecmwf.int/en/elibrary/81149-global-stratospheric-temperature-bias-and-other-stratospheric-aspects-era5-and).
Everything else is from the project docs and the repository as read today.
