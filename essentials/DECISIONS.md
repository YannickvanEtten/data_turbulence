# DECISIONS — what this code does, and why

Every choice in `essentials/` is listed here once, with its reason, its source,
and where it lives in the code. Longer histories are in the original repo's
`STATUS.md` (section numbers given as §); the appendix-oriented decision
register with page-level citations is `FORMULAS_AND_DECISIONS.md` (D1–D11) in
the Claude project. If this file and a docstring disagree, fix one of them —
do not leave both.

Diagnostic numbers (#1–#21) are always the **Williams & Joshi (2013) Table 1**
order, as in `config.DIAGNOSTICS`. Some older notes use the Williams (2017)
Table 2 order; use names when in doubt.

Contents

1. Data
2. The 21 diagnostics
3. The two conventions (baseline vs audit) — the open decision
4. Numerical choices shared by both conventions
5. Severity thresholds
6. Exceedance, frequencies and trends
7. Known limitations of the dataset
8. How correctness was established
9. What was deliberately left out of `essentials/`

---

## 1. Data

| Choice | Value | Why | Code |
|---|---|---|---|
| Source | ERA5 pressure-level product on CDS | the standard product; Lee et al. 2023 and Wong et al. 2025 use it | `config.DATASET` |
| Variables | u, v, T, z, and archived d, vo, pv | u, v, T, z are all Prosser needs; d, vo, pv are only used by the *baseline* convention (§3) | `config.VARIABLES` |
| Levels | 175, 200, 225 hPa; diagnostics evaluated at **200 hPa** | 175/200/225 is the finest stencil CDS offers around 200 hPa. Prosser's 188/197/206 hPa are model levels 73–75, MARS-only (§11.6). Other levels would need a full re-run | `config.PRESSURE_LEVELS_HPA`, `TARGET_LEVEL_HPA` |
| Time | 3-hourly, 1979–2020 | Prosser's period and sampling | `config.TIMES_3H`, `TREND_YEARS` |
| Grid | 0.25° | native | `config.RESOLUTION_DEG` |
| Trend region download | 30–60N, 75W–0 | a superset of Prosser's box (36–60N, 55–10W) | `config.DOMAINS` |
| Calibration region | global, year 2000, all days | Prosser calibrates on a global year; the full contiguous year replaced a 48-day sub-sample, which was systematically biased in the tail (§15.4) | `config.DOMAINS`, `CALIBRATION_YEAR` |
| Global request size | three ~10-day blocks per month, joined by concatenation | a full global month (5,208 fields) in one CDS request was never tested; GRIB messages are self-contained so `cat` gives a valid file (§15.7) | `download.py` |
| Integrity check | variables, levels, time steps **and grid size** | a global and a North Atlantic month have identical counts; only the grid tells them apart (§9a) | `download.check_file` |
| Concurrency | `%4` for downloads | CDS refused 8 of 12 requests at 8 concurrent (§11.9) | `jobs/1a`, `1b` |
| ERA5 vs ERA5.1 | ERA5 throughout | Prosser used ERA5.1 for 2000–2006; ERA5.1 is MARS-only (§11.7) | — |

## 2. The 21 diagnostics

Order and numbers as in Williams & Joshi (2013) Table 1; equation numbers are
Sharman et al. (2006) Appendix A. "conv." marks the diagnostics whose values
depend on the convention of §3 (F1 = `dy`, F12 = `fields`, F6 = `f2d_variant`).

| # | key | formula | eq. | conv. |
|---|---|---|---|---|
| 1 | `magnitude_pv` | \|PV\|; archived Ertel PV, or −g(ζ+f)∂θ/∂p | A18 | F12 |
| 2 | `colson_panofsky` | λ²(S_v² − N²/Ri_c), Ri_c = 0.5, λ = ICAO z(200) − z(225) | A4 | — |
| 3 | `brown1` | (0.3 ζ_a² + D_sh² + D_st²)^½, ζ_a = ζ + f | A13 | F1, F12 |
| 4 | `temperature_gradient` | \|∇_H T\| | A23 | F1 |
| 5 | `horizontal_divergence` | \|δ\|, δ = ∂u/∂x + ∂v/∂y | A33 | F12 |
| 6 | `vertical_wind_shear` | S_v = \|∂**V**/∂z\| | A3 | — |
| 7 | `endlich` | \|**V**\| · \|∂ψ/∂z\|, ψ = wind direction | A25 | (F10) |
| 8 | `deformation` | DEF = (D_sh² + D_st²)^½ | A17 | F1 |
| 9 | `wind_speed` | \|**V**\| | A24 | — |
| 10 | `ngm2` | \|∂T/∂z\| · DEF | A29 | F1 |
| 11 | `negative_richardson` | −Ri = −N²/S_v² | A1–A3 | — |
| 12 | `rva_magnitude` | \|u ∂ζ/∂x + v ∂ζ/∂y\| | W&J | F1, F12 |
| 13 | `ubf` | \|∇²Φ − 2J(u,v) − fζ + βu\| | A30 | F1, F12 |
| 14 | `nva` | max{−u ∂(ζ+f)/∂x − v ∂(ζ+f)/∂y, 0} | A37 | F1, F12 |
| 15 | `brown2` | Φ S_v² / 24 | A14 | F1, F12 |
| 16 | `vorticity_squared` | ζ² | A21 | F12 |
| 17 | `ti1` | S_v · DEF | A15 | F1 |
| 18 | `ngm1` | \|**V**\| · DEF | A28 | F1 |
| 19 | `ti2` | S_v · (DEF − δ) | A16 | F1, F12 |
| 20 | `f2d` | ½ D/Dt[(∂u/∂θ)² + (∂v/∂θ)²], signed or absolute | A9 | F1, F6 |
| 21 | `ncsu1` | max(u∂u/∂x + v∂v/∂y, 0) · \|∇ζ\| / max(Ri, 10⁻⁵) | A36 | F1, F12 |

Formula checks against the papers (all closed): all 21 transcribed and cited
(`AUDIT_diagnostics_vs_literature.md`); the 0.3 in `brown1` and its use of
*absolute* vorticity (A13); both clips in `nva`/`ncsu1` (A36, A37); the sign of
`ubf` (A30 under an absolute value); **all 21 are one-tailed upper criteria**
(Williams 2017 Table 2: every severity ladder increases, including the two
negative ones) — §18.16.

Seven of the 21 were hand-written in the original because rojak was wrong or
lacked them (Richardson N²/S_v instead of N²/S_v², hence also CP and NCSU1;
UBF sign, Laplacian and β; frontogenesis cross term; RVA missing; Brown2).
In `essentials/` all 21 are written out, and the 14 formerly taken from rojak
reproduce rojak @ 25b8685 exactly (§8).

## 3. The two conventions — the open decision

Three implementation choices differ between the two series on the share. No
CAT paper states its choice on any of them (§18.4).

**Where the decision stands (2026-09-15).** `FORMULAS_AND_DECISIONS.md`
records `audit` as the series to report and `baseline` as the documented
sensitivity. That is the working choice, not yet a final one: `essentials/`
therefore implements both, with `audit` as the default
(`config.DEFAULT_CONVENTION`), so either can be produced and compared at any
time. They are switches in
`config.CONVENTIONS`; every output store records the convention it was built
under, and thresholds and series of different conventions refuse to mix.

| switch | baseline (`derived/north_atlantic`, `thresholds_2026-09-07`) | audit (`derived/north_atlantic_audit`, `thresholds_2026-09-11`) |
|---|---|---|
| `dy` (F1) | `geodesic`: true meridian arc, then × a/M | `map`: a·Δφ, then × a/M |
| `fields` (F12) | `archived`: ERA5 vo, d, pv | `computed`: from u, v, T, z |
| `f2d_variant` (F6) | `C`: \|½ D/Dt Q\| | `A`: +½ D/Dt Q |
| `endlich` (F10) | `angle` | `angle` |

### What each switch is, and the case for each side

**F1 — meridional grid spacing.** MetPy (and rojak) measure dy as the true
meridian arc M·Δφ and then multiply by the meridional scale factor a/M, which
applies the ellipsoid correction twice: an error of +0.42 % at 30N to −0.27 % at
75N in every ∂/∂y. With dy = a·Δφ the error against exact ellipsoidal
gradients falls from 4.2×10⁻³ to 2×10⁻⁷ (§18.5, `tests/test_operators.py`).
*This is an accuracy fix rather than a trade-off*; the only argument for
`geodesic` is "same as MetPy". Reach: 13 of 21 change, 3 materially
(`brown1`, `deformation`, `ngm2`; flip rate 0.18–0.76 %).

**F12 — where ζ, δ and PV come from.**
- *For `computed`:* Prosser (2023, p. 2) downloads only u, v, T, z and computes
  all 21 from them; this matches his method.
- *For `archived`:* ERA5's ζ and δ are the IFS model's own prognostic spectral
  variables (u and v are derived from them), so recomputing them on a 0.25° grid
  is a lossy round trip; archived PV is full Ertel PV on the model's fine
  vertical grid, whereas A18 keeps only the vertical term on a 50 hPa stencil.
- Reach: 10 of 21 change, all materially (North Atlantic flip rate at p97):
  `magnitude_pv` 64 %, `rva_magnitude` 11.5 %, `nva` 10.3 %,
  `horizontal_divergence` 7.5 %, `vorticity_squared` 5.2 %, `ncsu1` 4.9 %,
  `ti2` 2.5 %, `brown1` 1.6 %, `ubf` 1.1 %, `brown2` 0.6 %. The change is
  latitude-structured: computed fields lose tail frequency in the tropics and
  gain it in the extratropics (predicted before the run, §18.12). All
  calibration thresholds of this family fall by 3–8 % (§18.14).
- `magnitude_pv` is really a change of quantity (full Ertel PV vs A18), not
  just of provenance; report it separately.

**F6 — the reading of Sharman A9 for `f2d`.**
- *For `A` (signed):* Williams & Storer (2022) Eq. (3), p. 1427 write
  F = D/Dt\|∂**u**/∂θ\|² with no sign change or absolute value — the only place
  in Prosser's lineage where it is written as an equation. Prosser's Figure 4
  panel (d) shows a flat, insignificant trend (+0.3 %), which `C` does not
  reproduce (+26 %).
- *For `C` (magnitude):* Williams (2017) Fig. 1 plots this diagnostic on a
  one-signed 0–300 axis, and Williams (2017) p97 / W&J (2013) median = 13.6,
  which a signed, zero-centred quantity cannot give. That is an inconsistency
  between two Williams papers, worth a sentence in any write-up.
- Reach: 1 diagnostic, but ρ(A, C) = 0.03 and 68–70 % of the tail flips.

### What the choice changes in the results (measured, §18.15)

- **Trends: essentially nothing.** Both series: 25/25 season × severity trends
  significant, 25/25 contain Prosser's value. Annual trend ratio
  (ours ÷ Prosser), LOG → SOG: baseline 0.96 / 0.96 / 0.98 / 0.98 / 0.99,
  audit 0.97 / 0.98 / 1.01 / 1.02 / 1.04.
- **Levels: audit is closer to Prosser.** Annual 1979 level ratio: baseline
  0.881 / 0.857 / 0.841 / 0.820 / 0.797, audit 0.904 / 0.888 / 0.880 / 0.865 /
  0.842 — about a fifth to a quarter of the deficit closes.
- Combined reach: 5 diagnostics identical in both series (`colson_panofsky`,
  `endlich`, `negative_richardson`, `vertical_wind_shear`, `wind_speed`),
  16 numerically different, 13 materially different (≥ 0.5 % of the
  exceedance set at some severity).

### How to decide, and how to try other combinations

The three switches are independent. A mixed convention (e.g. F1 on with archived
fields) is one new line in `config.CONVENTIONS`, but it needs its own global
year (Stage 2b, ~14 h), thresholds (Stage 3, ~2 h) and North Atlantic series
(Stage 2a, ~7 h, ~300 GB). Thresholds from one convention must never be applied
to a series of another: for `f2d` (A vs C) the result would be meaningless.
Test any combination first on one month with `main.py demo`.

**F10 (`endlich`) is a fourth, still-open sensitivity.** "Centred second-order
differences" (W&S 2022, p. 1428) can mean differencing the angle (rojak, used
in both series) or the components (closed form). They differ by 41–66 % of the
exceedance set, and `endlich`'s level ratio against Prosser is 0.51 although
W&S rank it among the most resolution-robust diagnostics (§18.16.6).

## 4. Numerical choices shared by both conventions

| Choice | What | Why / source | Code |
|---|---|---|---|
| Horizontal differences | `np.gradient` on cumulative grid distance: 2nd-order centred inside, 1st-order one-sided at domain edges | MetPy/rojak convention | `grid_difference` |
| Scale factors | k_x = a/(N cos φ), k_y = a/M from PROJ `+proj=latlon` | MetPy/rojak; PROJ agrees with the closed form to 10⁻¹⁰; at the poles PROJ gives a finite k_x ≈ 10⁵ | `grid_geometry` |
| Derivatives of u, v | with curvature terms (`wind_derivatives`) — never the scalar operator | a zonal flow has vorticity u tanφ/R without shear; without the terms ζ is 100 % wrong and ∇²Φ 99 % wrong (Stokes/divergence theorems, §4g) | `wind_derivatives`, `ubf` |
| ∇²Φ in UBF | divergence of ∇Φ with curvature terms | same | `ubf` |
| β in UBF | 2Ω cos φ / R, R = 6 371 008.77 m | rojak used f/R | `ubf` |
| Vertical derivative | ∂f/∂z = g (∂f/∂p)/(∂Φ/∂p), centred over 175–225 hPa | real geopotential thickness, not a standard atmosphere | `d_dz` |
| ∂/∂θ | (∂f/∂p)/(∂θ/∂p) on the pressure levels, no isentropic remap | same pattern | `f2d` |
| Wind direction difference | wrapped to the shorter arc; needs evenly spaced levels (asserted) | rojak `angles_gradient` | `d_angle_dp`, `prepare` |
| θ | T (p/1000)^−κ, κ = 2/7 | rojak | `potential_temperature` |
| DEF | un-squared (rojak's DEF diagnostic is DEF²) | A17 and Williams (2017) tabulate DEF; matters for extreme-value fits (the tail index doubles under squaring) | `deformation` |
| Ri floor in NCSU1 | max(Ri, 10⁻⁵) from below | A36 | `ncsu1` |
| CP length scale λ | ICAO height difference 200–225 hPa: a constant | a state-dependent λ would make the unit and stencil choices non-inert | `colson_panofsky` |
| Units left as computed | `brown2` in s⁻³ (published m² s⁻³ implies an unstated L²); A18 PV 100× SI (p in hPa); CP in m² s⁻² | constant factors cancel exactly in a percentile threshold | docstrings |
| Time derivative in f2d | centred 3-hourly, one-sided at the first and last step of **each file** (months are not stitched together) | as in the original production run; 3-hourly sampling damps a diurnal signal by 36 % (§4d) | `f2d`, `compute_month` |
| Chunking | 1 day at a time with a 1-step overlap, trimmed | memory only; result identical to a single pass (tested) | `compute_month` |
| Storage | float32 zarr, 200 hPa only, one store per month | percentiles need nothing finer; other levels require a re-run | `run_month` |

## 5. Severity thresholds

| Choice | Value | Why | Code |
|---|---|---|---|
| Percentiles | LOG 97.0, LMOG 99.1, MOG 99.6, MSOG 99.8, SOG 99.9 | Williams (2017) Table 1, used by Prosser | `config.SEVERITIES` |
| Sample | all cells × all 2,928 steps of global 2000 | Prosser recomputes thresholds from ERA5 (he publishes none) | `thresholds.calibrate` |
| Weight | cos φ per cell | a regular lat-lon grid over-samples the poles | `latitude_weights` |
| Quantile rule | weighted Hazen: P_i = (C_i − w_i/2)/W, linear interpolation, ties merged | reduces exactly to (i − ½)/n with equal weights | `weighted_percentile` |
| Method on 3×10⁹ values | keep values ≥ cut (≈ p88) plus W and W_below; identical result | a full sort needs > 100 GB (§15.2) | `scan_tail`, `thresholds_from_tail` |
| Cut sample | every 7th step | a multiple of 8 would sample a single time of day (§15.3) | `config.TAIL_CUT_STRIDE` |
| Guard | stop if the kept tail does not start below p97 − 0.005 | otherwise np.interp silently returns the cut as the threshold | `thresholds_from_tail` |
| File format | JSON, schema 1, same as the original + `convention` in provenance | both code bases read each other's files | `save_thresholds` |

## 6. Exceedance, frequencies and trends

| Choice | Value | Why | Code |
|---|---|---|---|
| Exceedance | D ≥ threshold, NaN stays undefined | a NaN must not count as "no turbulence" | `year_frequencies` |
| Ensemble | exceedance first, then mean over the diagnostics defined at that cell and time | Prosser's order; averaging raw values first is a different method | `year_frequencies` |
| Box | 36–60N, 55–10W | Prosser Fig. 2 caption (§11.4) | `config.PROSSER_BOX` |
| Weight | cos φ | area | `year_frequencies` |
| Seasons | DJF(Y) = Jan + Feb + Dec of year Y | the record starts 1 Jan 1979 | `config.SEASON_MONTHS` |
| Hours | frequency × 24 × days in season | Prosser Table 1's unit | `print_summary` |
| Trend | OLS on 42 annual values; change and 95 % interval from the fitted 1979 and 2020 values | a comparison of two raw years is not a test (§12.8) | `ols`, `trend_row` |
| Verdicts | "significant" (\|t\| ≥ t₀.₉₇₅,₄₀ = 2.021) and "Prosser inside our interval" — no other automatic verdicts | automatic judgements flipped on rounding (§15.6) | `trends.run` |

## 7. Known limitations of the dataset

- **Level deficit of 10–16 % (audit) against Prosser**, deepening with severity,
  in winter and spring but not summer. Attributed to the 50 hPa vertical stencil
  (Prosser: 18 hPa); four independent lines of evidence (§18.15). It does not
  affect the relative trends (≤ 0.05 in the ratio).
- `ncsu1`: level ratio 6.07, unexplained after three candidate causes (§18.3).
- `horizontal_divergence`: no trend in any season, as in Prosser; open
  pre-registered question whether it gains one under computed fields (§18.15.5).
- `f2d`: damped by 36 % for diurnal signals at 3-hourly sampling (§4d).
- ERA5 rather than ERA5.1 for 2000–2006 (§11.7); ERA5 observing-system changes
  across 1979–2020 are an unexamined source of time-varying bias (§5 risk 3).

## 8. How correctness was established

| Check | Result | Where |
|---|---|---|
| Diagnostics vs the original code (rojak + `2_diagnostics.py`), real ERA5 day (NA, 8 steps) | **all 21 bit-identical**, baseline, audit and baseline+F10 | `tests/compare_with_original.py` |
| Same, 3-day NA file (chunked, time derivative across chunks) | bit-identical, all three | same |
| Same, `era5_validation_subset.nc` (NetCDF, 61×61, 2 steps) | bit-identical, all three | same |
| Same, synthetic global grid (1°, poles, 0–360 longitudes) | bit-identical, baseline and audit | same |
| Chunk size 1, 2 days vs whole file | identical | `compute_month` |
| Thresholds vs original `ada/tail_thresholds.py` | identical JSON | `tests/compare_stages_3_4_with_original.py` |
| Tail method vs full sort | ≤ 2×10⁻⁸ relative | same, and `tests/test_operators.py` |
| Frequencies vs original `aggregate.py` + `full_trend_check.py` (3 synthetic years, all seasons, severities, 21 + ensemble) | ≤ 7×10⁻¹⁶ relative (summation order only) | `tests/compare_stages_3_4_with_original.py` |
| OLS vs original | identical | same |
| Geometry against exact ellipsoid, scale factors vs closed form, wrap, Hazen, guard | pass | `tests/test_operators.py` |
| **On ADA, against the stores on the share** | to run: `jobs/0_check_against_original.sbatch`, `jobs/0b_check_thresholds_and_trends.sbatch` | README |

Literature checks of the formulas themselves are in the original repo
(`AUDIT_diagnostics_vs_literature.md`, `FORMULA_AUDIT.md`, §18.16).

## 9. What was deliberately left out of `essentials/`

- rojak as a dependency (all formulas are written out; equality is tested).
- The single-diagnostic probes F2, F3, F11 (subsets of F12) and the ab-compare
  harness; the Miller-form frontogenesis; the W&J / Sharman / Lee magnitude
  comparison tables (constant factors make them smell tests only).
- The 48-day sub-sampled calibration and the whole-array calibration path.
- The `sign="-"` / `"either"` machinery (all 21 are `"+"`).
- Plotting and the Prosser figure reproductions (`ada/prosser_*.py`).
- Month-boundary stitching (`chunk_stitch.py`) — the production run never used it.
