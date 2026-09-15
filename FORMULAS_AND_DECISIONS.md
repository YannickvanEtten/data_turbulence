# FORMULAS AND DECISIONS — appendix source material

**Purpose.** Everything needed to write the methods appendix, with each claim
traced to a printed equation and page. This is *source material*, not the
appendix: it records what was implemented, what it is in line with, what was
chosen where the literature is silent, and what each choice was measured to be
worth. The prose comes later.

**Status: 2026-09-15.** Companion documents: `STATUS.md` §18 (the audit and the
re-derive), `AUDIT_diagnostics_vs_literature.md` (the 21 individual audits),
`RESULT_full_trend_42yr.md` (the baseline trend result).

**The reported dataset is `derived/north_atlantic_audit/`**, calibrated on
`calibration/thresholds_2026-09-11.json`. The baseline series
`derived/north_atlantic/` on `thresholds_2026-09-07.json` is retained as the
documented sensitivity. §D1 is the decision; §4 is the evidence.

---

## 0. READ THIS FIRST — the diagnostic numbering is ambiguous

**Two numbering schemes are in circulation and they disagree.** Both appear in
this project's own code and notes.

| | Williams & Joshi (2013) Table 1 order | Williams (2017) Table 2 order |
|---|---|---|
| `magnitude_pv` | **1** | 10 |
| `colson_panofsky` | **2** | 3 |
| `ubf` | **13** | 17 |
| `endlich` | **7** | 14 |
| `ncsu1` | **21** | 19 |
| `f2d` | **20** | 4 |

`2_diagnostics.py`'s `REFERENCE_TABLE["num"]` and its section headers use the
**W&J (2013)** scheme throughout. Two places deviate — the
`endlich_component_shear` header says "#14 (Williams 2017 order)" and the
`magnitude_pv_a18` header says "#1 (W&J) / #10 (Williams 2017)". Several job
scripts and audit notes written in September use the Williams-2017 scheme
("#19 `ncsu1`", "#14 `endlich`").

**For the appendix: use names, not numbers, or state the scheme once and use
W&J throughout** — it is what the code carries. Do not copy a "#N" from a note
without checking which scheme it came from.

---

## 1. THE DECISION REGISTER

Each entry: what was chosen, what it is in line with, and what the alternative
was measured to cost. "Flip rate" is the cos φ-weighted symmetric difference of
the exceedance set, W(A xor B)/W(A or B) — the only metric that can say whether
a change alters which cells are flagged, since thresholds are percentiles of
our own data and any constant factor cancels exactly.

### D1 — Field provenance: FOUR variables, not seven

**Chosen.** ζ, δ and PV are computed from u, v, T, z. ERA5's archived `vo`,
`d` and `pv` are not read. (`FixSet.four_variable_provenance`, F12.)

**In line with.** Prosser et al. (2023) p. 2 downloads exactly four fields —
"zonal and meridional wind speed, dry bulb temperature, and geopotential
height" — and states "The 21 turbulence diagnostics were then calculated from
the extracted reanalysis fields."

**The honest framing, which the appendix must carry.** This is *not* a
correction. ERA5's ζ and δ are the IFS's own prognostic spectral variables —
IFS Documentation CY49R1, Part III §2.2.7: *"The model state at time t is
defined by the spectral coefficients of ζ, D, T, q and ln ps"* — and u, v are
derived **from them** by inverse spectral transform. Recomputing ζ from
archived u, v on a 0.25° grid is a lossy round trip, not a refinement. The
archived field is closer to the model; the computed field is closer to Prosser.
Likewise ERA5's archived `pv` is diagnosed on the model's L137 vertical grid,
so its ∂θ/∂p is better resolved than anything obtainable from three pressure
levels. **Substituting Sharman A18 is worse physics and closer to the
replication target.**

**Measured cost.** Touches 10 of 21; **zero inert** — every diagnostic that
reads an archived field crosses the 0.5 % bar, every one that does not is
bit-identical. NA flip rates p97→p99.9: `magnitude_pv` 64→80 %,
`rva_magnitude` 11.5→17.4 %, `nva` 10.3→15.7 %, `horizontal_divergence`
7.5→12.5 %, `vorticity_squared` 5.2→8.3 %, `ncsu1` 4.9→7.5 %, `ti2` 2.5→5.0 %,
`brown1` 1.6→2.8 %, `ubf` 1.1→1.0 %, `brown2` 0.6→1.1 %.

**Mechanism, confirmed three independent ways.** A finite-difference ζ on a
0.25° grid is *damped* relative to ERA5's spectral ζ, and the excess
small-scale power is **tropical** (convection), not midlatitude. Confirmed by
(i) the flip rates, (ii) the latitude tilt — tropics negative, extratropics
positive, in 6 of 7 of the vorticity family, and (iii) **every threshold
moved down**, by 3–8 % for that family, which is what a thinner tail does to a
percentile.

### D2 — Horizontal metric: plate-carrée map spacing, one scale factor

**Chosen.** dy is the map distance a·dφ, so the meridional map-scale factor is
applied exactly once. (`meridional_metric_patch`, F1.)

**In line with.** MetPy's spherical convention for the *operator* —
curvature-corrected vector-component derivatives (`vector_derivative`), which
rojak inherits — and a **deliberate departure from MetPy for the spacing**.
MetPy's `nominal_lat_lon_grid_deltas` builds `dx = geod.a * Δλ` (nominal
equatorial map distance) but `dy` from `geod.inv` along the meridian (true
geodesic arc), then `geospatial_gradient` multiplies both by scale factors from
`Proj(CRS('+proj=latlon')).get_factors()` — `a/(N cos φ)` and `a/M`. For x that
is exact; for y, `M·dφ` is already ∂/∂y and `a/M` corrects it a second time.

**Evidence.** `tests/test_metric_terms.py` compares against two analytic fields
with exactly-known ellipsoidal gradients: MetPy/rojak's convention is off by
**4.2 × 10⁻³**, ours by **2 × 10⁻⁷** — four orders better, and the residual is
finite-difference truncation rather than a metric error. `a/M` runs 1.0068 at
the equator to 0.9966 near the pole, so the defect is a U symmetric in |φ|.

**Measured cost.** Touches 13 of 21, material in 3 (`brown1`, `deformation`,
`ngm2`); flip 0.18–0.76 %, ρ ≥ 0.99974. **Widest reach, smallest effect.**

**Appendix sentence.** *Horizontal gradients follow the MetPy spherical
convention (curvature-corrected vector-component derivatives), with the
meridional grid spacing taken as the plate-carrée map distance a·dφ so that the
meridional map-scale factor is applied exactly once; verified against exact
ellipsoidal gradients to 2 × 10⁻⁷.* **Worth reporting upstream to MetPy and
rojak** — the test is written.

### D3 — Frontogenesis: the signed material derivative (variant A)

**Chosen.** `F2D_DEFAULT_VARIANT = "A"`, the signed isentropic form, not its
magnitude.

**In line with.** Williams & Storer (2022) Eq. (3), p. 1427, which writes the
signed material derivative. Sharman (2006) A9 is the isentropic form; Sharman's
own Table B1 units (m² s⁻³ K⁻²) require the **squared** inner quantity, which
contradicts the printed exponent ½ in A9 — the table wins. Williams (2017)
Table 2 independently tabulates the frontogenesis function in
**10⁻⁹ m² s⁻³ K⁻²**, confirming the squared form.

**Measured cost.** One diagnostic, and the largest single change in the audit:
**ρ = 0.0275**, flip 68–70 %. Variant C (the magnitude) and variant A are
barely the same ranking.

**Counter-evidence retained, not suppressed.** The tests record what they used
to assert and why; `tests/test_analytic.py::test_variant_C_still_reproduces_the_archived_run`
keeps variant C reproducible.

**Known and irreducible.** A centred difference of a diurnal signal at 3-hourly
sampling returns the true derivative times sin(Wh)/(Wh) = **0.637**. F2D is
systematically damped by ~36 % for anything sub-diurnal. This is a property of
the sampling, not the code, and no change repairs it. **State it.**

### D4 — Deformation: DEF, not DEF²

**Chosen.** The square root is applied inside `compute_rojak_diagnostics`
(F4), so every store holds DEF.

**In line with.** Sharman (2006) A17: DEF = (D_SH² + D_ST²)^½. Williams (2017)
Table 2 tabulates flow deformation in **10⁻⁶ s⁻¹** (not s⁻²), independently
confirming the un-squared form.

**Why it was invisible and still matters.** Squaring is monotone on a
non-negative field, so the exceedance field is identical and no replication
number changes. It stops being harmless the moment magnitudes are used: a GPD
fitted to DEF² is not a GPD fitted to DEF (the tail index doubles). Relevant to
the econometrics, not to the replication.

### D5 — Grid assumptions asserted, not assumed

**Chosen.** `_assert_grid_assumptions` raises on non-uniform pressure levels or
a domain lying entirely south of ~3.14°N (F8).

**Why.** Two silent failures. (i) rojak's Endlich/DirectionalShear wrap the
wind-direction difference through `_WrapAroundAngleArray.__sub__`, which NumPy
only calls when the coordinate reduces to a scalar — i.e. only on uniform
spacing. On a non-uniform vertical grid the wrap-around stops working with no
error. (ii) rojak's `coriolis_parameter` decides degrees-vs-radians with
`latitude.max() > π`, using `max()` not `abs().max()`, so any domain entirely
south of ~3.14°N is read as radians and f is wrong by ~57×. 175/200/225 and
30–60°N are both safe; the assertions keep them so.

### D6 — Vertical stencil: 175/200/225 hPa at 200 hPa — a CONSTRAINT

**Chosen.** Not a choice. 188/197/206 hPa are ECMWF L137 model levels 73/74/75,
available only from the MARS tape archive; the CDS pressure-level product
offers 150/175/200/225/250 around 200 hPa and nothing between. **175/200/225 is
the finest stencil that exists at 200 hPa on CDS.** Model levels are also
hybrid sigma-pressure, i.e. terrain-following, so even with MARS access the
comparison would not be exact.

**Others use the same stencil.** Williams & Storer (2022) use ERA5 pressure
levels; Wong et al. (2025) — rojak's own author — uses 175/200/225 evaluated at
200 hPa. **Prosser is the outlier, not this project.**

**Measured cost, four independent lines:**
1. The deficit deepens with severity — the signature of a wide stencil
   under-resolving sharp shear maxima (§11.3).
2. The stencil-group split, pre-registered and passed: median 1979 level ratio
   **0.73** for the ten diagnostics that use the 175/225 stencil against
   **0.93** for the eleven that do not (§17.3).
3. Level/trend separation: recalibrating moved levels up to 4.6 % and the
   trend by ≤ 1 percentage point in all 25 cells (§15.11).
4. **Seasonal fingerprint** (new, §18.15): the severity-deepening deficit is a
   winter and spring phenomenon and is **absent in summer**. LOG−SOG level
   ratio spread: DJF +0.097, MAM +0.101, SON +0.077, **JJA −0.018**. Sharp
   200 hPa shear layers are a jet-stream phenomenon; the summer jet is weaker
   and more diffuse, so there is less structure for a wide stencil to miss.
   JJA instead shows a roughly uniform ~15 % deficit with no tail structure.
   **These are two different mechanisms and should not be averaged.**

**Declined option, costed:** ERA5 model levels 73–75 via MARS ("ERA5
complete") is the same data volume already held. The binding constraint is
tape-retrieval throughput — realistically weeks of queue for 42 years
3-hourly — and it would not change the trends. Storage is not the obstacle
(0.96 TB free).

**Appendix sentence.** *Diagnostics are evaluated at 200 hPa on a
175/200/225 hPa stencil, the finest available in the ERA5 pressure-level
product; Prosser et al. (2023) use model levels 73–75 (≈188/197/206 hPa) from
the MARS archive. The wider stencil biases the exceedance level by 10–16 %,
deepening with severity in winter and spring and flat in summer, and changes
the fitted 1979–2020 relative trend by at most one percentage point in any
season or severity.*

### D7 — Severity thresholds: cos φ-weighted percentiles of a full global year

**Chosen.** Williams (2017) Table 1 percentile ladder, verbatim:
**97.0 / 99.1 / 99.6 / 99.8 / 99.9**, applied as cos φ-weighted global
percentiles of the year-2000 reference year, computed on the full contiguous
year (not a sub-sample).

**In line with.** Williams (2017) Table 1, p. 580, derived from a log-normal
EDR distribution constrained by an observed 3.0 % light and 0.4 % moderate
probability. Prosser recomputes thresholds from ERA5 rather than reusing
Williams' numbers — confirmed from his own text — so latitude weighting is
required by the method, not an embellishment.

**Prosser publishes no threshold values.** Confirmed by reading his Supporting
Information: Figures S1–S5, no tables. The only published per-diagnostic
threshold table is Williams (2017) Table 2, and it is GFDL-CM2.1, so magnitudes
are not comparable — Williams says so himself.

**Measured cost of the sub-sample alternative.** 48 days → full year moves
thresholds a median 0.39 % and the fitted trend by **at most one percentage
point** in any of 25 cells, changing no verdict (§15.11). But the sub-sample was
systematically different, not merely noisier, and the difference grows with
severity; the two diagnostics moving opposite were both Richardson-based, which
localises it to **N², the static stability**. Report levels from the full-year
calibration for that reason independent of precision.

### D8 — Comparison box: 36–60°N, 55–10°W

**In line with.** Prosser (2023) Figure 2 caption. A strict subset of the
30–60°N / 75–0°W download, so Table 1 can be compared exactly. **Comparing on
the larger box gives systematically different frequencies and looks like an
error.**

### D9 — Colson–Panofsky's length scale λ

**Chosen.** λ = the vertical grid increment Δz taken from the **ICAO standard
atmosphere**, i.e. a function of pressure alone and therefore constant in space
and time.

**In line with.** Sharman (2006) A4: CP = λ² S_V² (1 − Ri/Ri_crit), *"where λ
is a length scale, taken as the local value of vertical grid increment Δz, and
Ri_crit is an empirical constant (≈0.5)"*.

**The ambiguity, stated.** On a pressure-level grid the "vertical grid
increment" in metres is not unique: the ICAO reading makes λ constant, while
using the actual geopotential thickness makes it a spatially and temporally
varying field — which would change ranks, not just scale. The constant reading
was chosen and is what makes the factor inert for a percentile calibration.
**This is the same family of ambiguity as D10 and should be reported the same
way.** (Audit §4.3, F8.)

**Implementation note.** CP is computed in the algebraically reduced form
λ²(S_V² − N²/Ri_crit) rather than through a Richardson round trip, because
rojak's `ColsonPanofsky` combines two different S_V² values that are supposed
to be the same number — they disagree by up to ~3000× at low-shear cells — so
the intended cancellation never happens.

### D10 — Endlich's directional shear — **OPEN, needs a decision**

**The only unresolved formula choice.** Not a missing citation: a literature
*silence*.

Sharman (2006) A25: *s |∂ψ/∂z|*, where ψ is the wind direction. Williams &
Storer (2022) Eq. (5) p. 1427 writes the same and specifies only "centred
second-order finite differences".

Two defensible readings of ∂ψ/∂z:
- **rojak's angle secant** (current production): form the wind direction at
  each level, difference it across the 50 hPa layer with a 2π wrap.
- **The closed form** (F10): |u ∂v/∂z − v ∂u/∂z| / (u² + v²), exact for the
  vector angle atan2(v, u), no arctangent and no 0/2π seam.

**Measured difference.** ρ = 0.9636 on the validation subset, median relative
difference 14.6 %, **41–66 % of the p97 exceedance set flips**.

**Why it matters for this dataset.** `endlich`'s 1979 level ratio against
Prosser is **0.51**, one of the worst of the 21, and it is *not* in the
stencil-sensitive group — so D6 does not explain it, and this is the only other
moving part it has. Williams & Storer (2022) Table 1 p. 1431 ranks wind speed ×
directional shear among their **most** resolution-robust diagnostics, which
makes the 0.51 harder to attribute to resolution.

**Recommendation, not yet adopted.** The closed form is the better estimator on
its own terms; the angle secant is the more likely reading of what Prosser's
stack did. Resolve it the way D1 was resolved: measure both on the full series,
report the one matching the method being replicated, appendix the other.
**One A/B run.** F10 was deliberately excluded from the production series for
exactly this reason.

### D11 — Excluded and why

- **F10 / `endlich`** — see D10. Documented sensitivity, not baked in.
- **F2, F3, F11** — retired. They are strict subsets of D1 (F12), and F3's
  variant computed A18 on the *archived* ζ, which is less coherent than
  computing it on the substituted ζ as F12 does.

---

## 2. THE 21 DIAGNOSTICS — formula, source, provenance

Numbering is **Williams & Joshi (2013)**, matching `REFERENCE_TABLE["num"]`.
Equation numbers are Sharman et al. (2006) Appendix A, pp. 282–285, unless
stated. "Computed" in the provenance column means built from u, v, T, z under
D1.

| # | key | Source eq. | Form as implemented | Provenance / notes |
|---|---|---|---|---|
| 1 | `magnitude_pv` | A18, p. 283 | \|PV\|, PV = −g ζ_a ∂θ/∂p | ζ computed; **A18 truncated form, not ERA5 Ertel PV** (D1) |
| 2 | `colson_panofsky` | A4, p. 282 | λ²(S_V² − N²/Ri_crit), Ri_crit = 0.5 | λ = ICAO Δz, constant (D9) |
| 3 | `brown1` | **A13** | (0.3 ζ_a² + D_SH² + D_ST²)^½ | **ζ_a = ζ + f, ABSOLUTE vorticity**; ζ computed |
| 4 | `temperature_gradient` | A23, p. 284 | [(∂T/∂x)² + (∂T/∂y)²]^½ | — |
| 5 | `horizontal_divergence` | A33 | \|∂u/∂x + ∂v/∂y\| | **computed**, not ERA5 `d` (D1) |
| 6 | `vertical_wind_shear` | A3 | [(∂u/∂z)² + (∂v/∂z)²]^½ | 50 hPa stencil (D6) |
| 7 | `endlich` | A25, p. 284 | s·\|∂ψ/∂z\| | **rojak angle secant — D10 OPEN** |
| 8 | `deformation` | A17 | (D_SH² + D_ST²)^½ | **un-squared** (D4) |
| 9 | `wind_speed` | A24, p. 284 | \|v\| | — |
| 10 | `ngm2` | A29 | \|∂T/∂z\|·DEF | 2.8 % truncation error, documented |
| 11 | `negative_richardson` | A1, p. 282 | −Ri, Ri = N²/S_V² | N² = (g/θ)(∂θ/∂z) (A2) |
| 12 | `rva_magnitude` | — | \|v·∇ζ\| | ζ computed; scalar gradient correct here |
| 13 | `ubf` | A30, p. 284 | \|−∇²Φ + 2J(u,v) + fζ − βu\| | ζ computed; **true spherical Laplacian and Jacobian** |
| 14 | `nva` | A37, p. 285 | MAX[−u ∂(ζ+f)/∂x − v ∂(ζ+f)/∂y, 0] | **one-sided clip is in the published formula** |
| 15 | `brown2` | **A14** | (1/24) Φ S_V² | **rank-only** — see §5 |
| 16 | `vorticity_squared` | A21, p. 283 | ζ² = \|∇ × v\|² | ζ computed |
| 17 | `ti1` | A15, p. 283 | S_V · DEF | — |
| 18 | `ngm1` | A28, p. 284 | \|v\| · DEF | — |
| 19 | `ti2` | A16, p. 283 | S_V (DEF − Δ_H) | **sign confirmed**: Ellrod & Knapp (1992) eq (9)/(10) define CVG = −(∂u/∂x + ∂v/∂y), so DEF + CVG = DEF − δ. δ computed (D1) |
| 20 | `f2d` | A9, p. 283 | signed isentropic material derivative, variant A | D3; 0.637 diurnal damping |
| 21 | `ncsu1` | A36, p. 284 | [1/MAX(Ri,10⁻⁵)]·MAX(u ∂u/∂x + v ∂v/∂y, 0)·\|∇ζ\| | ζ computed; **both MAX are in the published formula** |

**Verbatim quotes worth having in the appendix:**

> **A13** — Φ = (0.3 ζ_a² + D_SH² + D_ST²)^(1/2), *"where the shearing
> deformation D_SH = ∂v/∂x + ∂u/∂y, the stretching deformation
> D_ST = ∂u/∂x − ∂v/∂y, absolute vorticity ζ_a = ζ + f, with ζ = ∂v/∂x − ∂u/∂y
> and f is the Coriolis frequency."*

> **A14** — ε = (1/24) Φ S_V²

> **A4** — CP = λ² S_V² (1 − Ri/Ri_crit), *"where λ is a length scale, taken as
> the local value of vertical grid increment Δz, and Ri_crit is an empirical
> constant (≈0.5)."*

> **A18** — PV = −g ζ_a ∂θ/∂p

> **A30** — UBF = −∇²Φ + 2J(u,v) + fζ − βu, *"where Φ is geopotential, J is the
> Jacobian operator, and β is the Coriolis frequency gradient."*

> **A37** — NVA = MAX{[−u ∂(ζ+f)/∂x − v ∂(ζ+f)/∂y], 0}

---

## 3. THE SIGN CONVENTION — verified for all 21

**All 21 diagnostics are one-tailed upper**: exceedance is `value ≥ threshold`,
sign `"+"`. **There is no two-tailed diagnostic among them and no lower-tail
one.**

**Evidence: Williams (2017) Table 2, p. 580.** It tabulates every diagnostic's
onset threshold at all five severities, and **every one of the 21 ladders is
monotonically increasing Light → Severe** — including the two whose thresholds
are negative:

```
                                          Light  L-to-M    Mod   M-to-S  Severe
Negative Richardson number                -15.4    -9.8   -7.9    -6.7    -5.9
Vertical shear of horizontal wind (1e-3)    5.3     6.6    7.4     7.9     8.4
Colson-Panofsky index (1e3 kt^2)          -29.3   -27.0  -25.2   -23.7   -22.2
Frontogenesis function (1e-9 m2 s-3 K-2)    770    1280   1660    1980    2340
Brown index (1e-6 s-1)                       99     106    110     113     118
Brown energy dissipation rate (1e-6)        870    1370   1730    2030    2330
Ellrod TI1 (1e-9 s-2)                       195     292    360     419     472
Ellrod TI2 (1e-9 s-2)                       184     282    356     419     477
Flow deformation (1e-6 s-1)                50.9    60.9   66.9    71.8    76.3
Magnitude of potential vorticity (PVU)     8.33    8.73   8.98    9.19    9.41
Relative vorticity squared (1e-9 s-2)      2.46    3.74   4.70    5.50    6.24
|Horizontal temperature gradient| (1e-6)   14.7    17.6   19.4    20.8    22.0
Wind speed (m s-1)                         40.9    48.4   52.4    55.3    58.5
Wind speed x directional shear (1e-3)      3.21    3.94   4.39    4.72    5.08
Deformation x wind speed (1e-3 m s-2)      1.65    2.29   2.76    3.17    3.54
Deformation x vertical T gradient (1e-9)     53      84    106     127     151
|Residual of nonlinear balance| (1e-12)    1230    1840   2270    2610    2960
|Horizontal divergence| (1e-6 s-1)         11.9    15.7   18.2    20.4    22.5
NCSU1 (1e-18 s-3)                          1200    3600   6300    9300   13000
Negative absolute vorticity adv. (1e-9)    1.33    1.86   2.23    2.56    2.93
|Relative vorticity advection| (1e-9)      1.44    1.99   2.34    2.66    3.00
```

A monotone increasing ladder is a one-tailed upper criterion by definition.
Williams explains the negative Colson–Panofsky values directly: *"the
Colson–Panofsky index is proportional to 1 − Ri/0.5, and the Richardson number
(Ri) is rarely less than 0.5 in the GFDL-CM2.1 model."*

**Why this needed checking.** A flipped sign is invisible in every magnitude
comparison — a ratio against a published median is unchanged — and fatal in
rank, which is all the exceedance counting uses.

**This table is GFDL-CM2.1 at ~2°.** It establishes **sign and units only**;
magnitudes are not comparable and Williams says so.

---

## 4. THE MEASURED RESULT — what the conventions are worth

### 4a. Reach versus magnitude

| | magnitude | diagnostics changed | material (flip ≥ 0.5 %) |
|---|---|---|---|
| D2 (F1) metric | 0.4–0.7 % on ∂/∂y | 13 of 21 | 3 |
| D1 (F12) provenance | 1.6–11 % on the fields | 10 of 21 | 10 |
| D3 (F6) f2d variant | ρ = 0.0275 | 1 | 1 |
| **union** | | **16 of 21 differ** | **13 of 21 material** |

**5 of 21 are bit-identical** under all three: `colson_panofsky`, `endlich`,
`negative_richardson`, `vertical_wind_shear`, `wind_speed` — the pure vertical
and point quantities, with neither a horizontal ∂/∂y nor an archived field.
Confirmed independently by the threshold comparison, which returns **exactly
0.00 %** for those five at all five severities.

**The 16-versus-13 gap belongs in the appendix.** Ten of D2's thirteen moved by
less than the exceedance-set noise floor. *"16 of 21 changed"* invites a reader
to think the dataset was overhauled.

### 4b. The trends are invariant

Annual trend ratio against Prosser, baseline → audit:
**0.96→0.97, 0.96→0.98, 0.98→1.01, 0.98→1.02, 0.99→1.04.** Maximum shift five
hundredths. **25/25 significant and 25/25 containing Prosser** in both.

> The fitted 1979–2020 trends are invariant to the spherical metric
> convention, the vorticity/divergence/PV provenance, and the frontogenesis
> variant — three implementation choices that no paper in the CAT literature
> specifies, that together change 16 of the 21 derived fields, and that move
> the calibration thresholds by up to 8 %. No significance verdict changes and
> no trend ratio moves by more than five hundredths.

### 4c. The levels improved

Annual 1979 exceedance, ours ÷ Prosser:

| | baseline | audit | deficit before → after | gap closed |
|---|---|---|---|---|
| LOG | 0.881 | **0.904** | 11.9 % → 9.6 % | 19 % |
| LMOG | 0.857 | **0.888** | 14.3 % → 11.2 % | 22 % |
| MOG | 0.841 | **0.880** | 15.9 % → 12.0 % | 25 % |
| MSOG | 0.820 | **0.865** | 18.0 % → 13.5 % | 25 % |
| SOG | 0.797 | **0.842** | 20.3 % → 15.8 % | 22 % |

**This was not guaranteed, and the two-sided prediction is on record** (audit
§6 F3: A18 moves `magnitude_pv` into the stencil-sensitive group and may verify
worse). It improved at all five severities, monotonically.

### 4d. Level ratio by season — the D6 fingerprint

| season | LOG | LMOG | MOG | MSOG | SOG | LOG−SOG |
|---|---|---|---|---|---|---|
| DJF | 0.957 | 0.930 | 0.910 | 0.884 | 0.859 | +0.097 |
| MAM | 0.916 | 0.897 | 0.881 | 0.860 | 0.815 | +0.101 |
| **JJA** | 0.843 | 0.836 | 0.845 | 0.857 | 0.861 | **−0.018** |
| SON | 0.897 | 0.882 | 0.874 | 0.860 | 0.820 | +0.077 |
| Annual | 0.904 | 0.888 | 0.880 | 0.865 | 0.842 | +0.062 |

---

## 5. UNITS AND DIMENSIONAL CAVEATS — all constants, all inert

**The governing fact:** thresholds are percentiles of our own data, so for any
strictly increasing transformation, `D' ≥ percentile_p(D')` ⟺
`D ≥ percentile_p(D)`. **A constant multiplicative factor cannot change a
single exceedance decision.** Every item below is therefore a
magnitude-comparison caveat, not a correctness one.

- **`brown2` — the missing length².** A14 gives ε = (1/24)Φ S_V², which is
  **s⁻³**; Williams tabulates 10⁻⁶ J kg⁻¹ s⁻¹ = 10⁻⁶ m² s⁻³. **The dimensional
  gap is in Sharman's own equation**, not in this pipeline. No length scale is
  invented. **This explains `brown2` sitting "five orders" from its published
  value** — it is the missing L², it is expected, and it is harmless. Report
  `brown2` as rank-only; do not compare its magnitude to any published
  J kg⁻¹ s⁻¹ figure.
- **`magnitude_pv` — ×100 and PVU.** A18 differentiates θ against a pressure
  coordinate carried in hPa, giving 100× the SI value, and published tables are
  in PVU = 10⁻⁶ SI on top of that. Neither factor is applied. Divide out the
  ×100 and **archived full Ertel PV and A18 differ by only ~2.8 % at the
  threshold** — against a 64–80 % flip rate. The two definitions produce
  similarly *distributed* values and rank cells very differently: worth saying,
  because 64–80 % alone reads as though the diagnostic was replaced wholesale.
- **`colson_panofsky` — 10³ kt².** Native m² s⁻²; the conversion factor is
  1/(1000 × 0.2646526) = 3.778538 × 10⁻³. Getting this wrong is the entire
  reason CP read "1–2 orders off" in early comparison tables.
- **`f2d` — signed thresholds.** Under D3 the thresholds cross zero, so a
  *relative* difference between two threshold sets is meaningless for `f2d`
  (the same trap as `colson_panofsky`). Read `f2d` from flip rates only.
- **`deformation` — see D4.** Inert for exceedance, **not** inert for EVT.

---

## 6. INDEPENDENT MAGNITUDE ANCHORS

Three published tables, so that agreement with Prosser is not carrying the
whole argument:

| Source | Data | What it establishes | Result |
|---|---|---|---|
| Williams (2017) Table 2 | GFDL-CM2.1 ~2° | **sign and units only** | consistent; **this is the sign evidence, §3** |
| Lee et al. (2023) Table 1 | **ERA5** | closest match in data source | all five within ±23 % |
| Sharman (2006) Table B1 | 20-km RUC | the only table covering the hand-written seven | six of seven within a factor of 3 |

Sharman B1 detail: TI1 1.02, \|∇_H T\| 1.03, UBF 1.85, CP 2.20, NCSU1 2.75,
−Ri 0.60, F2D 35.6. This **cleared UBF and NCSU1**, which looked alarming
against the GCM and are simply resolution-sensitive. F2D's 35.6× is consistent
with the 0.637 diurnal damping plus a 20-km-versus-0.25° gap. **`brown2` is
absent from B1**, which is why §5's caveat matters for it.

---

## 7. STILL OPEN

1. **D10 — `endlich`.** The one unresolved formula choice. Needs a decision and
   one A/B run.
2. **Provenance sweep of the 504 audit stores.** `jobs/15` reads every store's
   `.zattrs`; it hardcodes the baseline directories and has never been run on
   `north_atlantic_audit`. Count, size and one file's attributes are verified;
   all 504 are not.
3. **Per-diagnostic scorecard on the audit series.** §17.5's investigate list
   (`ncsu1` 6.07, `ubf` 0.19, `colson_panofsky` 0.14, …) was computed on the
   **baseline** series and has never been rebuilt. This is the table that says
   which individual indicators can be trusted.
4. **`horizontal_divergence`** — pre-registered: if it acquires a significant
   trend under D1, ERA5 observing-system changes are implicated; if it stays
   flat, that hypothesis loses its best candidate.
5. **Sharman & Pearson (2017) Part I**, JAMC 56, 317–337 — nice-to-have. The
   file in `Articles/` under that name is **Part II** (Pearson and Sharman,
   nowcasting).
6. **Brown (1973)**, *Meteorol. Mag.* **102**, 347–361 — provenance only; A13
   is what the replication chain uses.
7. **Peaks-over-threshold archive for the North Atlantic 1979–2020.** The tail
   archives are the **global reference year only**. EVT work on the NA series
   needs per-gridpoint excesses across 42 years (~23 GB, §6 phase 5). Not
   needed to call the indicator dataset final; needed before the econometrics.

---

## 8. PRIMARY SOURCES

All in `Articles/` unless noted.

- Sharman, R., et al. (2006). *Wea. Forecasting* **21**, 268–287. Appendix A
  (the diagnostic definitions), Appendix B Table B1 (magnitude anchors).
- Williams, P. D. (2017). *Adv. Atmos. Sci.* **34**, 576–586. Table 1 (the
  percentile ladder), Table 2 p. 580 (**the sign evidence**).
- Williams, P. D., and M. M. Joshi (2013). *Nature Clim. Change* **3**,
  644–648. Table 1 (the diagnostic numbering and published medians).
- Williams, P. D., and L. N. Storer (2022). *Q. J. R. Meteorol. Soc.* **148**,
  1424–1435. Eq. (3) p. 1427 (**D3**), Eq. (5) p. 1427 (**D10**), Eq. (7)
  p. 1428, Table 1 p. 1431 (resolution sensitivity).
- Prosser, M. C., et al. (2023). *Geophys. Res. Lett.* **50**. p. 2 (**D1**,
  the four fields), Figure 2 caption (**D8**, the box), Figure 4 (per-diagnostic
  values), Supporting Information (no threshold tables).
- Ellrod, G. P., and D. I. Knapp (1992). *Wea. Forecasting* **7**, 150–165.
  Eq. (9)/(10) — the TI2 convergence sign.
- Koch, S. E., and F. Caracena (2002) — the UBF residual form.
- Storer, L. N., et al. (2017) — dropped \|PV\| from this basket "because it
  was found to give unrealistic results". Relevant context for D1's
  `magnitude_pv` result.
- Lee, J. H., et al. (2023) — ERA5 magnitude anchor; Eq. 6 p. 3 is A18's
  hydrostatic equivalent.
- Wong, W. K., et al. (2025) — rojak's author; uses 175/200/225 at 200 hPa
  (**D6**).
- Jaeger, E. B., and M. Sprenger (2007). *J. Geophys. Res.* **112**. Reference
  list gives Brown (1973) as *Meteorol. Mag.* 102, **347–361**.
- **IFS Documentation CY49R1, Part III**, §2.2.7 (**D1**, ζ and D are
  prognostic). Not in `Articles/`; ECMWF, open access.
- **MetPy** `src/metpy/calc/tools.py` — `nominal_lat_lon_grid_deltas`,
  `parse_grid_arguments`, `geospatial_gradient`, `vector_derivative`
  (**D2**). Not a paper; cite the version.
- **rojak** pinned at rev `25b8685c670401883bf6d186a522ccfd4561c908`. The pin
  is the reproducibility record.
