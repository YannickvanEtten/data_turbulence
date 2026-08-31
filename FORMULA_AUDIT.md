# FORMULA AUDIT — the 21 diagnostics against their source equations

Created 2026-08-29. Companion to `STATUS.md` §4 and §7 and to
`CALIBRATION_REFERENCE.md`.

> **§10 was added the same evening and reports the measured results. Where it
> contradicts §4.4, §7.2 or §9 below, §10 wins** — two of the claims in those
> sections were read off a colour figure by eye and did not survive checking.
> The rest of the document stands.

**Division of labour between the three documents.** `STATUS.md` §4 verifies
*code against docstring* (does the implementation compute the formula it claims
to?). `CALIBRATION_REFERENCE.md` compares *magnitudes* against published
tables. This file closes the gap between them: *docstring against literature* —
does each diagnostic compute the formula the source paper actually states?

Sources read for this audit: Sharman et al. (2006) Appendix A (eqs A1–A37) and
Appendix B Table B1, in full; rojak's own source at the pinned rev `25b8685`
for all 14 diagnostics it supplies; Williams (2017) §2 and Table 2; Williams &
Joshi (2013) Table 1; Prosser (2023) Supporting Information Figure S5.

---

## 1. Headline

**Eighteen of the 21 are now closed against the literature.** Nine `STATUS.md`
§7 items that were open this morning are closed below, four of them by reading
equations that were already on disk in `Articles/`.

**One diagnostic does not match its source and is the highest-value thing in
this document: `f2d`.** Two independent lines of evidence — the printed form of
A9, and the shape of the published frontogenesis distribution — say our
implementation is not the quantity Williams and Prosser computed. §4.

**Two more are flagged, both of which the 39 months already on disk can
settle:** `ubf` (formula confirmed correct, so the §12.6 anomaly is elsewhere)
and `deformation` (saved squared — harmless for the replication, not harmless
for phase 5).

> **Measured 2026-08-29 evening, §10.** `f2d` is confirmed as the single
> outlier and variant A is refuted quantitatively. `ubf` is **cleared** — none
> of the three proposed mechanisms survives measurement, and the "fragile
> residual" framing repeated from `STATUS.md` §4g is itself wrong: the residual
> is 73 % of the largest term, not a near-cancellation. `deformation` is fixed
> and confirmed three ways. The Figure S5 readings in §4.4 and §7.2 were done
> by eye and one of them compared a relative change against a map of absolute
> change; §10.6 corrects both.

---

## 2. What Sharman Appendix A actually says, term by term

Equation numbers are Sharman's. "rojak" means the implementation at rev
`25b8685`, read directly, not its docstring.

| # | Diagnostic | Source eq. | Source statement | Implementation | Verdict |
|---|---|---|---|---|---|
| 1 | `magnitude_pv` | A18 | \|PV\| | `np.abs(pv)`, ERA5's own PV | ✅ |
| 2 | `colson_panofsky` | A4 | CP = λ²S_V²(1 − Ri/Ri_crit), λ = Δz, Ri_crit ≈ 0.5 | `(Δz)²(S_V² − N²/Ri_crit)` — algebraically identical, since S_V²·Ri = N² | ✅ |
| 3 | `brown1` | A13 | Φ = (0.3ζ_a² + D_SH² + D_ST²)^½ | `sqrt(0.3·ζ_a² + D_sh² + D_st²)`, ζ_a = ζ + f | ✅ **0.3 confirmed** |
| 4 | `temperature_gradient` | A23 | \|∇_H T\| = [(∂T/∂x)² + (∂T/∂y)²]^½ | `magnitude_of_geospatial_gradient(T)` | ✅ |
| 5 | `horizontal_divergence` | A33 | Δ_H = ∂u/∂x + ∂v/∂y | `np.abs(divergence)`, ERA5's own | ✅ |
| 6 | `vertical_wind_shear` | A3 | S_V = \|∂**v**/∂z\| = [(∂u/∂z)² + (∂v/∂z)²]^½ | same | ✅ |
| 7 | `endlich` | A25 | s·\|∂ψ/∂z\|, ψ = wind direction | `wind_speed × abs(dψ/dz)`, angle-aware gradient | ✅ |
| 8 | `deformation` | A17 | DEF = (D_SH² + D_ST²)^½ | rojak's DEF **diagnostic** returns DEF². See §5. | ⚠️ |
| 9 | `wind_speed` | A24 | s = \|**v**\| | same | ✅ |
| 10 | `ngm2` | A29 | \|∂T/∂z\|·DEF | same, with un-squared DEF | ✅ |
| 11 | `negative_richardson` | A1, A2, A3 | Ri = N²/S_V² | hand-written; `N²/S_V²` | ✅ |
| 12 | `rva_magnitude` | — | not in Sharman; W&J-specific | \|u∂ζ/∂x + v∂ζ/∂y\| | ✅ (no source eq. to check) |
| 13 | `ubf` | A30 | UBF = −∇²Φ + 2J(u,v) + fζ − βu | code computes −(A30) then takes \|·\| → identical | ✅ **sign closed** |
| 14 | `nva` | A37 | max{[−u ∂(ζ+f)/∂x − v ∂(ζ+f)/∂y], 0} | rojak uses **absolute** vorticity ζ+f and clips at 0 | ✅ **clipping closed** |
| 15 | `brown2` | A14 | ε = (1/24)Φ S_V² | `(1/24)·brown1·S_V²` | ✅ (see §6 for the units) |
| 16 | `vorticity_squared` | A21 | ζ² = \|∇×**v**\|² | ERA5's ζ, squared | ✅ |
| 17 | `ti1` | A15 | TI1 = S_V·DEF | same, un-squared DEF | ✅ |
| 18 | `ngm1` | A28 | \|**v**\|·DEF | same, un-squared DEF | ✅ |
| 19 | `ti2` | A16 | TI2 = S_V(DEF − Δ_H) | `vws·(DEF + (−divergence))` | ✅ |
| 20 | `f2d` | A9 | see §4 | **does not match A9 as printed** | ❌ |
| 21 | `ncsu1` | A36 | [1/max(Ri, 10⁻⁵)]·max(u∂u/∂x + v∂v/∂y, 0)·\|∇ζ\| | identical, including both clips | ✅ **clipping closed** |

### 2.1 `STATUS.md` §7 items closed by this table

- **`brown1`'s 0.3 coefficient.** §7 says Brown (1973) *Meteorological Magazine*
  102 "is not in `Articles/` and must be obtained". It does not need to be.
  **Sharman A13 prints the coefficient directly**: Φ = (0.3ζ_a² + D_SH² +
  D_ST²)^½, with ζ_a defined in the same paragraph as absolute vorticity ζ + f.
  rojak matches exactly. Closed — no PDF hunt required.
- **`nva` and `ncsu1` clipping.** §7 asks whether A36 and the NVA definition
  really are one-sided. **Both are, explicitly.** A37 is written with an outer
  `MAX{·, 0}` and A36 with an inner `MAX(·, 0)`; A36 also prints the
  `MAX(Ri, 10⁻⁵)` floor that the code implements. Confirmed independently by
  Sharman Table B1, where NCSU1's null threshold T1 is exactly `0.0` — only
  meaningful for a quantity clipped at zero. Closed.
- **`ubf`'s sign.** A30 reads UBF = −∇²Φ + 2J(u,v) + fζ − βu. The code computes
  `∇²Φ − 2J − fζ + βu` and then takes the absolute value, so it is −1 × A30
  under a magnitude — identical. Closed.
- **Every `sign` entry in `REFERENCE_TABLE`.** §7 called this "a cheap,
  high-value re-derivation". Sharman Table B1 does it for free, because it
  prints the null→extreme threshold ladder for the upper levels and the
  direction of that ladder *is* the sign convention:
  - `−Ri`: −20, −2.0, −0.6, −0.3, **+0.5** — monotone increasing, one-tailed
    upper. Confirms `"+"`. (Note T5 > 0: severe requires Ri < −0.5, i.e. static
    instability.)
  - `CP`: **0**, 1000, 5000, 12000, 30000 kt² — increasing from zero, one-tailed
    upper. Confirms `"+"`, and independently of the Williams (2017) argument
    already recorded in `CALIBRATION_REFERENCE.md` §4.2.
  - `NCSU1`: 0.0, 10⁻¹³, 3.5×10⁻¹³, 1.5×10⁻¹², 4×10⁻¹² — same.
  - `UBF`, `TI1`, `|∇_H T|`, `F_θ`: all positive and increasing.
  All 21 `"+"` entries stand. Closed.

### 2.2 One correction to `STATUS.md` §7's framing

§7 lists CP's *scale factor* as having been wrong (`SCALE_TO_TABLE`, §11.5) but
its unit right. Sharman Table B1 prints CP in **kt²**, Williams (2017) Table 2
in **10³ kt²**. Both are in circulation; `REFERENCE_TABLE` uses Williams'
10³ kt², which is the right choice since Williams is the calibration lineage.
Worth a one-line comment in the code so the factor of 1000 is never
re-litigated.

---

## 3. What rojak's 14 actually compute

The §4 analytic suite in `STATUS.md` verifies rojak's output against
hand-derived expected values. That establishes rojak computes *its docstring*.
This audit read rojak's source at `25b8685` to check the docstrings themselves,
and to check the wiring — which is where a factory can silently pass the wrong
field.

**The wiring is correct.** The one thing worth checking was whether the
`total_deformation` handed to TI1, TI2, NGM1 and NGM2 is DEF or DEF²; if it
were squared, those four would be `S_V·DEF²` rather than `S_V·DEF`, which is
**not** a monotone transform and would change ranks. It is not squared:
`CATData.total_deformation()` calls `magnitude_of_vector(..., is_squared=False)`,
and only the standalone `DeformationSquared` diagnostic squares it. All four
products are correct. See §5 for the consequence that *does* remain.

`NegativeVorticityAdvection` was the other wiring risk, since advecting relative
rather than absolute vorticity drops the −vβ term, which at v ≈ 10 m s⁻¹ is
≈ 1.6 × 10⁻¹⁰ s⁻² — the same order as W&J's published NVA median of
2.05 × 10⁻¹⁰. It uses `absolute_vorticity(ζ)`. Correct.

---

## 4. `f2d` — the one that does not match

### 4.1 What A9 says

Sharman defines the frontogenesis function as F = (D/Dt)|∇θ|, then writes:

> Expanding on a constant-θ surface and invoking continuity gives
>
> F_θ ∝ **−** D/Dt [ (∂u/∂θ)² + (∂v/∂θ)² ]^½
>   = |∂**v**/∂θ|⁻¹ [ (∂u/∂θ)·D/Dt(∂u/∂θ) + (∂v/∂θ)·D/Dt(∂v/∂θ) ]    (A9)

and immediately after, the constant-pressure form used at midlevels:

> F_p ∝ (∂u/∂p)·D/Dt(∂u/∂p) + (∂v/∂p)·D/Dt(∂v/∂p)    (A10)

We compute

```python
a9 = 0.5 * D/Dt[ (du/dθ)² + (dv/dθ)² ]
```

**That is A10's algebraic form evaluated in θ rather than p.** It is not A9.
The two differ in two ways:

1. **A normalisation.** A9 carries a `|∂v/∂θ|⁻¹` factor; ours does not. Since
   D/Dt(Q^½) = D/Dt(Q)/(2Q^½), our quantity equals A9's times |∂**v**/∂θ| — a
   *spatially varying field*, not a constant. **It changes ranks**, which is
   `STATUS.md` §5's category-1 error, the one the percentile calibration cannot
   absorb.
2. **A sign.** A9's left-hand side carries a leading minus that its own
   right-hand side does not. That is an internal inconsistency in the printed
   paper — the two sides are not equal as written. The surrounding text ("…
   invoking continuity gives F_θ ∝ −D/Dt…") puts the minus on the physics, so
   the RHS is most likely the typo. Physically the minus is right: on an
   isentropic surface, strengthening ∇θ corresponds to *weakening* |∂**v**/∂θ|,
   so frontogenesis is a *decrease* of the bracket. **Our implementation has
   no minus.** If the minus belongs, we are flagging frontolysis, not
   frontogenesis — a disjoint set of grid cells.

### 4.2 Why the units argument settles less than `STATUS.md` §7 thinks

§7 records this as closed: *"`f2d`'s units argument — CONFIRMED. Sharman's own
Table B1 units require the squared form, which contradicts the printed `^(1/2)`
in A9. The table wins."*

Table B1 does support the squared form, and this audit confirms it: B1 lists
F_θ (A9) in **m² s⁻³ K⁻²**, which are the units of our un-normalised
expression, whereas literal A9 with its `|∂v/∂θ|⁻¹` has units m s⁻² K⁻¹. B1 also
lists F_p (A10) in m² s⁻³ Pa⁻² — the same algebraic shape one coordinate over,
which is exactly the pattern our code follows. So the *algebra* is defensible.

**But units are invariant to sign.** The conclusion "the table wins" closes
question (1) and says nothing at all about question (2). §7 already notes the
leading minus separately as unchecked; this audit confirms it is still open and
raises its priority, because it is a §5 category-2 error — invisible in every
magnitude comparison made so far and fatal in rank.

### 4.3 A third, independent piece of evidence — and it is the strongest

Williams (2017) and Williams & Joshi (2013) publish, for the *same* model, box,
season and sampling (GFDL-CM2.1, DJF, 200 hPa, 50–75°N/10–60°W, **daily-mean
fields** — confirmed in Williams 2017 §2), both a median and a p97 threshold for
every diagnostic. Their ratio is dimensionless and therefore transfers across
models far better than either number alone. For frontogenesis:

    p97 / median = 770 / 56.6 = 13.6

**A signed material derivative cannot produce that.** D/Dt of a positive
quantity in a statistically stationary atmosphere is centred on zero; its median
is near zero and the ratio diverges. A ratio of 13.6 is what a *positive-definite*
quantity looks like — and it fits the family pattern of the other 20 exactly
(§7.1). Our own most recent full-field statistics (`cat_outputs/summary_stats.csv`,
January 1979, North Atlantic) give p99/median ≈ 101 for `f2d` against
p99/median ≈ 6–15 for its algebraic peers. Our distribution is centred; theirs
is not.

That is not a resolution effect and it is not a units effect. It says Williams'
frontogenesis diagnostic is one-signed and ours is two-signed.

### 4.4 And a fourth: Prosser's Figure S5

S5 (extracted from the Supporting Information, `image10.jpeg`) maps the MOG
absolute change for each of the 21. The **frontogenesis function** panel is
near-white to slightly *negative* over the North Atlantic and strongly negative
in the tropics. Our `f2d` is the **second-largest positive relative trend** of
the 21 in `STATUS.md` §12.6 (+226 %). Sign and pattern both disagree.

### 4.5 What to do about it

Three candidate forms, all cheap to compute from the 39 months already on disk,
none needing a single new CDS request:

| variant | expression |
|---|---|
| **A** (current) | `+½ D/Dt[(∂u/∂θ)² + (∂v/∂θ)²]` |
| **B** (A9 sign) | `−½ D/Dt[…]` |
| **C** (magnitude) | `\|½ D/Dt[…]\|` |

Compute all three on the global year-2000 calibration set, and for each report
(i) the median and the p97/median ratio, against the published 13.6, and (ii)
the resulting DJF 1979 → 2020 MOG change over Prosser's box, against S5. Variant
**C** is the one predicted to reproduce a finite ratio near 13.6; if it also
lands closer to S5's near-zero North Atlantic change than variant A does, that
is three independent lines agreeing and the question is settled.

Until then `f2d` should be reported with a caveat, and it is worth checking how
much the 21-member ensemble mean moves when it is excluded — with 21 members, a
single misbehaving diagnostic contributes about 5 % of the ensemble, so this is
almost certainly *not* what drives the §11 and §12 headline results. It matters
for the per-diagnostic story, not for the replication.

---

## 5. `deformation` is saved squared

rojak's `DEF` diagnostic returns DEF² (`DeformationSquared`). `3_pipeline.py`
knows this and applies `sqrt` in `PRETRANSFORM` — **but only on the path that
builds the comparison table.** `compute_all_21()` does not, so
`ada/diagnostics_global.py` writes DEF² into every production zarr.

- **Harmless for the replication.** Squaring is strictly increasing on a
  non-negative field, so by `STATUS.md` §5.5 the exceedance field is *identical*.
  Nothing in §11 or §12 is affected.
- **Not harmless for phase 5.** §6's plan is peaks-over-threshold and GPD
  fitting on retained magnitudes. A GPD fitted to DEF² is not a GPD fitted to
  DEF — the shape parameter ξ doubles under squaring. Every magnitude-based
  result for this diagnostic would silently be about the wrong variable.
- **Fix it before phase 4, not after.** It is a one-line `np.sqrt` in
  `compute_all_21`, and the alternative is re-deriving 287 GB. The same audit
  should confirm nothing downstream assumes the squared form.

---

## 6. `brown2`'s five orders are a length scale, and they are harmless

`STATUS.md` §7 lists `brown2` as sitting "five orders from its published value
in a way that resolution does not explain". It is explained.

Sharman A14 gives ε = (1/24)Φ S_V², with Φ in s⁻¹ and S_V² in s⁻², so **A14 is
dimensionally s⁻³**. Williams & Joshi tabulate the Brown energy dissipation rate
in J kg⁻¹ s⁻¹ = **m² s⁻³**. The published quantity is therefore A14 times a
length squared, and the code's docstring already says so ("no length scale
invented").

Solving for that length from the measured gap:

| comparison | implied L² | implied L |
|---|---|---|
| vs Williams (2017) MOG threshold (`CALIBRATION_REFERENCE.md` §10.6) | 1.32 × 10⁶ m² | **1147 m** |
| vs W&J (2013) median, Jan-1979 NA | 2.33 × 10⁶ m² | **1526 m** |

**L is of order 1–1.5 km — a model vertical grid increment.** That is the same
λ = Δz that appears explicitly in A4 for Colson–Panofsky, and it is the right
order for GFDL-CM2.1's spacing near 200 hPa. The missing factor is not a formula
error; it is the same length scale GTG uses elsewhere, which W&J applied and
Sharman's A14 leaves out.

It is a **constant**, so by §5.5 it cancels exactly from the exceedance field.
Combined with §12.6 (level 2nd of 21, trend +52 %, both mid-pack) and S5 (Brown
energy dissipation rate is strongly positive over the North Atlantic, matching
ours), `brown2` should be moved out of the §7 suspect list. Note for the
write-up: if a magnitude for `brown2` is ever reported rather than a rank,
multiply by Δz² and say so.

---

## 7. Route 2 — what the other papers can and cannot do, and a better use of them

### 7.1 The shape-ratio test — the strongest free check available

`CALIBRATION_REFERENCE.md` §4.1 correctly concludes that the two Williams
tables cannot anchor magnitude, because a 2° GCM and 0.25° ERA5 differ in every
gradient. It then sets them aside. **That gives up too much.** The two tables
are drawn from the *same* model, box, season and daily-mean sampling, so their
**ratio** is a pure distribution-shape statistic — dimensionless, and far more
transferable than either number.

Computed for all 21 (Williams 2017 Table 2 "light" = p97; W&J 2013 Table 1 =
median):

| diagnostic | p97/median | | diagnostic | p97/median |
|---|---|---|---|---|
| `ncsu1` | 108.1 | | `ti1` | 6.19 |
| `f2d` | 13.60 | | `rva_magnitude` | 6.18 |
| `vorticity_squared` | 11.13 | | `horizontal_divergence` | 4.22 |
| `ubf` | 7.64 | | `endlich` | 3.37 |
| `brown2` | 7.50 | | `vertical_wind_shear` | 2.82 |
| `ngm1` | 6.57 | | `wind_speed` | 2.74 |
| `nva` | 6.49 | | `deformation` | 2.74 |
| `ngm2` | 6.49 | | `temperature_gradient` | 2.56 |
| `ti2` | 6.39 | | `brown1` | 1.28 |
| | | | `magnitude_pv` | 1.22 |

*(`negative_richardson` 0.12 and `colson_panofsky` 0.84 are excluded: both are
negative-valued, so a ratio is not a shape statistic. Use
(p97 − median)/IQR for those two.)*

**The table has structure, and the structure is the test.** The ratio tracks
algebraic degree almost exactly:

- degree 1 (`wind_speed`, `deformation`, `temperature_gradient`,
  `vertical_wind_shear`): **2.5–2.8**
- degree 2 — every product of two degree-1 quantities (`ti1`, `ti2`, `ngm1`,
  `ngm2`, `nva`, `rva_magnitude`): **6.2–6.6**, and 2.7² = 7.3 ✓
- `vorticity_squared` (a square): 11.1 ✓
- `brown2` = `brown1` × S_V²: 1.28 × 2.82² = 10.2, measured 7.50 ✓
- `ncsu1` (cubic, and clipped at zero, which crushes the median): 108 ✓

This is `CALIBRATION_REFERENCE.md` §5.1's composition-ratio idea generalised
from medians to distribution shape — and it is strictly stronger, because it
constrains the tail rather than the bulk, which is where the econometrics is
going to live.

**How to run it, using only what is on disk.** Take `derived/global/` (26 GB,
12 months, year 2000), average to **daily means** to match Williams' sampling,
subset to **50–75°N / 10–60°W** and **DJF** to match his box, and compute the
median and p97 of each of the 21. Report our ratio against the table above. No
downloads, one short job.

**What it will show, predicted in advance so the test can fail:**
- Most diagnostics should land within a factor of ~2 of the published ratio.
- The shear-based ones should come in **systematically low**, because the 50 hPa
  stencil damps the tail (§11.3) — giving a *second, independent* measurement of
  the stencil effect, from distribution shape rather than from exceedance
  frequency.
- `f2d` should come in **far too high** (hundreds, not 13.6) if §4.3 is right.

### 7.2 Prosser's Figure S5 — 21 targets where there was 1

S5 is extracted and readable (`image10.jpeg` inside the Supporting Information
`.docx`; it is the last of the ten embedded images). It maps the MOG **absolute
change in percentage points**, globally, one panel per diagnostic — the only
published per-diagnostic replication target that exists.

To compare, convert `STATUS.md` §12.6's *relative* changes to *absolute*
percentage-point changes, then rank the 21 and compare that ranking against the
panels over the North Atlantic. Two caveats to respect: S5 is annual and global
while §12.6 is DJF over Prosser's box, and S5's colour scale saturates at
±0.5 pp while our DJF `vertical_wind_shear` change alone is +1.14 pp. So this is
a **sign-and-ranking** test, not a magnitude test. It is still 21 constraints
where the project currently has one.

A first read of the panels over the North Atlantic, against §12.6:

| | S5 over the NA | ours (§12.6) | |
|---|---|---|---|
| VWS, TI1, TI2, DEF, NGM1, NGM2, brown1, **brown2**, \|∇T\|, ζ², NVA, RVA, endlich | clearly positive | positive, mid-to-high | ✅ agree |
| `magnitude_pv` | negative globally, **positive patch over the NA/Europe** | +132 % | ✅ agree |
| `negative_richardson` | ≈ 0 at midlatitudes | +176 %, but on a 0.022 % base → +0.039 pp | ✅ agree in absolute terms |
| **`ubf`** | **strongly positive over the NA and Europe** | **+11 %, last of 21** | ❌ **disagree** |
| **`f2d`** | ≈ 0 to slightly negative over the NA | +226 %, 2nd of 21 | ❌ disagree — see §4 |
| `colson_panofsky` | ≈ 0 to slightly negative over the NA | +148 % | ⚠️ check |

**This independently corroborates `STATUS.md` §12.6's suspicion about `ubf`, and
sharpens it.** §12.6 flags UBF as anomalous because it breaks the internal
pattern. S5 says more: UBF is not merely supposed to be *somewhere* in the pack,
it is supposed to be near the *top* over exactly the box we are measuring. Since
§2 of this document confirms the formula now matches A30 exactly, the fault is
not in the equation. The remaining candidates are the inputs and the
near-cancellation itself: ERA5's own `vo` versus the ζ implied by our
`vector_derivatives`, and float32 precision in a residual of terms that cancel
to a few parts in a thousand.

### 7.3 Pointwise composition identities — free, and stronger than the median version

`CALIBRATION_REFERENCE.md` §5.1 tests TI1 = VWS × DEF by comparing *median
ratios* and gets agreement to 0.15 %. That is a good smell test, but medians of
products are not products of medians, so it can only ever be approximate. The
saved zarr allows the exact version, cell by cell:

| identity | must equal |
|---|---|
| `ti1 / (vertical_wind_shear × √deformation)` | 1 |
| `ti2 / (vertical_wind_shear × (√deformation − δ))` | 1 |
| `ngm1 / (wind_speed × √deformation)` | 1 |
| `ngm2 / (\|∂T/∂z\| × √deformation)` | 1 |
| `brown2 / (brown1 × S_V²)` | 1/24 |

*(the `√` because of §5 — which is also why running this test would have caught
the squared-DEF issue immediately.)*

Any deviation beyond float32 rounding means a mis-wired argument, a level
mismatch, or a chunk-boundary error. This is minutes of compute on one month and
it verifies five diagnostics against fields that are themselves already
verified. `brown2` in particular becomes fully closed: A14 confirmed by §2,
composition confirmed pointwise, magnitude explained by §6.

### 7.4 The level question — one premise to correct, and one to keep

> *"I have a bit of a higher level around 200 than Prosser… the general finding
> should not depend on that tweak right? and others are likely to use my levels"*

**The second half is right, and stronger than stated.** 200 hPa is not a
deviation from the literature — it *is* the literature. Williams & Joshi (2013)
and Williams (2017) are both at 200 hPa; Lee et al. (2023) is at 250. Prosser's
188/197/206 are ERA5 *model* levels 73/74/75 and are the outlier
(`CALIBRATION_REFERENCE.md` §11). So for every cross-paper check in §7.1–7.2 the
level is an asset, and the only place it costs anything is the Prosser
replication specifically. That is worth saying plainly in the write-up rather
than apologising for it.

**The first half needs qualifying.** "The general finding should not depend on
that tweak" is true for the *qualitative* result and false for the
*quantitative* one, and `STATUS.md` §12.5 has already measured which is which:

| | level ratio vs Prosser | trend ratio vs Prosser |
|---|---|---|
| LOG | 0.87 | 1.03 |
| MOG | 0.79 | 1.16 |
| SOG | 0.69 | **1.43** |

Sign, severity ordering and significance pattern are robust — those do not
depend on the stencil. But the fitted *relative* trend is inflated by the
stencil, monotonically in severity, reaching **+43 % at SOG**. The mechanism is
understood (a wider stencil damps the tail, so we sit lower on a convex
exceedance curve, so the same physical shift produces a larger relative change).

That matters here more than it would in most projects, because §1's stated
direction is **tail-behaviour modelling**. The stencil bias is largest precisely
in the tail the econometrics wants to model, and it is not a nuisance constant —
it is severity-dependent. It should be carried as a quantified bias, not
dismissed as a tweak.

### 7.5 The stencil experiment, for free

`CALIBRATION_REFERENCE.md` §11.1 proposes downloading 150/175/200/225/250 to get
two points on the stencil-width curve — 15 CDS requests and about a day.

**There is a version that costs no downloads at all.** The existing files
already hold three levels. Recompute the vertical derivatives **one-sided across
200→175 (25 hPa)** instead of centred across 225→175 (50 hPa). That is a second
stencil width from the same bytes, and 25 hPa is much closer to Prosser's 18 hPa
than 50 hPa is.

- **Cost:** a stencil option in `altitude_derivative_on_pressure_level`'s
  callers, then re-run the 12 global calibration months and the 27 NA months —
  about a day of ADA time at the measured rates, and **zero** CDS requests.
- **What it buys:** if the MOG level ratio moves from 0.76 toward 1.0 as the
  stencil narrows, §11.3's mechanism is confirmed rather than merely plausible,
  and the slope lets Prosser's 18 hPa be *extrapolated*. "Why are your levels
  different" stops being an apology and becomes a number.
- **Caveats to state:** a one-sided difference is O(h) rather than O(h²)
  accurate, and its effective evaluation height is nearer 187 hPa than 200. It
  is a **sensitivity probe**, not a better estimate — do not adopt it as the
  production configuration. The five-level download in §11.1 remains the clean
  version if a referee pushes.

---

## 8. Sequencing — this does not gate the download

The single most useful scheduling fact: **the raw download is
diagnostic-independent.** `01_download.sbatch` fetches u, v, t, z, d, vo, pv at
175/200/225 and knows nothing about the 21. Every question in this document
affects only the *derived* stage.

Measured costs from `STATUS.md` §13.2: download **~36 h**, CDS-bound and not
compressible; diagnostics **~7 h**, and re-runnable.

So the answer to *"should we get better insight before we do this all on the 40
years?"* is **yes for the derived stage, no for the raw stage** — and the two
can run at the same time:

1. **Start the 504-month download now.** It is 36 hours of wall clock that
   nothing here can shorten, it is resumable, and no finding below can invalidate
   a byte of it.
2. **While it runs**, on the 39 months already on disk:
   - §7.3 pointwise composition identities — minutes.
   - §7.1 shape-ratio test across all 21 — one short job.
   - §7.2 S5 per-diagnostic comparison — no compute, just the conversion to
     absolute percentage points and a careful look.
   - §4.5 the three `f2d` variants.
   - §7.2 the `ubf` investigation: ERA5 `vo` versus derived ζ, and float32 in
     the residual.
3. **Decide before submitting the diagnostics array** (not before the download):
   the `f2d` variant, the `deformation` √ fix (§5), and whether to persist
   magnitudes for phase 5 (`STATUS.md` §6).
4. **§7.5's stencil probe** afterwards, as a separate deliverable — it answers a
   referee, not a correctness question, and nothing depends on it.

The one genuinely blocking item is unchanged and is not in this document: the
`unlimited` QOS question to `itvo.ucit@vu.nl` (`STATUS.md` §8).

---

## 9. Revised `STATUS.md` §7 suspect list

| was | now |
|---|---|
| `ubf` — top of list, empirical | **stays top**, and S5 (§7.2) makes it worse: it should be near the *top* of the 21 over the NA, not the bottom. Formula confirmed correct against A30, so look at inputs and precision. |
| `f2d`'s leading minus sign | **promoted**. Not one issue but three: sign, the missing `\|∂v/∂θ\|⁻¹`, and a published distribution shape that ours does not reproduce. §4. |
| `brown1`'s 0.3 coefficient | **closed** — Sharman A13 prints it. §2.1. |
| `nva` / `ncsu1` clipping | **closed** — A36 and A37 print the clips; Table B1's `T1 = 0.0` confirms. §2.1. |
| `brown2`'s five orders | **closed** — an implied Δz² with Δz ≈ 1.1–1.5 km, constant, cancels. §6. |
| every `sign` entry | **closed** — Sharman Table B1's threshold ladders. §2.1. |
| — | **new:** `deformation` is persisted squared. Harmless now, wrong for phase 5. §5. |
| — | **new:** `colson_panofsky`'s +148 % trend against a near-zero S5 panel. §7.2. |

---

## Sources

- Sharman, R., C. Tebaldi, G. Wiener, & J. Wolff (2006). An integrated approach
  to mid- and upper-level turbulence forecasting. *Weather and Forecasting*, 21,
  268–287. **Appendix A eqs A1–A37 and Appendix B Table B1**, read in full.
- Williams, P. D. (2017). *Advances in Atmospheric Sciences*, 34(5), 576–586.
  §2 (daily-mean sampling) and Table 2.
- Williams, P. D., & Joshi, M. M. (2013). *Nature Climate Change*, 3, 644–648.
  Table 1.
- Prosser, M. C., et al. (2023). *GRL*, 50, e2023GL103814. Supporting
  Information **Figure S5**.
- rojak, ImperialCollegeLondon/rojak at rev `25b8685c670401883bf6d186a522ccfd4561c908`
  — `src/rojak/turbulence/diagnostic.py`, `src/rojak/turbulence/calculations.py`,
  `src/rojak/core/data.py`, read directly.

---

## 10. MEASURED — 2026-08-29 evening

Jobs `cat-tests`, `cat-audit`, `cat-ubf`, `cat-f2d` on `node013`, repo
`2f711f5`. Logs in `logs/tests-*.out`, `audit-*.out`, `ubf-*.out`, `f2d-*.out`.
Total cost: about nine minutes of compute, no CDS requests, alongside a running
production download.

### 10.1 Headline

**Twenty of the 21 have distribution shapes consistent with the published
ones.** Every composition identity holds pointwise. `f2d` is the single
outlier and the evidence against its current form is now quantitative rather
than interpretive. `ubf` is **cleared** on all three mechanisms proposed for
it — and two claims made earlier in this document, both read off Figure S5 by
eye, do not survive being checked. §10.6.

### 10.2 Composition identities — all five PASS

On `diagnostics_na_1979-01.zarr`, all 9,032,408 cells, no sub-sampling:

| identity | median rel. err | p99 | |
|---|---|---|---|
| I1a `ti1/vws == ngm1/wind_speed` | 2.54e-08 | 8.35e-08 | PASS |
| I1b `deformation == (ti1/vws)^2` | 4.24e-08 | 1.28e-07 | PASS |
| I2 `brown2/(brown1·vws²) == 1/24` | 4.95e-08 | 1.75e-07 | PASS |
| I3 `\|(ti1-ti2)/vws\| == \|divergence\|` | 6.78e-08 | 4.58e-06 | PASS |
| I4 `(ζ_a-f)² == vorticity_squared` | 2.31e-07 | 1.24e-05 | PASS |

All at the float32 floor, in the tail as well as the bulk. `I5` gives a median
`|∂T/∂z|` of **1.64e-03 K/m**, exactly the order expected near the tropopause.

**`brown2` is now fully closed**: A14 confirmed against the paper (§2), its
published magnitude explained as an implied Δz² (§6), and its composition
confirmed pointwise here. It comes off the §7 suspect list.

**The DEF² finding is confirmed on real data.** I1b fits the squared form and
not the un-squared one, independently of reading rojak's source.

### 10.3 Shape ratios — 20 of 21 in family

**The daily-mean column is the valid comparison**, and the instantaneous run
proves why rather than merely asserting it. Two diagnostics are clipped at zero
(`nva` via A37, `ncsu1` via A36), so instantaneously their median sits at zero
and the ratio explodes — `nva` reads **215** instantaneous against **7.03**
daily-mean, versus a published 6.49. Averaging is not a cosmetic adjustment
here; it is what makes those two comparable at all, and Williams' own fields
are daily means.

Daily-mean, Williams' box, DJF 2000, ours/published:

| in family (0.5–1.7) | | outlier |
|---|---|---|
| `magnitude_pv` 0.96, `brown1` 1.08, `vertical_wind_shear` 1.02, `wind_speed` 0.86, `temperature_gradient` 0.88, `endlich` 0.82, `horizontal_divergence` 0.74, `ti1` 1.08, `ti2` 1.27, `ngm1` 0.70, `ngm2` 0.89, `nva` 1.08, `rva_magnitude` 0.85, `ubf` 0.67, `brown2` 1.60, `ncsu1` 1.71, `vorticity_squared` 0.50 | | **`f2d` 22.3** |

**`deformation`'s apparent 1.93 is not an anomaly — it is the DEF² storage
seen from a third direction.** Quantiles transform monotonically, so
p97(DEF²)/p50(DEF²) = (p97/p50)². Taking the square root of the measured 5.27
gives **2.30 against a published 2.74, i.e. 0.84** — squarely in family with
the other degree-1 diagnostics. Three independent confirmations of §5 now: the
source, identity I1b, and this.

**Prediction 2 FAILED, and the failure is informative.** The shear family came
in at **1.08**, not below 1.0. The 50 hPa stencil does *not* thin the bulk
distribution's shape. That is consistent with §11.3 rather than against it:
the stencil's effect is on **tail position relative to a fixed threshold**,
which is what exceedance counting sees, not on the shape of the distribution
itself. Worth stating in the write-up — it narrows what the stencil bias is.

### 10.4 `f2d` — variant A is refuted, variant C indicated

| variant | median | p97 | p97/median | frac > 0 |
|---|---|---|---|---|
| **A** `+½ D/Dt[Q]` (current) | 7.01e-08 | 5.29e-05 | **754** | 50.9 % |
| B `−½ D/Dt[Q]` | −7.01e-08 | 5.79e-05 | −825 | 49.1 % |
| **C** `\|½ D/Dt[Q]\|` | 3.98e-06 | 9.03e-05 | **22.7** | 100 % |
| D `−D/Dt[√Q]` (literal A9) | −4.73e-07 | 8.05e-05 | −170 | 49.0 % |

against a published **13.6**. Prosser's box gives the same verdict (C at 20.5,
A at 599). §4.3's prediction — stated before running — held exactly: a signed
material derivative is centred on zero (50.9 % positive), its median is noise,
and the ratio runs to the hundreds. C is the only positive-definite reading
among the four and lands within 1.5–1.7× of the published shape, in family with
the other twenty.

Two supporting numbers:

- **A vs B tail overlap = 0.000**, exactly, in both boxes. The sign choice
  selects two disjoint populations of grid cells. It cannot be left to a
  default, which is the whole argument of §4.1.
- **Spearman ρ(A, D) = −0.94.** Negative because D carries A9's minus; the
  magnitude 0.94 says the `|∂v/∂θ|⁻¹` normalisation does reorder, but modestly.
  The sign is the large effect, the normalisation the small one.

**Not yet tested: a clipped variant**, `max(±½ D/Dt[Q], 0)`. Under daily-mean
smoothing that could also produce a finite ratio — `nva` and `ncsu1` show
exactly that behaviour. It is the obvious fifth candidate if C's trend does not
hold up.

### 10.5 `ubf` — all three mechanisms refuted

| | 1979-01 | 2020-01 | |
|---|---|---|---|
| cancellation ratio, \|residual\|/largest term | **0.729** | **0.751** | not a near-cancellation |
| Spearman ρ(float32, float64) | **1.000000** | **1.000000** | precision is a non-issue |
| flagged-set overlap at p99.6 | 1.0000 | 1.0000 | identical tails |
| \|f(ζ_ERA5 − ζ_derived)\| / \|residual\| | **0.012** | **0.011** | vorticity inconsistency is ~1 % |

The decomposition matched `2_diagnostics.ubf()` at ρ = 1.000000, so the script
was measuring the real thing.

**The most useful result is the first row, and it corrects the project's own
framing.** `STATUS.md` §4g, §7 and §12.6 all describe UBF as "a residual of
near-cancelling large terms" and "the most fragile of the 21 by construction",
and this document repeated it. **Measured, the residual is 73 % of the largest
term.** It is not a near-cancellation, nothing is amplified by the reciprocal
of a small number, and the fragility argument that has been used to justify
suspicion of UBF does not hold. Its shape ratio (0.67 daily-mean, 0.98
instantaneous) is also healthy.

Both months agree to within 3 %, so nothing here grows across the record and
none of it could produce a *trend* anomaly even if it were large.

**Conclusion: there is no longer any measured evidence that `ubf` is wrong.**
Its +11 % trend in §12.6 may simply be what ERA5 gives. It comes off the top of
the §7 list and becomes an open question rather than a suspected defect.

### 10.6 Two corrections to §4.4 and §7.2 — my own readings of Figure S5

Both panels were re-examined at 3× magnification after the numbers came back.

- **UBF.** §7.2 says S5 shows UBF "strongly positive over the NA and Europe"
  and that it "should be near the TOP of the 21". Enlarged, the North Atlantic
  patch is **positive and moderate**, comparable to the other positive panels;
  the saturated red is in the equatorial band, not the North Atlantic. The
  discrepancy against our +11 % is milder than stated — perhaps an order of
  magnitude on absolute change, not a reversal.
- **Frontogenesis.** §4.4 says S5's panel is "≈ 0 to slightly negative over the
  NA" while ours is "+226 %, 2nd of 21", and calls that a disagreement. **That
  comparison was invalid** — it puts a *relative* change against a map of
  *absolute* change, the exact error §7.2 warns about two paragraphs earlier.
  In absolute terms our f2d moves 0.026 % → 0.086 %, i.e. **+0.06 pp**, which
  on S5's ±0.5 pp scale is very nearly white. Prosser's NA frontogenesis is
  also nearly white. **They do not visibly disagree.**

Neither correction weakens the case against variant A, because that case rests
entirely on §10.4's shape measurement, which is quantitative and independent of
any figure. But it removes S5 as *corroboration*, and it means the eyeball
comparison in §7.2 should not be quoted until it is done properly — by sampling
the panels against their colourbars rather than by looking at them.

### 10.7 Revised suspect list

| diagnostic | status after measurement |
|---|---|
| `f2d` | **the only open defect.** Variant A refuted quantitatively; C indicated; trend test outstanding |
| `ubf` | **cleared.** No measured mechanism survives; the "fragile residual" framing is wrong |
| `brown2` | **closed.** Formula, magnitude and composition all confirmed |
| `deformation` | fix applied; three independent confirmations of the DEF² storage |
| the other 17 | shapes consistent with published, identities hold, formulas confirmed against Sharman |

### 10.8 What decides `f2d`, and the one thing to be careful about

The remaining test is whether variant C's DJF exceedance and trend behave like
its siblings. A concrete prediction to register before running: under A, `f2d`'s
DJF 1979 MOG exceedance is **0.026 %**, the second-lowest of the 21 and roughly
25× below the median of 0.642 % (`STATUS.md` §12.6) — a global-calibrated
threshold that a North Atlantic winter almost never crosses. **Under C it
should rise into the 0.1–1 % band where its siblings sit.** If it does not, C
is wrong too and the clipped variant is next.

**Do not run this by re-running all 21 into the existing paths.**
`jobs/04` writes `derived/global/diagnostics_glob_2000-MM.zarr` and `jobs/07`
writes `derived/north_atlantic/diagnostics_na_YYYY-MM.zarr` — the same paths
the §11 and §12 results were computed from. Overwriting them under a different
f2d variant would destroy the ability to compare, for a question that concerns
one diagnostic out of 21. A focused `f2d`-only calibrate-and-fit, writing
nothing into `derived/`, answers it at a fraction of the cost and with no risk
to the existing record.
