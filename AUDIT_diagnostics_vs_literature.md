# AUDIT — the 21 CAT diagnostics against their primary literature

**Scope.** For each of the 21 clear-air-turbulence diagnostics produced by this
repository, is the quantity the executing code computes the quantity the paper
defines?

**Date of audit:** 2026-09-08.
**Code audited:** `2_diagnostics.py` @ mtime 2026-08-30, and `rojak-cat` at the
pinned revision `25b8685c670401883bf6d186a522ccfd4561c908` (the rev in
`pixi.toml`), obtained fresh from
`github.com/ImperialCollegeLondon/rojak` and read at that commit — not from
any description of it in this repository.
**Production entry point:** `ada/diagnostics_global.py` →
`diag.compute_all_21(catdata, target_level=200, f2d_variant=<flag or
F2D_DEFAULT_VARIANT>)`. Inputs: ERA5 pressure levels **175 / 200 / 225 hPa**,
0.25°, 3-hourly, variables `u v t z d vo pv`
(`download_plan.PRESSURE_LEVELS`, `download_plan.VARIABLES`).

**Rules of evidence applied.** Nothing in this repository was treated as
evidence about the literature. `FORMULA_AUDIT.md`, `STATUS.md`, docstrings,
comments and variable names were not read as authority; where a docstring makes
a claim it is re-derived here or marked unverified. Every verdict carries an
equation number and page. Magnitude agreement was never used to settle a
formula question; where a magnitude is quoted it is quoted as a *derivation
check* (an independently predicted ratio), never as a verdict.

**Sources on disk and used:**

| Source | Used for |
|---|---|
| Sharman, Tebaldi, Wiener & Wolff (2006), *Wea. Forecasting* **21**, 268–287, Appendix A (pp. 281–284) and Appendix B Table B1 (p. 285) | reference of record for 17 of the 21 |
| Ellrod & Knapp (1992), *Wea. Forecasting* **7**, 150–165 | TI1, TI2, DEF, DST, DSH (Eqs. 1, 2, 3, 9, 10, pp. 151–152) |
| Williams (2017), *Adv. Atmos. Sci.* **34**, 576–586, Table 1 (p. 580), Table 2 (p. 580), Fig. 1 (pp. 578–579) | canonical 21-item list, canonical order, units, severity ladder, and the **sign convention for all 21** |
| Williams & Joshi (2013), *Nature Clim. Change*, Table 1 | medians, units, and the repository's `REFERENCE_TABLE` ordering |
| Koch & Caracena (2002), §2 | the nonlinear-balance-equation residual (UBF) |
| Kaplan et al. (2005), *Meteor. Atmos. Phys.* **88**, 129–153 | NCSU lineage — **does not contain the NCSU1 equation** (see Unresolved) |
| J. H. Lee et al. (2023), *JGR Atmos.*, Eqs. 1–6, p. 3 | independent ERA5 implementation of Ri, N², VWS, DEF, TI1, TI2 |
| **Kaplan et al. (2004), NASA/CR-2004-213025, Eq. 3, p. 3** — retrieved 2026-09-08 | the primary definition of NCSU1, which A36 modifies |
| **Ko, H.-C. et al. (2023), *Atmos. Chem. Phys.* **23**, 12589** — retrieved 2026-09-08 | radiosonde-vs-ERA5 evidence on how vertical resolution damages shear and N² differently |
| **Wong, H. L. et al. (2025), *JOSS* **10**(116), 9282** — on disk, previously unread | rojak's own paper: confirms it uses 175/200/225 hPa at 200 hPa, and that its methodology follows Williams & Storer (2022) |
| **Williams, P. D. & Storer, L. N. (2022), *Q. J. R. Meteorol. Soc.* **148**, 1424–1438, Eqs. (1)–(7), pp. 1427–1428; Tables 1 and 2, p. 1431** — obtained 2026-09-08 | **the only paper in Williams' own lineage that writes these diagnostics as equations.** Settles frontogenesis and NCSU1; independently reproduces the thermodynamic/dynamical split; and reports a published coarse-vs-fine resolution experiment |
| **Storer, Williams & Joshi (2017), *GRL* **44**, 9976–9984** — obtained 2026-09-08 | writes **no** equations ("the same basket of CAT diagnostics indices as Williams and Joshi (2013) and Williams (2017)"), but records that this lineage **dropped \|PV\| and used 20 diagnostics, not 21** — see §4.10 |
| **Jaeger & Sprenger (2007), *JGR* **112**, D20106, Eq. (3)** — obtained 2026-09-08 | uses only four indicators (Ri, PV, N², TI) and defines **no** Brown index; its PV is the truncated `(1/ρ)(ζ+f)∂θ/∂z`, a third independent statement of the form A18 uses |
| ~~Pearson & Sharman (2017), *JAMC* **56**, 339–351~~ — obtained 2026-09-08, **wrong paper** | this is **Part II** (nowcasting). It confirms an EDR remapping exists and cites Part I for it, but does not give it. **Part I** — Sharman & Pearson (2017), *JAMC* **56**, 317–337, "…Part I: Forecasting Nonconvective Turbulence" — is still the one needed |
| Prosser et al. (2023), *GRL* **50**, e2023GL103814 | the study being replicated: levels, method, Fig. 4 |

**Sources NOT on disk:**

* **Brown, R. (1973), *Meteor. Mag.* **102**, 347–360** — sole authority for the
  0.3 coefficient in the Brown index and for the 1/24 in the Brown EDR.
* ~~**Kaplan et al. (2004), NASA CR-2004-213025**~~ — **OBTAINED 2026-09-08**
  from NASA NTRS (`ntrs.nasa.gov/citations/20040110976`). Its Eq. (3), p. 3
  defines NCSU1 as `(U·∇U)|∇ζ|/Ri`, which differs from Sharman A36 in three
  ways. See §4.19 and §6 F9. Kaplan et al. (2005) *on disk* is the 44-case
  synoptic study and does **not** define NCSU1.
* ~~**STILL NEEDED**~~ **Williams & Storer (2022)** — **OBTAINED 2026-09-08.**
  It defines seven of the 21 as equations (Eqs. 1–7, pp. 1427–1428) and it
  closed three open items: frontogenesis (§7 Q2, now DIFFERS not ambiguous),
  the NCSU1 variant (§6 F9, withdrawn), and UBF's absolute value and β. It also
  supplied the published coarse-vs-fine resolution experiment used in §5.3. It
  does **not** define Brown, Brown EDR, |PV|, deformation, NGM1/2, |∇T|, |δ|,
  ζ², NVA, RVA or TI2, so §7 Q1 and F3 are unaffected by it.

---

## 1. Executive summary

Ranked by blast radius. One sentence each; full derivations below.

1. **`spatial_gradient`'s meridional derivative carries a spurious factor
   `a/M(φ)` — a latitude-structured error of +0.42 % at 30 °N falling to −0.27 %
   at 75 °N — because `nominal_grid_spacing` returns the *true geodesic*
   meridional arc `M Δφ` while `spatial_gradient` then also multiplies by the
   meridional map-scale factor `a/M`, double-correcting; the zonal derivative is
   exact.** Reaches every horizontal derivative in the system: **13 of the 21**
   (measured 2026-09-08 by running both versions — the 8 it leaves
   bit-identical are exactly the 8 with no horizontal ∂/∂y: `−Ri`,
   `vertical_wind_shear`, `colson_panofsky`, `endlich`, `wind_speed`, and the
   three that take an ERA5 archived field whole, `magnitude_pv`,
   `vorticity_squared`, `horizontal_divergence`). (§3, B1/B2; measured, not
   asserted, in §3.1.)
2. **The vertical stencil is 175/225 hPa (Δp = 50 hPa) where Prosser uses
   188/206 hPa (Δp = 18 hPa), and the damage is confined to the diagnostics that
   differentiate a *thermodynamic or directional* field vertically — not to the
   whole stencil group.** Prosser (2023, p. 5) names the stencil as the reason
   his results differ from Lee et al. (2023). Phase 3 sharpens it: diagnostics
   that differentiate only `u` and `v` vertically replicate as well as ones that
   differentiate nothing (median level ratio 0.873 versus 0.934), while the six
   that differentiate `θ`, `T` or wind direction sit at a median of **0.502**.
   The mechanism is that a stencil change acts as a near-uniform *rescaling* of
   each individual derivative — measured here as 1.45x on `Sv²` against 1.06x on
   `N²` — and a rescaling cancels exactly out of a percentile calibration. So it
   is harmless for diagnostics that are pure PRODUCTS of derivatives and decisive
   for those that combine them as a RATIO or a DIFFERENCE, which partitions the
   ten stencil diagnostics exactly. Corroborated against radiosondes by Ko et al.
   (2023), who find ERA5 underestimates shear by ~50 % while its N² is good to
   ~11 %. This is a *data* difference, not a formula defect, and it is not
   fixable: 175/200/225 is the finest stencil the CDS pressure-level product
   offers — and it is the stencil rojak's own paper uses (Wong et al. 2025).
   (§3.3, §5.3, §6 F7.)
3. **Brown index (#5 Williams order) is, by construction, ≈ √0.3 · f — a
   latitude field with a small flow-dependent perturbation** — and Brown EDR
   (#6) inherits it. A13 is transcribed exactly, and because Prosser computes it
   the same way the f-dominance cancels: both replicate cleanly (level ratios
   1.087 and 0.913). Recorded not as a defect but because it means these two of
   the 21 carry very little flow information. (§4.5, §5.2.)
4. **Frontogenesis (#4) computes `|½ D/Dt Q|`, and the published definition is
   the signed `D/Dt Q`.** Williams & Storer (2022) Eq. (3), p. 1427 — obtained
   after the first draft of this audit, and the only equation for this
   diagnostic anywhere in Williams' lineage — has no absolute value, no leading
   minus and no normalisation. The project's variant **A** is right and the
   current default **C** is wrong. This was recorded as AMBIGUOUS-IN-SOURCE and
   is now a defect with a citation. Note the unresolved tension it leaves inside
   Williams' own papers. (§4.4, §7 Q2, §6 F6.)
5. **NCSU1's variant question is closed, in favour of the code.** Kaplan et al.
   (2004) Eq. 3 — the origin, retrieved 2026-09-08 — defines NCSU1 with the full
   advection `U·∇U`, a plain Ri and no clipping; Sharman A36 adds a
   half-contraction, a `max(Ri,10⁻⁵)` floor and an outer `max(·,0)`. **Williams
   & Storer (2022) Eq. (7) reproduces A36 exactly**, so A36 is what the
   replicated lineage uses and what the code correctly implements. NCSU1's 6×
   disagreement is therefore *not* a variant mismatch, and that candidate is
   eliminated. (§4.19, §5.6, §6 F9 withdrawn.)
6. **Gradient diagnostics are unbounded at the poles**: PROJ's
   `parallel_scale` for `+proj=latlon` is 2.28 × 10² at 89.75 °N and
   9.97 × 10⁴ at 90 °N, so ∂/∂x is inflated by those factors in the two polar
   rows of the *global calibration* domain. cos φ weighting suppresses but does
   not remove them; the severe (p99.9) thresholds are the exposed ones. (§3.2.)
7. **This project takes ζ, δ and PV from ERA5's archive; Prosser computed all
   three from `u, v, T, z`** — his Methods list (2023, p. 2) contains no
   vorticity, divergence or potential vorticity. Ten of the 21 are affected. The
   axis does not explain the ensemble, but it is decisive for the two
   worst-replicating diagnostics that have *no* vertical derivative: **UBF**
   (archived ζ against computed `J` and `∇²Φ`, inside a residual of
   near-cancelling terms — the worst possible place for two operators) and
   **|PV|** (ERA5's full Ertel PV where Sharman A18 defines the vertical term
   only). Both were flagged from the literature in Phase 2 and both are flagged
   empirically. (§3.4, §4.10, §4.17, §5.5.)
8. **Colson–Panofsky's λ is the ICAO-standard-atmosphere layer thickness
   (816.7 m), constant in space and time, while its Sv² and N² are centred over
   the full 50 hPa stencil** — a factor-of-two stencil mismatch against
   Sharman's "local value of vertical grid increment". Because λ is constant the
   whole thing is a multiplicative constant and *exactly* cancels out of the
   exceedance field. MATCHES A4 algebraically; recorded for the record. (§4.3.)

**No defect was found in the algebra of:** −Ri, VWS, CP, Brown1, Brown2, TI1,
TI2, DEF, ζ², |∇ₕT|, |v|, Endlich, NGM1, NGM2, UBF, |δ|, NCSU1, NVA. Their
verdicts are MATCHES against a cited equation. **|PV| is DIFFERS (defensibly)**
and **RVA is UNVERIFIABLE** — no primary equation for it exists in any source on
disk. The substantive open questions are frontogenesis's absolute value, Brown's
0.3 (unverifiable), NCSU1's variant (unverifiable), and the shared-block items
above.

**What Phase 3 changed.** Nothing in §4 was revised after the empirical
comparison was opened. Two exec-summary items were re-ranked: Brown's
f-dominance was demoted (both studies share it, so it cancels), and the vertical
stencil was promoted and sharpened once the group A/B/C split in §5.3 showed the
damage is confined to thermodynamic and directional vertical derivatives. Of the
eight diagnostics whose exceedance level disagrees with Prosser, **all eight are
accounted for by two mechanisms Phase 2 identified from the literature alone**
— the stencil (six) and archived-versus-computed ERA5 fields (two) — with the
single exception of `ncsu1`'s *sign*, which remains open (§5.6).

---

## 2. The dependency graph (Phase 0)

Derived from the code, not assumed. Numbering is **Williams (2017) Table 2
order** throughout this report; the repository's `REFERENCE_TABLE["num"]` is
Williams & Joshi (2013) Table 1 order (ranked by % change) and is *different* —
see §4.0.

### 2.1 Which code path executes

`compute_all_21` runs 14 diagnostics through `rojak.DiagnosticFactory` and 7
through local functions, then applies two post-hoc corrections. Verified by
reading the dispatch, not the comments:

| # (W2017) | key | executing path |
|---|---|---|
| 1 | `negative_richardson` | **hand** `richardson(ds, negative=True)` |
| 2 | `vertical_wind_shear` | rojak `VerticalWindShear` |
| 3 | `colson_panofsky` | **hand** `colson_panofsky(ds)` |
| 4 | `f2d` | **hand** `frontogenesis_isentropic(..., variant=C)` |
| 5 | `brown1` | rojak `BrownIndex1` |
| 6 | `brown2` | **hand** `brown2(ds, brown1)` where `brown1` is a *fresh* rojak `BrownIndex1` |
| 7 | `ti1` | rojak `TurbulenceIndex1` |
| 8 | `ti2` | rojak `TurbulenceIndex2` |
| 9 | `deformation` | rojak `DeformationSquared`, **then `sqrt(abs(·))` applied in `compute_all_21`** |
| 10 | `magnitude_pv` | rojak `MagnitudePotentialVorticity` |
| 11 | `vorticity_squared` | rojak `VerticalVorticitySquared` |
| 12 | `temperature_gradient` | rojak `HorizontalTemperatureGradient` |
| 13 | `wind_speed` | rojak `WindSpeed` |
| 14 | `endlich` | rojak `Endlich` |
| 15 | `ngm1` | rojak `NestedGridModel1` |
| 16 | `ngm2` | rojak `NestedGridModel2` |
| 17 | `ubf` | **hand** `ubf(ds, target_level=200)` |
| 18 | `horizontal_divergence` | rojak `HorizontalDivergence` |
| 19 | `ncsu1` | **hand** `ncsu1(ds, target_level=200)` |
| 20 | `nva` | rojak `NegativeVorticityAdvection` |
| 21 | `rva_magnitude` | **hand** `rva(ds, target_level=200)` |

Two of the hand-written diagnostics are hand code *wrapping rojak primitives*:
`brown2` calls `DiagnosticFactory(...).create(BROWN1).computed_value` and then
multiplies by its own `Sv²`; `ncsu1` and `ubf` call
`rojak.core.derivatives.vector_derivatives`. `f2d`'s isentropic derivative uses
a local `theta_derivative_on_pressure_level` modelled on rojak's
`altitude_derivative_on_pressure_level`.

`frontogenesis_2d()` (the Miller / constant-pressure form) is **dead code** — it
is defined in `2_diagnostics.py` but the `builders` dict routes `"f2d"` to
`frontogenesis_isentropic`. It is audited below only to record that it is not
the executing path.

### 2.2 Shared building blocks and their blast radii

| Block | Definition in code | Consumed by (W2017 #) | Blast radius if wrong |
|---|---|---|---|
| **B1 — scalar horizontal gradient** `spatial_gradient(f, DEG, GEOSPATIAL)` = `np.gradient(f, cumsum(nominal Δ))` × PROJ scale factor | `derivatives.py:305` | 4 (∇Q), 12 (∇T), 17 (∇Φ), 19 (∇ζ), 20 (∇ζₐ), 21 (∇ζ) | 6 direct. Latitude-structured if the metric term is wrong. |
| **B2 — velocity-gradient tensor** `vector_derivatives(u,v,DEG)`, metpy convention | `derivatives.py:384`, cached on `CATData` | 5, 6, 7, 8, 9, 15, 16 via DEF; 17 via J; 19 via `du_dx, dv_dy` | **9 of 21.** Largest single block. |
| **B3 — vertical derivative** `g (∂f/∂p)(∂Φ/∂p)⁻¹` | `calculations.py:122` | 1, 2, 3, 6, 7, 8, 14, 16, 19 (via Ri), 4 (θ-analogue) | **10 of 21.** |
| **B4 — Sv² = (∂u/∂z)² + (∂v/∂z)²** (built on B3) | `calculations.py:317` | 1, 2, 3, 6, 7, 8, 19 | 7 of 21. |
| **B5 — N² = (g/θ)(∂θ/∂z)** (built on B3, B7) | `diagnostic.py:660` | 1, 3, 19 | 3. |
| **B6 — Ri = N²/Sv²** (B4 ÷ B5) | `richardson()` | **1, 19** | 2 direct. CP (#3) does *not* route through Ri in the executing code — it uses N² and Sv² directly, algebraically identically (§4.3). |
| **B7 — θ = T(1000/p)^κ**, κ = 0.285714… | `calculations.py:209` | 1, 3, 4, 19 (via N²) | 4. |
| **B8 — DEF = hypot(D_sh, D_st)** (built on B2) | `data.py:204` | 5, 6, 7, 8, 9, 15, 16 | 7 of 21. |
| **B9 — ERA5 archived ζ** (`vo`) | `data.py:146` | 5, 6, 11, 17, 19, 20, 21 | 7 of 21. |
| **B10 — ERA5 archived δ** (`d`) | `data.py:128` | 8, 18 | 2. |
| **B11 — ERA5 archived PV** (`pv`) | `data.py:143` | 10 | 1. |
| **B12 — f = 2Ω sin φ**, Ω = 7.292115 × 10⁻⁵ s⁻¹ | `calculations.py:231` | 5, 6, 17, 20 | 4. **Latitude-structured by construction.** |
| **B12b — β = 2Ω cos φ / R_E** | hand `ubf()` only | 17 | 1. |
| **B13 — ICAO altitude coordinate** | `_icao_altitude()` | 3 (λ only) | 1, and only as a constant. |
| **B14 — Eulerian ∂/∂t**, 3-hourly | `frontogenesis_isentropic` | 4 | 1. Only diagnostic with a time derivative. |

**Direction of the effect, per block:**

* B1/B2 too large in ∂/∂y at low latitude and too small at high latitude ⇒ every
  gradient-derived diagnostic has a weak *equator-ward* bias in the 30–75 °N
  band, monotone in φ, ≤ 0.7 % peak-to-peak.
* B3/B4 with a 50 hPa stencil instead of 18 hPa ⇒ vertical derivatives are
  **smoothed**: extremes reduced, distribution narrowed, upper tail thinner. All
  10 consumers biased the same way, in the same direction.
* B12 ⇒ Brown1, Brown2 and (weakly) NVA inherit a monotone-in-|φ| field. In a
  50–75 °N or 30–60 °N box measured against a *global* threshold, this alone
  sets the exceedance level.

---

## 3. Phase 1 — the shared building blocks, audited

### 3.1 Horizontal metric terms — the cos φ question, settled by measurement

`spatial_gradient(array, units, GradientMode.GEOSPATIAL)` does two things
(`derivatives.py:305–348`):

1. `grid_deltas = nominal_grid_spacing(lat, lon, units)`, which computes
   `dx` as the **geodesic distance along the equator** between consecutive
   longitudes — i.e. `a Δλ`, the plate-carrée map distance — and `dy` as the
   **geodesic distance along a meridian** between consecutive latitudes — i.e.
   the true meridional arc `M(φ) Δφ`, **not** `a Δφ`.
2. `first_derivative` = `np.gradient(f, cumsum(Δ), axis)`, then multiplies by
   `parallel_scale` (for x) or `meridional_scale` (for y) from
   `pyproj.Proj(CRS("+proj=latlon")).get_factors`.

Measured directly (pyproj 3.x, WGS84):

```
lat      parallel_scale     a/(N cos φ)      meridional_scale     a/M(φ)
30       1.15373388         1.15373388       1.00421324           1.00421
45       1.41184476         1.41184476       1.00168911           1.00169
60       1.99497290         1.99497290       0.99916709           0.99917
75       3.85161817         3.85161817       0.99732219           0.99732
```

So `parallel_scale = a/(N cos φ)` and `meridional_scale = a/M(φ)` exactly.

**x-direction — correct.** Raw derivative = `(1/(a Δλ)) ∂f/∂λ`; the correct
value is `(1/(N cos φ)) ∂f/∂λ`; the required factor is `a/(N cos φ)` =
`parallel_scale`. The 1/cos φ metric term **is present and is exact on the
ellipsoid.**

**y-direction — double-corrected.** Raw derivative is already
`(1/(M Δφ)) ∂f/∂φ = ∂f/∂y`, because `dy` is the *true* arc. Multiplying by
`meridional_scale = a/M` therefore introduces a spurious factor `a/M(φ)`.
The two directions are not in the same coordinate system: `dx` is a map
distance (`a Δλ`) and `dy` is a ground distance (`M Δφ`), so only `dx`
composes correctly with its scale factor.

**Verified numerically** by running rojak's own `spatial_gradient` at this rev
on two analytic fields with exactly-known gradients on the WGS84 ellipsoid
(30–80 °N × 60–0 °W, 0.25°):

```
field f = eastward arc along the parallel = N(φ) cos φ · (λ − λ₀); true ∂f/∂x = 1
  lat 30: 1.000000   err +0.0000 %
  lat 45: 1.000000   err +0.0000 %
  lat 60: 1.000000   err +0.0000 %
  lat 75: 1.000000   err +0.0000 %

field f = meridional arc = ∫M dφ; true ∂f/∂y = 1
  lat 30: 1.004213   err +0.4213 %
  lat 45: 1.001689   err +0.1689 %
  lat 60: 0.999167   err −0.0833 %
  lat 75: 0.997322   err −0.2678 %
```

The x error is machine zero. The y error is exactly `a/M(φ) − 1`, monotone
decreasing in latitude, spanning **0.69 % across 30–75 °N**.

This is a **latitude-structured error**, of exactly the class §6 asks about,
and it is invisible in a global median. It is also small: 0.7 % peak-to-peak.
`divergence()` is unaffected in the sense that the *sum* `DU_DX + DV_DY` was
verified exact (machine zero) for an analytic non-divergent solid-body flow,
because the `DV_DY` half is the only one carrying the factor and for that flow
it is identically zero; for a general flow the `∂v/∂y` half carries it.

**Fix (proposed, not applied):** in `nominal_grid_spacing`, `dy` should be the
map distance `a Δφ` (equatorial radius × angle), matching `dx`'s construction,
so that `meridional_scale = a/M` converts it correctly. Equivalently, leave
`dy` as the geodesic and do not multiply by `meridional_scale`. Either makes
the y-direction exact. This is upstream code; the local mitigation is a
monkey-patch or a local `spatial_gradient` wrapper.

### 3.2 Pole behaviour

`parallel_scale` at |φ| → 90° is **not** clamped to anything physical:

```
lat        88        89       89.75        90
p_scale  28.56     57.11     228.42    99664.72
```

Every `∂/∂x` in the global calibration domain (`GLOBAL_BOX = [90,−180,−90,180]`)
is therefore multiplied by up to 10⁵ in the polar rows. `calibration.py` applies
cos φ weighting (`lat_weights`, `weighted_percentile`), which gives the two
polar rows a combined weight fraction of order 2 × 10⁻⁵ — negligible against the
3 % light tail, but only ~50× below the 0.1 % severe tail. **Recorded as a risk
to the severe thresholds specifically, not as a demonstrated defect.**

There is **no** longitude-wrap defect: after the −180…180 conversion and
`sortby`, `np.gradient` uses one-sided differences at the two edge columns.
These are first-order rather than second-order accurate but they are not wrong —
no false gradient is manufactured across the date line.

### 3.3 Vertical discretisation and the hydrostatic conversion

`altitude_derivative_on_pressure_level(f, Φ)` returns
`g · (∂f/∂p) / (∂Φ/∂p)` (`calculations.py:122–148`). Since `Φ = gz`,
`∂f/∂z = g ∂f/∂Φ = g (∂f/∂p)(∂Φ/∂p)⁻¹`. **Correct**, and the pressure units
cancel identically, so the hPa-vs-Pa question does not arise. `xarray
.differentiate` on a 3-level uniform coordinate gives a centred difference at
200 hPa over **175 → 225 hPa** and one-sided at 175 and 225 (which are
discarded).

Against the replication target: **Prosser et al. (2023), p. 2** — "fields on the
188 and 206 hPa levels were also extracted" — a 18 hPa centred stencil at
197 hPa. This project uses a **50 hPa** stencil at 200 hPa. Prosser (p. 5)
attributes the disagreement between his results and Lee et al. (2023)
specifically to this: "we used input fields at 206 hPa and 188 hPa … In
contrast, J. H. Lee et al. (2023) appear to have used input fields at 200 hPa
and 300 hPa … This means we have used much finer (and therefore more accurate)
vertical finite differences," and notes that the two diagnostics *without* a
vertical derivative (deformation and divergence) are the two that agree between
the studies.

This is the single most consequential input difference in the replication and it
is not a formula defect. It affects B3, B4, B5, B6, B7 and therefore
diagnostics 1, 2, 3, 4, 6, 7, 8, 14, 16, 19 — and it does **not** cancel out of a
percentile calibration, because coarsening a derivative changes the *shape* of
the distribution (thinner tail), not merely its scale.

The **ICAO altitude coordinate** (`_icao_altitude`) is used for exactly one
thing: Colson–Panofsky's λ. It is a function of pressure alone, so it is
constant in space and time (816.7 m for 175→200 hPa; 737.8 m for 200→225 hPa).
Everywhere else the vertical coordinate comes from ERA5's own geopotential.

### 3.4 Archived versus computed ERA5 fields

`CATData` returns **archived** ERA5 fields for divergence (`d`), vorticity
(`vo`) and potential vorticity (`pv`) (`data.py:128, 143, 146`); everything else
is computed from `u, v, t, z`. Both are defensible. The choice is not recorded
anywhere in the repository, and three expressions **mix** the two:

* **TI2** = `Sv (DEF − δ)` uses a locally computed `DEF` and ERA5's archived `δ`.
  ERA5's divergence is derived from the spectral representation of the wind
  field, so `δ_archived ≠ ∂u/∂x + ∂v/∂y` as computed here. The two enter the
  same bracket.
* **UBF** = `|∇²Φ − 2J(u,v) − fζ + βu|` uses ERA5's archived `ζ` for the `fζ`
  term but locally computed derivatives for `J`. UBF is a near-cancellation of
  large terms, so a mixed-provenance vorticity is the *worst* place to have one.
  (The hand implementation's own docstring makes this point; it is re-derived
  here independently and it is correct.)
* **Brown1** = `√(0.3 ζₐ² + D_sh² + D_st²)` uses archived `ζ` and computed
  deformation.

There is a further inconsistency worth naming: **ERA5's archived `pv` is the
full Ertel PV** (it includes the horizontal-vorticity × horizontal-θ-gradient
terms), whereas **Sharman A18/A19 define PV = −g ζₐ ∂θ/∂p** (the vertical term
only). Lee et al. (2023) Eq. 6, p. 3 writes the same truncated form,
`PV = (1/ρ)(ζ+f) ∂θ/∂z`, which is hydrostatically identical to A18. So the
implemented |PV| is a *different, more complete* quantity than the published
definition. Defensible-but-different.

### 3.5 Constants

| Constant | Value in code | Authority | Verdict |
|---|---|---|---|
| g | 9.80665 m s⁻² | standard | ok |
| Ω | 7.292115 × 10⁻⁵ s⁻¹ | NIST, cited in `calculations.py:24` | ok |
| R_E | 6 371 008.7714 m | mean Earth radius | ok; note the *derivatives* use WGS84 `a` = 6 378 137 via pyproj, so β (which uses R_E) and ∂/∂y (which uses `a`) rest on different radii. 0.1 % — recorded, not material. |
| κ = R/c_p | 0.28571428571428564 (= 2/7) | Wallace & Hobbs Eq. 3.54 | ok |
| p₀ | 1000 hPa | Poisson | ok |
| Ri_crit | 0.5 | Sharman A4, p. 281: "Ri_crit is an empirical constant (≈ 0.5)" | ok |
| 0.3 | Brown1 | Sharman A13, p. 283, attributed to Brown (1973) | **unverifiable at source** — see §7 Q1 |
| 1/24 | Brown2 | Sharman A14, p. 283 (and A11 SCATR) | **unverifiable at source** — see §7 Q1 |
| 10⁻⁵ | NCSU1 Ri floor | Sharman A36, p. 284: `MAX(Ri, 10⁻⁵)` | ok, verbatim |
| ½ | frontogenesis | *not* in A9; introduced by the `½ D/Dt Q` reading. A constant, so immaterial to ranks. | ok as a choice |

The hand `frontogenesis_2d()` (dead code) uses `κ = 287/1004 = 0.285857`, which
differs from `2/7` in the fifth decimal. It does not execute; noted only so the
discrepancy is not rediscovered.

---

## 4. Phase 2 — the 21 audits

Written before `data_prosser/per_diagnostic_vs_prosser.csv`, `STATUS.md` §17 or
any list of suspect diagnostics was opened.

### 4.0 A note on ordering, and on the sign convention for all 21

The repository's `REFERENCE_TABLE["num"]` is **not** Williams (2017) Table 2
order — it is Williams & Joshi (2013) Table 1 order, which is *ranked by
percentage change in the median* and therefore has no canonical status. Prosser
(2023) Figure 4's panel order is Williams (2017) Table 2's row order. This
report uses Williams Table 2 order. The `REFERENCE_TABLE` units and W&J medians
were checked line-by-line against W&J Table 1 and are **all correct**, including
the two that differ between the two papers (`rva_magnitude` and `nva`: W&J gives
10⁻¹⁰ s⁻², Williams 2017 gives 10⁻⁹ s⁻²; the table follows W&J, which is where
the medians come from — consistent).

**Question 5 — which tail indicates turbulence, derived from the source.**
Williams (2017) Table 2 (p. 580) gives, for each of the 21, the onset thresholds
for Light → Light-to-moderate → Moderate → Moderate-to-severe → Severe. In
**every one of the 21 rows** the five values are **monotonically increasing**,
including the two negative-valued rows:

* Negative Richardson number: −15.4 → −9.8 → −7.9 → −6.7 → **−5.9**
* Colson–Panofsky index: −29.3 → −27.0 → −25.2 → −23.7 → **−22.2**

A monotonically increasing ladder is a one-tailed **upper**-tail criterion by
construction, because Table 1 (p. 580) assigns those categories to the
97.0 / 99.1 / 99.6 / 99.8 / 99.9 percentiles — i.e. progressively further into
the *right*-hand tail. Sharman's Table B1 (p. 285) independently confirms it for
the two signed ones: −Ri runs −20 → −2.0 → −0.6 → −0.3 → 0.5 and CP runs
0 → 1000 → 5000 → 12 000 → 30 000 kt², both increasing.

**Verdict: for all 21, the turbulent tail is the upper one.** The
`REFERENCE_TABLE` `"sign": "+"` entries are therefore correct for all 21,
including `colson_panofsky` and `negative_richardson`. This is an independent
re-derivation, from Williams Table 1 + Table 2 and Sharman Table B1, of the two
entries the project record says were previously wrong. It agrees with the
current values — but it was derived from the papers, not read off the table.

---

### 4.1 (#1) Negative Richardson number

| field | value |
|---|---|
| **source** | Sharman (2006) **A1**, p. 281, with **A2** and **A3**, p. 281. Independently: Lee et al. (2023) **Eq. 5**, p. 3, and **Williams & Storer (2022) Eq. (1), p. 1427** — the last being Williams' own lineage, i.e. the one Prosser inherits. W&S write it in pressure coordinates, `Ri = −(1/ρθ)(∂θ/∂p)/(∂**u**/∂p)²`, which is **algebraically identical** to A1 under hydrostatic balance: substituting `∂p = −ρg ∂z` gives `N² = −(ρg²/θ)(∂θ/∂p)` and `Sv² = (ρg)²(∂**u**/∂p)²`, whose ratio is exactly W&S Eq. (1). **Four independent statements of the same formula, and the code matches all four.** |
| **published** | `Ri = N²/Sv²`; `N² = (g/θ)(∂θ/∂z)` (A2); `Sv = |∂**v**/∂z| = [(∂u/∂z)² + (∂v/∂z)²]^{1/2}` (A3). Diagnostic is `−Ri`. |
| **implemented** | `richardson()`: `n2 = (g/θ)·d(θ)/dz`; `sv_squared = (du/dz)² + (dv/dz)²`; `ri = n2/sv_squared`; returns `−ri`. |
| **variant** | A2 offers `g/θ ∂θ/∂z` **or** `g/θₑ ∂θₑ/∂z` (equivalent-potential-temperature form). The dry form is implemented. Williams & Joshi / Prosser give no indication which; the dry form is the GTG default and the only one computable from the seven variables downloaded (no `q` is used — it is stubbed with zeros). Recorded as a variant choice with no evidence either way. |
| **constants** | g = 9.80665; κ = 2/7. |
| **sign** | Upper tail. Derived: Sharman p. 281 — "regions of small Ri should be favored regions of turbulence" ⇒ small Ri ⇒ large −Ri; and Table B1's −Ri ladder is increasing. |
| **clipping** | None. Ri is unbounded where Sv² → 0; no floor, no mask. A2/A3 specify none. |
| **inputs** | u, v, t, z at 175/200/225; computed. |
| **units** | dimensionless ✓ (Williams Table 2: "—"). |
| **verdict** | **MATCHES** Sharman (2006) A1 with A2, A3, p. 281 — including the square on Sv, which is where rojak's `GradientRichardson._compute` (`diagnostic.py:754`, `ri = brunt_vaisala / vws`) departs. That rojak defect was **re-confirmed at rev 25b8685** by reading the source, and it is not on the executing path. |

Note: rojak's own docstring for `GradientRichardson` states `Ri = N²/Sv²`
correctly while the code divides by `Sv`. This is exactly why docstrings were
excluded as evidence.

---

### 4.2 (#2) Magnitude of vertical shear of horizontal wind

| field | value |
|---|---|
| **source** | Sharman **A3**, p. 281. Lee et al. (2023) Eq. 1 (VWS term), p. 3. |
| **published** | `Sv = [(∂u/∂z)² + (∂v/∂z)²]^{1/2}` |
| **implemented** | rojak `VerticalWindShear` → `vertical_wind_shear(u, v, geopotential=z)` → `np.hypot(du_dz, dv_dz)` with `du_dz = g(∂u/∂p)(∂Φ/∂p)⁻¹`. `is_abs_velocities=False`, `is_vws_squared=False`. |
| **variant** | A3 as printed writes `|∂v/∂z|` on the left with `(|∂u/∂z|² + |∂v/∂z|²)^{1/2}` on the right; the absolute values inside are redundant under the square. Ellrod & Knapp (1992, p. 152, Eq. 6) instead use a **bulk layer shear** `ΔV/Δz` over a 50–100 mb layer. The pointwise derivative is implemented, which is Sharman's and Lee's form. |
| **constants** | g only. |
| **sign** | Upper tail (Williams Table 2: 5.3 → 8.4, increasing). |
| **clipping** | None. |
| **inputs** | u, v, z; computed. Vertical stencil 175↔225 (§3.3). |
| **units** | s⁻¹ ✓ (Williams: 10⁻³ s⁻¹). |
| **verdict** | **MATCHES** Sharman A3, p. 281. |

---

### 4.3 (#3) Colson–Panofsky index

| field | value |
|---|---|
| **source** | Sharman **A4**, p. 281 (from Colson & Panofsky 1965, not on disk). |
| **published** | `CP = λ² Sv² (1 − Ri/Ri_crit)`, λ = "a length scale, taken as the local value of vertical grid increment Δz", `Ri_crit ≈ 0.5`. |
| **implemented** | `colson_panofsky()`: `cp = λ² (Sv² − N²/Ri_crit)`, with λ = `altitude.diff("pressure_level", label="upper")` and Sv², N² computed on the full three-level stencil and then selected onto λ's levels. |
| **algebra** | `λ²Sv²(1 − Ri/Ri_crit) = λ²Sv² − λ²Sv²·(N²/Sv²)/Ri_crit = λ²(Sv² − N²/Ri_crit)`. **Identical** to A4, given A1. Not an approximation. |
| **variant** | Sharman also cites Laikhtman & Al'ter-Zalik's A5 with a different length scale; not implemented, correctly. **Williams & Storer (2022) Eq. (2), p. 1427** restates A4 verbatim with `Ri_crit = 0.5` and defines Δz as "the vertical grid spacing" — confirming that λ is the grid increment and therefore, on a fixed pressure-level set, a constant. |
| **constants** | Ri_crit = 0.5 ✓ A4, p. 281. |
| **sign** | Upper tail — see §4.0. CP is almost everywhere negative at 200 hPa because `N²/0.5 ≫ Sv²`; Williams (2017, p. 580) says exactly this: "the Colson–Panofsky thresholds in Table 2 are negative, because the Colson–Panofsky index is proportional to 1−Ri/0.5, and the Richardson number (Ri) is rarely less than 0.5". |
| **clipping** | None. A4 specifies none. |
| **inputs** | u, v, t, z + the ICAO altitude coordinate. |
| **units** | Implemented in **m² s⁻²**, tagged `"units": "m2 s-2"`. Sharman Table B1 (p. 285) tabulates **kt²**; Williams Table 2 tabulates **10³ kt²**. **No conversion is applied.** 1 m² s⁻² = 3.7785 kt², so `CP[10³ kt²] = CP[m² s⁻²] × 3.7785 × 10⁻³`. |
| **verdict** | **MATCHES** Sharman A4, p. 281, algebraically and in constants. Two recorded discrepancies, neither affecting the exceedance field: (i) units are native SI, un-converted; (ii) λ is a *constant*. |

**Two stencil observations.** λ is the 175→200 ICAO thickness = **816.7 m**
(or 737.8 m for 200→225, depending on whether cfgrib returns the level
coordinate ascending or descending — the code works either way but picks a
different λ). Sv² and N² are *centred* over 175→225, i.e. an effective layer
twice as thick as λ. Sharman's "local value of vertical grid increment" is
ambiguous between the two; the mismatch is recorded, not resolved. Because λ is
a function of pressure alone it is **constant in space and time**, so `CP` is
`const × (Sv² − 2N²)` and every percentile, rank and exceedance is *exactly*
that of `(Sv² − 2N²)`.

**Derivation check (not a verdict).** Inverting Williams' light threshold:
−29.3 × 10³ kt² = −7755 m² s⁻²; with typical lower-stratospheric N² ≈ 4 × 10⁻⁴ s⁻²
and Sv² ≈ 4 × 10⁻⁶ s⁻², A4 gives λ ≈ 3113 m — a plausible GFDL-CM2.1 layer
thickness at 200 hPa. Substituting this project's λ = 816.7 m into the same
expression gives −2.01 × 10³ kt². The predicted ratio is `(3113/816.7)² = 14.5`,
i.e. the entire magnitude offset against Williams is a λ² difference and nothing
else. This is offered as a consistency check on the *reading* of A4; it is not
evidence about the code, and the ratio cancels from every percentile.

---

### 4.4 (#4) Frontogenesis function

| field | value |
|---|---|
| **source** | Sharman **A9**, p. 282 (isentropic, "the form used in GTG at upper levels"), with **A10**, p. 282 (the constant-pressure midlevel form) and **Table B1**, p. 285. |
| **published (A9, verbatim)** | `F_θ ∝ − D/Dt[(∂u/∂θ)² + (∂v/∂θ)²]^{1/2} = |∂**v**/∂θ|⁻¹ [ (∂u/∂θ)·D/Dt(∂u/∂θ) + (∂v/∂θ)·D/Dt(∂v/∂θ) ]` |
| **implemented** | `frontogenesis_isentropic(..., variant="C")`: `Q = (∂u/∂θ)² + (∂v/∂θ)²`; `D/Dt = ∂/∂t + u ∂/∂x + v ∂/∂y`; result `= |½ D/Dt Q|`. |
| **the printed inconsistency** | Write `Q = (∂u/∂θ)² + (∂v/∂θ)²`. The bracket on A9's right-hand side is `½ D/Dt(Q)`, and `|∂**v**/∂θ| = Q^{1/2}` (bold **v** is the horizontal wind vector, exactly as in A3). So A9's RHS `= ½ D/Dt(Q) / Q^{1/2} = D/Dt(Q^{1/2})`, while A9's LHS is `−D/Dt(Q^{1/2})`. **The two sides of A9 as printed differ by a sign.** This is a defect in the source, not in the code, and it is the origin of the whole question. |
| **the algebra, settled by A10 + Table B1** | Table B1 (p. 285) gives F_θ (A9) units **m² s⁻³ K⁻²**. `D/Dt(Q^{1/2})` has units (m s⁻¹ K⁻¹)/s = **m s⁻² K⁻¹** — *not* Table B1's. `D/Dt(Q)` has units (m s⁻¹ K⁻¹)²/s = **m² s⁻³ K⁻²** — exactly Table B1's. The confirmation is A10: the midlevel form `F_p ∝ (∂u/∂p)D/Dt(∂u/∂p) + (∂v/∂p)D/Dt(∂v/∂p)` **is** `½ D/Dt(Q_p)`, un-normalised, and Table B1 lists it as **m² s⁻³ Pa⁻²** — the exact p-analogue. So GTG's F_θ is the un-normalised `½ D/Dt(Q)` in θ, and A9's printed `|∂**v**/∂θ|⁻¹` is spurious. Williams (2017) Table 2 lists **10⁻⁹ m² s⁻³ K⁻²** and W&J Table 1 lists **10⁻⁹ m² s⁻³ K⁻²**, both agreeing with Table B1 and therefore with the un-normalised form. **The algebra is settled: variants A/B/C are right and D is wrong.** |
| **the sign / absolute value — NOT settled** | Nothing in Sharman states an absolute value. A9's leading minus puts frontogenesis at *decreasing* Q, which is physically right on an isentropic surface (strengthening \|∇θ\| ⇒ weakening \|∂**v**/∂θ\|). Table B1's thresholds for F_θ are all **positive and increasing** (10⁻⁵ → 7.1 × 10⁻⁴), which tells us GTG uses a positive-tail quantity but not whether it got there by a minus sign, a clip, or an absolute value. Williams (2017) Fig. 1 plots this diagnostic on **0…300 × 10⁻⁹**, anchored at zero — in the same figure where −Ri is plotted on −300…0 and CP on −45…−25, so Williams does show signed diagnostics on negative ranges and frontogenesis is not one of them. That establishes non-negativity; it does not distinguish `\|X\|` from `max(X,0)`. |
| **variant** | Four readings are enumerated in `F2D_VARIANTS`; **C = `\|½ D/Dt Q\|`** is `F2D_DEFAULT_VARIANT` and is what `ada/diagnostics_global.py` passes unless overridden. Prosser (2023) does not state which reading he used; he inherits "the 21 diagnostics used in Williams and Joshi (2013) and Williams (2017)" (p. 2) and neither paper writes the equation. |
| **constants** | ½ (a choice, immaterial to ranks); κ = 2/7 for θ. |
| **sign** | Upper tail (Williams Table 2: 770 → 2340, increasing). |
| **clipping** | `np.abs`. **Not specified by any source.** |
| **inputs** | u, v, t, z. `∂/∂θ` via chain rule `(∂/∂p)/(∂θ/∂p)` on the existing 175/200/225 levels — *not* an isentropic remap. `∂/∂t` via `xarray.differentiate("time", datetime_unit="s")`, centred over ±3 h. |
| **units** | m² s⁻³ K⁻² ✓ Table B1, Williams Table 2, W&J Table 1. |
| **RESOLVED 2026-09-08 — Williams & Storer (2022) Eq. (3), p. 1427** | The replicated lineage does write this one down, and it reads, in full: `F_θ = D/Dt \|∂**u**/∂θ\|²` — "where t is time, D/Dt denotes the Lagrangian time derivative, and the partial derivative is taken at fixed horizontal position using potential temperature as a vertical coordinate". Since `**u** = (u, v)`, `\|∂**u**/∂θ\|² = Q`. So the published quantity is **`D/Dt[Q]`: signed, un-normalised, with NO leading minus, NO absolute value and NO clip.** Up to the constant ½ that is immaterial to ranks, that is **variant A** — the project's pre-2026-08-30 default. |
| **verdict** | **DIFFERS.** The executing variant is **C**, `\|½ D/Dt Q\|`. The absolute value has no authority in Sharman (2006), Williams (2017), Williams & Joshi (2013) **or Williams & Storer (2022)**, and the last of those states the formula explicitly. The algebra (un-normalised `D/Dt Q`) was already established from A10 + Table B1 and is confirmed. **The sign treatment is not ambiguous any more: it is wrong.** See §6 F6. |

**Two consequences of choosing `|·|` that should be on the record.**

1. It **conflates frontogenesis with frontolysis**: a strongly *weakening* front
   scores identically to a strongly strengthening one. The reading most faithful
   to A9's printed physics is `max(−½ D/Dt Q, 0)`, which would keep only
   frontogenetic cells. Williams' Fig. 1 first bin (≈ 18 %) is the project's
   stated evidence against a clip, on the grounds that a clip would put ~50 % of
   the mass in the first bin. That argument is sound **only if** `½ D/Dt Q` is
   symmetric about zero, which is an assumption about the atmosphere, not a
   reading of the paper. It is a good argument; it is not a citation.
2. The **discrete material derivative is part of the definition in practice.**
   `∂/∂t` is centred over 6 h and the advection term uses instantaneous `u,v`.
   For a feature advected at 50 m s⁻¹, a parcel traverses ~1080 km — about 40
   grid points — in that interval. The two terms of `D/Dt` are individually
   large and nearly cancelling, so the result is a small residual of two big
   numbers evaluated on a coarse stencil. No source specifies the sampling
   interval. **Recorded as an unresolved discretisation choice**, and as the
   reason F2D is the only one of the 21 whose value depends on the time
   chunking (`ada/diagnostics_global.py` handles this with a one-step overlap
   buffer; the mechanism is correct, though its docstring's claim that
   `differentiate` "drops the first and last timestep" is wrong — it returns a
   one-sided estimate there, not NaN).

**Not the executing path:** `frontogenesis_2d()` implements the Bluestein
constant-pressure form `F = −|∇θ|⁻¹[θx² u_x + θy² v_y + θx θy (v_x + u_y)]`.
Against **rojak's** `Frontogenesis2D._compute` (`diagnostic.py:291–296`) which
computes `θx² u_x + θy² v_y + θx θy v_y + θx θy u_y` — the third term's
coefficient is `dv_dy` where the docstring's own LaTeX says `∂v/∂x`. **That
rojak defect is confirmed at rev 25b8685** and is not on the executing path.
Neither form is the right diagnostic for #4: Table B1 gives F_θ in m² s⁻³ K⁻²,
and the constant-pressure kinematic form has units K² m⁻² s⁻¹, which cannot
produce it.

---

### 4.5 (#5) Brown index

| field | value |
|---|---|
| **source** | Sharman **A13**, p. 283, attributed to Brown (1973), *Meteor. Mag.* **102**, 347–360 — **not on disk**. |
| **published** | `Φ = (0.3 ζₐ² + D_SH² + D_ST²)^{1/2}`, with `D_SH = ∂v/∂x + ∂u/∂y`, `D_ST = ∂u/∂x − ∂v/∂y`, `ζₐ = ζ + f`, `ζ = ∂v/∂x − ∂u/∂y`, f the Coriolis frequency (Sharman p. 283, defining all four in the same paragraph). |
| **implemented** | rojak `BrownIndex1._compute`: `sqrt(0.3·ζₐ² + D_sh² + D_st²)` with `ζₐ = ζ_ERA5 + 2Ω sin φ`, `D_sh = dv_dx + du_dy`, `D_st = du_dx − dv_dy` from `vector_derivatives`. |
| **variant** | None known. |
| **constants** | **0.3** — see §7 Q1. |
| **sign** | Upper tail (Williams: 99 → 118). |
| **clipping** | None; the square root of a sum of squares is non-negative by construction. |
| **inputs** | u, v computed (deformation); ζ archived; f from latitude. |
| **units** | s⁻¹ ✓ (Williams: 10⁻⁶ s⁻¹). |
| **verdict** | **MATCHES** Sharman A13, p. 283, transcribed exactly, including the deformation definitions which independently match Ellrod & Knapp (1992) Eqs. 1 and 2, p. 152. |

**The finding that matters here is not a coding error.** At 200 hPa,
`0.3 ζₐ² ≈ 0.3 f²` dominates the deformation terms by an order of magnitude:
at 60 °N, `0.3 f² = 4.9 × 10⁻⁹` against `D² ≈ (2 × 10⁻⁵)² = 4 × 10⁻¹⁰`.
So `Brown1 ≈ √0.3 · f` to within ~5 %:

```
lat      √0.3·f (10⁻⁶ s⁻¹)
30       39.9
45       56.5
50       61.2
62       70.5
75       77.2
```

W&J's pre-industrial median over 50–75 °N is **77.1 × 10⁻⁶ s⁻¹** — consistent
with `√0.3 · f` at the poleward end of the box plus a deformation contribution.
**Brown1 is, to first order, a map of latitude.** This is a property of the
published diagnostic, not of the implementation, but it has a direct
consequence for this project's method: a diagnostic whose value is set by `f`
has a *region-versus-globe contrast* determined entirely by the box's latitude
band, so its exceedance frequency against a global threshold is near-degenerate.
Carried to §6.

---

### 4.6 (#6) Brown energy dissipation rate

| field | value |
|---|---|
| **source** | Sharman **A14**, p. 283: `ε = (1/24) Φ Sv²`, with Φ from A13. Compare **A11** (SCATR), p. 282: `SCATR = (1/24) Φ Sv²` with Φ from **A12**. |
| **published** | `ε_brown = (1/24) Φ Sv²`. |
| **implemented** | `brown2(ds, brown1)`: `(1/24) · brown1 · sv_squared`, where `brown1` is a freshly-created rojak `BrownIndex1` and `sv_squared = (du/dz)² + (dv/dz)²` computed locally on the full three-level stencil. Identical in form to rojak's `BrownIndex2._compute` (`(1/24)·brown1·square(vws)`). |
| **variant** | rojak's *docstring* states a piecewise form — `ε = (1/24)ΦSv²` if `dRi/dt > 0`, else `0` — attributed to Brown (1973). **Neither rojak's code nor this repository's implements the piecewise clip**, and Sharman A14 as printed has no condition. Since Φ = A13 is a square root and therefore always ≥ 0, the condition `dRi/dt > 0` is never binding *if* `dRi/dt` is identified with A13's Φ; the clip is vacuous under Sharman's A13. Recorded because rojak's docstring asserts a source (Brown 1973) that cannot be checked. |
| **constants** | 1/24 — see §7 Q1. |
| **sign** | Upper tail (Williams: 870 → 2330). |
| **clipping** | None. See "variant". |
| **inputs** | Brown1 (archived ζ + computed DEF) × locally computed Sv². |
| **units** | **s⁻³** as implemented (Φ in s⁻¹ × Sv² in s⁻²). Williams Table 2 and W&J Table 1 both give **10⁻⁶ J kg⁻¹ s⁻¹ = 10⁻⁶ m² s⁻³**. A14 as printed cannot produce m² s⁻³ — it yields s⁻³. There is an **unstated length² factor** in the published units. The code keeps native s⁻³ and does not invent one. |
| **checked 2026-09-08, still open** | Three further papers were searched for the missing length²: **Storer et al. (2017)** writes no equations; **Jaeger & Sprenger (2007)** uses four indicators and none of them is Brown's; **Pearson & Sharman (2017)** is Part II and only cites Part I for the EDR remapping. The gap is unchanged. The two candidates that could still close it are **Sharman & Pearson (2017) Part I**, *JAMC* **56**, 317–337 (where GTG's diagnostic-to-EDR remapping is described) and **Brown (1973)** itself. |
| **verdict** | **MATCHES** Sharman A14, p. 283, as printed. The units disagreement with Williams/W&J is **AMBIGUOUS-IN-SOURCE**: A14 is dimensionally s⁻³ and the published tables are m² s⁻³, and nothing on disk supplies the missing L². Since an unknown *constant* L² cancels from every percentile, this is inert for the exceedance field and fatal only for magnitude comparison — which is how the code treats it. Correct handling; the ambiguity is real and belongs in the record. |

---

### 4.7 (#7) Ellrod turbulence index, variant 1

| field | value |
|---|---|
| **source** | **Ellrod & Knapp (1992) Eq. 9, p. 152**: `TI1 = VWS × DEF`; **Sharman A15**, p. 283: `TI1 = Sv DEF`; independently **Lee et al. (2023) Eq. 1, p. 3**. |
| **published** | `TI1 = Sv · DEF`, `DEF = (D_ST² + D_SH²)^{1/2}` (E&K Eq. 3, p. 152 = Sharman A17, p. 283). |
| **implemented** | rojak `TurbulenceIndex1._compute`: `vertical_wind_shear(u,v,geopotential) * CATData.total_deformation()`, where `total_deformation` is `magnitude_of_vector(D_sh, D_st, is_squared=False)` = `hypot`. **Un-squared** — confirmed by reading `data.py:204`. |
| **variant** | E&K compute VWS as a bulk layer shear over 300–400 mb (p. 153); the pointwise A3 derivative is used here, consistent with Sharman and Lee. |
| **constants** | none. |
| **sign** | Upper tail (195 → 472). |
| **clipping** | None. |
| **inputs** | u, v, z; all computed. |
| **units** | s⁻² ✓ (Table B1 p. 285: s⁻²; Williams: 10⁻⁹ s⁻²). |
| **verdict** | **MATCHES** Ellrod & Knapp (1992) Eq. 9, p. 152 and Sharman A15/A17, p. 283. |

---

### 4.8 (#8) Ellrod turbulence index, variant 2

| field | value |
|---|---|
| **source** | **Ellrod & Knapp (1992) Eq. 10, p. 152**: `TI2 = VWS × [DEF + CVG]`, with **CVG = convergence = −(∂u/∂x + ∂v/∂y)** defined on p. 152 immediately under Eq. 4. **Sharman A16**, p. 283: `TI2 = Sv(DEF − Δ_H)` with `Δ_H = ∂u/∂x + ∂v/∂y` (A33, p. 284). **Lee et al. (2023) Eq. 2, p. 3**: `TI2 = VWS × (DEF − DIV)`. |
| **published** | `TI2 = Sv (DEF − δ)`. All three sources agree. |
| **implemented** | rojak `TurbulenceIndex2._compute`: `convergence = -self._divergence; return vws * (total_deformation + convergence)` = `Sv (DEF − δ)`. |
| **variant** | E&K note TI2 is the AFGWC variant and TI1 the NMC variant (p. 152). Both are computed here as separate diagnostics, correctly. |
| **constants** | none. |
| **sign** | Upper tail (184 → 477). **The divergence sign is the load-bearing one and it is right**: turbulence is favoured by *convergence*, i.e. δ < 0, which makes `−δ > 0` and raises TI2. Derived from E&K's own words on p. 152 ("The TI used at AFGWC retains the CVG term") plus their definition CVG = −(u_x + v_y). |
| **clipping** | None. |
| **inputs** | **Mixed**: DEF computed from `vector_derivatives`; δ is **ERA5's archived `d`** (`data.py:128`). §3.4. |
| **units** | s⁻² ✓. |
| **verdict** | **MATCHES** Ellrod & Knapp (1992) Eq. 10, p. 152; Sharman A16, p. 283. Mixed field provenance recorded as a defensible-but-unrecorded choice. |

---

### 4.9 (#9) Flow deformation

| field | value |
|---|---|
| **source** | Sharman **A17**, p. 283; **Ellrod & Knapp (1992) Eqs. 1–3, pp. 151–152**. |
| **published** | `DEF = (D_SH² + D_ST²)^{1/2}`; `D_ST = ∂u/∂x − ∂v/∂y` (E&K Eq. 1); `D_SH = ∂v/∂x + ∂u/∂y` (E&K Eq. 2). |
| **implemented** | rojak's `DEF` diagnostic is `DeformationSquared` — it returns `square(total_deformation)`, i.e. **DEF²**. `compute_all_21` then applies `np.sqrt(np.abs(·))`, restoring `DEF`. Verified in the dispatch (`diagnostic.py:1378`) and in `2_diagnostics.py:999–1014`. |
| **variant** | none. |
| **constants** | none. |
| **sign** | Upper tail (50.9 → 76.3). |
| **clipping** | `np.abs` before the square root. On a non-negative input this is a float-noise guard only; it cannot change a value. |
| **inputs** | u, v; computed via `vector_derivatives`. |
| **units** | s⁻¹ ✓ (Williams: 10⁻⁶ s⁻¹). Note that **DEF² would have units s⁻²**, so the units tag alone distinguishes the two, and Williams' 10⁻⁶ s⁻¹ settles it in favour of DEF. |
| **verdict** | **MATCHES** Sharman A17, p. 283, *on the production path only*. `compute_rojak_diagnostics()` alone returns DEF²; anything that calls it directly (e.g. `tests/`) gets the squared quantity. Recorded as a live footgun, not a defect in the output. |

---

### 4.10 (#10) Magnitude of potential vorticity

| field | value |
|---|---|
| **source** | Sharman **A18**, p. 283: `|PV|`, with `PV = −g ζₐ ∂θ/∂p` (p. 283). Lee et al. (2023) Eq. 6, p. 3: `PV = (1/ρ)(ζ+f) ∂θ/∂z` — hydrostatically identical. |
| **published** | `|PV|` with PV the *vertical-term-only* Ertel PV. |
| **implemented** | rojak `MagnitudePotentialVorticity` → `np.abs(CATData.potential_vorticity())` → `np.abs(ds["potential_vorticity"])`, i.e. **ERA5's archived `pv`**. rojak *has* a `potential_vorticity()` function implementing A18 exactly (`calculations.py:281`), and it is **not used** by this diagnostic. |
| **variant** | **This is a real variant difference.** ERA5's archived PV is the full Ertel PV `−g η·∇θ`, including the horizontal-vorticity × horizontal-θ-gradient contributions that A18 drops. At 200 hPa in the jet the dropped terms are small but not zero and are largest exactly where the shear is largest — i.e. in the diagnostic's own upper tail. **Defensible-but-different.** |
| **constants** | none. |
| **sign** | Upper tail (8.33 → 9.41 PVU). |
| **clipping** | `np.abs`, specified by A18's `|PV|`. |
| **inputs** | archived `pv`. |
| **units** | ERA5 `pv` is K m² kg⁻¹ s⁻¹ (SI); Williams and W&J tabulate **PVU = 10⁻⁶ SI**. No conversion applied; a constant factor, inert for percentiles. |
| **verdict** | **DIFFERS (defensibly)** from Sharman A18, p. 283: the implemented quantity is the full Ertel PV, not `−g ζₐ ∂θ/∂p`. Not a bug; an undocumented substitution. A18's truncated form is stated three times independently — Sharman p. 283, Lee et al. (2023) Eq. 6 p. 3, and Jaeger & Sprenger (2007) Eq. (3) — so there is no doubt what the published diagnostic is. |

**A finding that reframes this one rather than fixing it.** Storer, Williams &
Joshi (2017), p. 2 states: "we calculate the same basket of CAT diagnostics
indices as Williams and Joshi (2013) and Williams (2017), **except that we
exclude the potential vorticity diagnostic because it was found to give
unrealistic results**" — and they work with **20** diagnostics, not 21. So
Williams' own group, in the paper sitting between W&J (2013) and Williams
(2017), judged this diagnostic unreliable and dropped it. Prosser (2023) put it
back. That is worth knowing before F3 is run: `magnitude_pv` is the one
diagnostic in the set whose own originators removed it, it is the one that costs
this replication its 17th significant trend (t = 1.97 against 2.021), and a
large flip under F3 is at least as likely to be telling you something about the
diagnostic as about this pipeline. **F3 still stands** — A18 is what Prosser
must have computed, and matching his definition is the point — but the result
should be read with Storer's sentence beside it.

---

### 4.11 (#11) Relative vorticity squared

| field | value |
|---|---|
| **source** | Sharman **A21**, p. 283: `ζ² = \|∇ × **v**\|²`, motivated on p. 283 as a curvature measure with maxima in troughs and minima in ridges. |
| **published** | `ζ²`. |
| **implemented** | rojak `VerticalVorticitySquared` → `np.square(ds["vorticity"])` — **ERA5's archived `vo`**. |
| **variant** | Archived vs computed. ERA5's `vo` is the true spherical relative vorticity from the spectral wind representation; the locally computed `dv_dx − du_dy` from `vector_derivatives` would be the alternative and is *not* used here. Archived is the better choice. |
| **constants** | none. |
| **sign** | Upper tail (2.46 → 6.24). Note A21 squares precisely so that **both** troughs (ζ > 0) and ridges (ζ < 0) map to the upper tail — Sharman p. 283 states this intent explicitly. |
| **clipping** | none (the square is the operator). |
| **inputs** | archived `vo`. |
| **units** | s⁻² ✓ (Williams: 10⁻⁹ s⁻²). |
| **verdict** | **MATCHES** Sharman A21, p. 283. |

---

### 4.12 (#12) Magnitude of horizontal temperature gradient

| field | value |
|---|---|
| **source** | Sharman **A23**, p. 283. |
| **published** | `\|∇_H T\| = [(∂T/∂x)² + (∂T/∂y)²]^{1/2}` |
| **implemented** | rojak `HorizontalTemperatureGradient` → `magnitude_of_geospatial_gradient(temperature)` → `hypot(dfdx, dfdy)` from `spatial_gradient(..., GEOSPATIAL)`. |
| **variant** | A23 is on temperature, not potential temperature. Temperature is used ✓. |
| **constants** | none. |
| **sign** | Upper tail (14.7 → 22.0). |
| **clipping** | none. |
| **inputs** | `t` at 200 hPa; computed gradient. |
| **units** | K m⁻¹ ✓ (Table B1 p. 285: K m⁻¹; Williams: 10⁻⁶ K m⁻¹). |
| **verdict** | **MATCHES** Sharman A23, p. 283 — subject to the shared **B1** meridional-metric defect of §3.1, which inflates the `∂T/∂y` component by `a/M(φ)`. Since `\|∇T\|` at 200 hPa is dominated by the meridional component in the jet, this diagnostic is among the most exposed to that ≤ 0.7 % latitude-structured error. |

---

### 4.13 (#13) Wind speed

| field | value |
|---|---|
| **source** | Sharman **A24**, p. 283: `s = \|**v**\|`. |
| **implemented** | rojak `WindSpeed` → `np.hypot(u, v)`. |
| **verdict** | **MATCHES** Sharman A24, p. 283. Units m s⁻¹ ✓ (Table B1 p. 285 and Williams both m s⁻¹). Upper tail (40.9 → 58.5). No clipping, no constants, no derivatives — the only one of the 21 that consumes no shared building block at all. |

---

### 4.14 (#14) Wind speed × directional shear (Endlich)

| field | value |
|---|---|
| **source** | Sharman **A25**, p. 283, after Endlich (1964): `s \|∂ψ/∂z\|`, ψ the wind direction. |
| **published** | `\|**v**\| · \|∂ψ/∂z\|`. |
| **implemented** | rojak `Endlich._compute`: `wind_direction(u,v)` (meteorological convention, `π/2 − atan2(−v,−u)` wrapped to [0,2π]); vertical gradient of ψ via `angles_gradient` with a `__sub__` override that takes the **shorter** of `\|Δψ\|` and `2π − \|Δψ\|`; then `g·(∂ψ/∂p)/(∂Φ/∂p)`; then `speed × abs(·)`. |
| **the wrap-around, verified** | The `_WrapAroundAngleArray` override only takes effect if `np.gradient` actually calls `__sub__`. NumPy reduces a uniformly-spaced coordinate array to a scalar and then uses `(f[i+1] − f[i−1])/(2Δ)` — which *does* call `__sub__`. The project's level coordinate 175/200/225 is uniform (Δp = 25), so the override **is** active. On a **non-uniform** vertical coordinate NumPy would instead use the three-point weighted form `a·f + b·f + c·f`, which contains no subtraction, and the wrap-around would silently stop working. **Recorded as a latent fragility, not a current defect.** |
| **variant** | The override returns `\|Δψ\|` ≥ 0, so the *sign* of ∂ψ/∂z is destroyed. Immaterial: A25 takes the absolute value. |
| **constants** | g. |
| **sign** | Upper tail (3.21 → 5.08). |
| **clipping** | `np.abs`, specified by A25. |
| **inputs** | u, v, z. |
| **units** | rad s⁻¹ × m s⁻¹. Williams and W&J both label it **10⁻³ rad s⁻¹** — which is dimensionally inconsistent with `s·∂ψ/∂z` (that would be m rad s⁻² ). The published unit label is wrong in both papers; the *quantity* is unambiguous from A25. Recorded. |
| **verdict** | **MATCHES** Sharman A25, p. 283, and **Williams & Storer (2022) Eq. (5), p. 1427** (`√(u²+v²)·\|∂ψ/∂z\|`), verbatim. |

**But its empirical disagreement is now unexplained.** W&S Table 1, p. 1431 puts "wind speed × directional shear" among the *most robust* of their seven to a large coarsening of both the horizontal grid and the vertical levels — inter-dataset standard deviation 0.30–0.40 %, second-best of the seven. This project's `endlich` sits at a level ratio of **0.51** against Prosser, worse than its stencil peers. **So Endlich's disagreement is probably NOT the vertical stencil**, and §5.3's group C should not have contained it. The leading new suspect is the only part of this diagnostic that is not a plain centred difference: rojak's `angles_gradient` and its `_WrapAroundAngleArray.__sub__`, which returns the *shorter* of `\|Δψ\|` and `2π − \|Δψ\|` and therefore discards the sign. W&S compute ψ's derivative with "centred second-order finite differences" (p. 1428) and say nothing about wrap-around handling. **New open item; see §6 F10.**

---

### 4.15 (#15) Flow deformation × wind speed (NGM1)

| field | value |
|---|---|
| **source** | Sharman **A28**, p. 284, from Reap (1996) — not on disk. |
| **published** | `\|**v**\| DEF`. |
| **implemented** | rojak `NestedGridModel1` → `wind_speed(u,v) * CATData.total_deformation()` (un-squared DEF ✓). |
| **sign** | Upper tail (1.65 → 3.54). |
| **units** | m s⁻² ✓ (Table B1 p. 285: m s⁻²; Williams: 10⁻³ m s⁻²). |
| **verdict** | **MATCHES** Sharman A28, p. 284. |

---

### 4.16 (#16) Flow deformation × vertical temperature gradient (NGM2)

| field | value |
|---|---|
| **source** | Sharman **A29**, p. 284, from Reap (1996). |
| **published** | `\|∂T/∂z\| DEF`. |
| **implemented** | rojak `NestedGridModel2` → `abs(altitude_derivative_on_pressure_level(T, Φ)) * total_deformation()`. |
| **sign** | Upper tail (53 → 151). |
| **clipping** | `np.abs`, specified by A29's `\|·\|`. |
| **units** | K m⁻¹ s⁻¹ ✓ (Williams: 10⁻⁹ K m⁻¹ s⁻¹). |
| **verdict** | **MATCHES** Sharman A29, p. 284. |

---

### 4.17 (#17) Magnitude of residual of the nonlinear balance equation (UBF)

| field | value |
|---|---|
| **source** | Sharman **A30**, p. 284; **Koch & Caracena (2002) §2**, which writes `R = −∇²Φ + 2J(u,v) + fζ − βu ≠ 0` — character-for-character the same expression. |
| **published** | `UBF = −∇²Φ + 2J(u,v) + fζ − βu`, β the Coriolis frequency gradient (Sharman p. 284). Neither Sharman nor Koch & Caracena prints an absolute value; Williams (2017) Table 2 and W&J Table 1 name the diagnostic "**Magnitude** of residual of nonlinear balance equation", and **Williams & Storer (2022) Eq. (6), p. 1427 prints it explicitly: `UBF = \|−∇²Φ + 2J(u,v) + fζ − βu\|`, with β stated as "the latitudinal derivative of the Coriolis parameter"** — confirming both the absolute value and that β is `∂f/∂y`, i.e. `2Ω cos φ / R`, which is what the hand implementation uses and what rojak's `latitudinal_derivative` gets wrong. |
| **implemented** | hand `ubf()`: `residual = ∇²Φ − 2J − fζ + βu`, then `np.abs`. This is `−(A30)`; **under the absolute value it is identical to `\|A30\|`.** |
| **∇²Φ** | Computed as `div(grad Φ)`: `grad Φ` via the scalar `spatial_gradient` (correct — Φ is a scalar), then its divergence via `vector_derivatives(dΦ/dx, dΦ/dy)`. Derivation: `div(A,B) = (a cos φ)⁻¹∂A/∂λ + a⁻¹∂B/∂φ − B tan φ / a`; with `A = Φ_x, B = Φ_y` this is `Φ_xx + Φ_yy − (tan φ/a)Φ_y`, which **is** the spherical Laplacian. Correct. rojak's `spatial_laplacian` at rev 25b8685 does the same thing (`derivatives.py:365–373`), so the two agree — the pin note's claim that upstream fixed this is confirmed by reading the source. |
| **β** | `2Ω cos φ / R_E` ✓ the true `∂f/∂y`. rojak's `latitudinal_derivative` (`calculations.py:249–262`) returns `coriolis_param / R_E = 2Ω sin φ / R_E` while its own docstring states `2Ω cos φ / R_E`. **That rojak defect is confirmed at rev 25b8685** and is not on the executing path. |
| **J(u,v)** | `u_x v_y − u_y v_x` from `vector_derivatives`. **AMBIGUOUS:** A30 is written in Cartesian notation and does not say how the Jacobian is to be evaluated on a sphere. metpy's convention (which rojak follows) distributes the curvature term `u tan φ/a` asymmetrically into `DU_DX` and `DV_DX`, so that `div` and `curl` come out right but the individual partials are not the "true" spherical partials. At 60 °N, `u tan φ/a ≈ 8 × 10⁻⁶ s⁻¹` against a typical `u_x ≈ 1 × 10⁻⁵ s⁻¹` — **comparable, not negligible**, and latitude-structured. A different (equally defensible) distribution would give a materially different J. Recorded as an implementation choice with no authority in the sources. |
| **variant** | A30's `−∇²Φ` versus the code's `+∇²Φ`: immaterial under `\|·\|`, since the whole expression is negated. rojak's `UBF._compute` (`diagnostic.py:657`) computes `\|mass + inertial − β u\|` with `mass = +∇²Φ`, which is **not** `\|A30\|` — it flips the sign of one term only. Confirmed at rev 25b8685; not on the executing path. |
| **constants** | Ω, R_E. |
| **sign** | Upper tail (1230 → 2960). |
| **clipping** | `np.abs`, authorised by Williams (2017) Table 2's name, not by A30. |
| **inputs** | **Mixed**: archived ζ; computed J and ∇²Φ; u archived-none. §3.4. |
| **units** | s⁻² ✓ (Table B1 p. 285: s⁻²; Williams: 10⁻¹² s⁻²). |
| **verdict** | **MATCHES** Sharman A30, p. 284 and Koch & Caracena (2002) §2, up to an overall sign that the absolute value removes. Two recorded ambiguities: the spherical Jacobian convention, and the mixed provenance of ζ inside a residual. |

---

### 4.18 (#18) Magnitude of horizontal divergence

| field | value |
|---|---|
| **source** | Sharman **A33**, p. 284: `Δ_H = ∂u/∂x + ∂v/∂y`, listed in §n "Unbalanced flow" as a measure of unbalanced flow. |
| **published** | `Δ_H`. Williams Table 2 and W&J both take the **magnitude**. |
| **implemented** | rojak `HorizontalDivergence` → `np.abs(ds["divergence_of_wind"])` = `\|ERA5 archived d\|`. |
| **variant** | Archived vs computed. Archived is arguably the better choice (spectral, consistent with ERA5's own dynamics) and it is also the *only* one of the four gradient-derived quantities in this set that is exact in the §3.1 test. |
| **sign** | Upper tail (11.9 → 22.5). |
| **clipping** | `np.abs`, authorised by the published name, not by A33. |
| **units** | s⁻¹ ✓ (Williams: 10⁻⁶ s⁻¹). |
| **verdict** | **MATCHES** Sharman A33, p. 284, with `\|·\|` from Williams (2017) Table 2. |

---

### 4.19 (#19) NCSU index, version 1

| field | value |
|---|---|
| **source** | Sharman **A36**, p. 284. Sharman attributes it to **Kaplan et al. (2004), NASA/CR-2004-213025**, which was **retrieved from NASA NTRS on 2026-09-08** and is now a source of record for this audit. Kaplan et al. (2005) *is* on disk and does **not** contain this equation; it is the 44-case synoptic study, and the only NCSU index it mentions (p. 134, item 23) is "NCSU modification of the Ellrod index (Ellrod index ÷ ipv)" — **a different index**. Any project record citing "Kaplan (2005)" for the NCSU1 equation is citing the wrong paper. |
| **published (A36)** | `NCSU1 = [1/MAX(Ri, 10⁻⁵)] · MAX(u ∂u/∂x + v ∂v/∂y, 0) · \|∇ζ\|` |
| **implemented** | hand `ncsu1()`: `advection = (u·du_dx + v·dv_dy).clip(min=0)`; `grad_zeta = hypot(dζ/dx, dζ/dy)` from the scalar `spatial_gradient`; `ri_safe = ri.clip(min=1e-5)`; `out = advection · grad_zeta / ri_safe`. `du_dx, dv_dy` come from `vector_derivatives` (projection-corrected), **not** from the scalar gradient — correct, since u and v are vector components. |
| **variant** | **SETTLED 2026-09-08, and not in favour of the alternative.** Three sources are now in hand and they do not agree with each other. **Kaplan et al. (2004) Eq. (3), p. 3** — the origin — reads `NCSU1 = (U·∇U)·|∇ζ| / Ri`: full advection, plain Ri, no floor, no clip. **Sharman A36, p. 284** introduces all three departures. **Williams & Storer (2022) Eq. (7), p. 1428** reproduces Sharman A36 *exactly*, including both `max()` operators: `NCSU1 = [1/max(Ri,10⁻⁵)]·max(u ∂u/∂x + v ∂v/∂y, 0)·|∇ζ|`. Since W&S is Williams' own lineage — the lineage Prosser inherits — **A36 is definitively the right form for this replication, and the implemented code matches it.** The Kaplan form, though it is the primary definition and differs materially (note that `u u_x + v v_y` is not any natural scalar reduction of `U·∇U`, whose components are `(u u_x + v u_y, u v_x + v v_y)` — it takes one term from each and drops the cross terms), is **not** what is being replicated. The proposed F9 test is therefore withdrawn: NCSU1's 6× disagreement is not a variant mismatch. |
| **constants** | 10⁻⁵ ✓ verbatim from A36. |
| **sign** | Upper tail (1200 → 13 000). |
| **clipping** | **Both clips are specified by the source.** `MAX(u u_x + v v_y, 0)` is printed in A36; `MAX(Ri, 10⁻⁵)` is printed in A36. This answers Question 3 for NCSU1: not an implementation choice. |
| **inputs** | u, v computed derivatives; ζ archived; Ri from B6 (and hence from the 50 hPa vertical stencil). |
| **units** | s⁻³ ✓ (Table B1 p. 285: s⁻³; Williams: 10⁻¹⁸ s⁻³). |
| **verdict** | **MATCHES** Sharman A36, p. 284, verbatim — including the floor's asymmetry. |

**Flagged consequence, faithful to the source.** `MAX(Ri, 10⁻⁵)` is a *lower*
bound, not `\|Ri\|`. Wherever `N² < 0` — convective instability, which ERA5 does
resolve at 200 hPa in deep-convective columns — Ri is negative, the floor
returns 10⁻⁵, and `1/Ri_safe = 10⁵`. Every such cell is amplified by five orders
of magnitude and lands in the extreme tail. This is what A36 says. It means
NCSU1's upper tail is, in part, a map of where the reanalysis has negative
static stability on a 50 hPa stencil. It is not a bug and should not be "fixed";
it is a property of the published diagnostic that the vertical-stencil choice of
§3.3 interacts with directly, since a coarser stencil changes *how often* N² < 0.

---

### 4.20 (#20) Negative absolute vorticity advection

| field | value |
|---|---|
| **source** | Sharman **A37**, p. 284, after Bluestein (1992, p. 335). |
| **published** | `NVA = MAX{ [ −u ∂(ζ+f)/∂x − v ∂(ζ+f)/∂y ] , 0 }` |
| **implemented** | rojak `NegativeVorticityAdvection._compute`: `abs_vorticity = ζ_ERA5 + 2Ω sin φ`; `x = u · ∂ζₐ/∂x`; `y = v · ∂ζₐ/∂y`; `nva = −x − y`; `nva.clip(min=0)`. |
| **variant** | The gradient is taken of the **full** `ζ + f` field, so the `β` contribution enters through `∂ζₐ/∂y` automatically rather than being added separately. Correct and equivalent. |
| **constants** | Ω. |
| **sign** | Upper tail (1.33 → 2.93). |
| **clipping** | `MAX(·, 0)` — **specified by A37**, p. 284. Answers Question 3 for NVA. |
| **inputs** | archived ζ, computed gradient, f from latitude. |
| **units** | s⁻² ✓ (Williams: 10⁻⁹ s⁻²; W&J: 10⁻¹⁰ s⁻²). |
| **verdict** | **MATCHES** Sharman A37, p. 284. |

---

### 4.21 (#21) Magnitude of relative vorticity advection

| field | value |
|---|---|
| **source** | **No primary equation exists on disk.** Sharman (2006) Appendix A contains A37 (negative *absolute* vorticity advection) but no relative-vorticity-advection entry. Williams & Joshi (2013) Table 1 lists "Magnitude of relative vorticity advection" and refers only to "their usual definitions¹,¹⁴⁻¹⁶,²⁰"; Williams (2017) Table 2 repeats the name; Prosser (2023) inherits the list wholesale (p. 2). |
| **published** | Reconstructible from the name only: `\|**v** · ∇ζ\|`, ζ **relative** (the name says "relative", distinguishing it from #20 which says "absolute"). |
| **implemented** | hand `rva()`: `np.abs(u · ∂ζ/∂x + v · ∂ζ/∂y)` with ζ = ERA5 archived `vo` and the scalar geospatial gradient. |
| **variant** | The sign of the advection is immaterial under `\|·\|`, so `\|**v**·∇ζ\|` and `\|−**v**·∇ζ\|` are the same field. Nothing else is in doubt. |
| **constants** | none. |
| **sign** | Upper tail (1.44 → 3.00). |
| **clipping** | `np.abs`, **authorised by the published name** ("Magnitude of…") and by nothing else. **This answers Question 3 for RVA: it is neither Sharman-specified nor an arbitrary implementation choice — the only authority is Williams' and W&J's title for the diagnostic, and that authority does specify a magnitude rather than a one-sided clip.** Note the contrast with #20, where the source *does* specify a one-sided `MAX(·,0)`; the two neighbouring diagnostics genuinely use different operators, and the code correctly uses different operators. |
| **inputs** | archived ζ, computed gradient. |
| **units** | s⁻² ✓. |
| **verdict** | **UNVERIFIABLE** against a primary equation — none exists in any source on disk. The implementation is the only reading consistent with the published name and with its neighbour A37, and its units match. Recorded as unverifiable rather than as matching, because there is no equation to match against. |

---

## 5. Phase 3 — reconciliation with the empirical evidence

`data_prosser/per_diagnostic_vs_prosser.csv` and `STATUS.md` §17.5 were opened
only after §1–§4 above were written and saved (snapshot kept, timestamped, in
the session scratchpad). Prosser's Figure 4 prints `rel`, `abs` and `p` per
panel; his 1979 *levels* are read off the position of his red crosses and are
good to ~±0.02 — so `level_ratio` is a sign-and-magnitude quantity, not an exact
one, and a ratio of 0.87 versus 0.93 means nothing. A ratio of 0.14 or 6.07
means a great deal.

### 5.1 The whole table, sorted by level ratio

```
colson_panofsky        0.138    trend ratio    1.37
ubf                    0.193                   0.63
f2d                    0.425                  87.4
negative_richardson    0.494                  −5.41
endlich                0.511                   1.36
magnitude_pv           0.613                   0.51
ngm2                   0.615                   0.83
--------------------------------------------------- everything below is fine
vertical_wind_shear    0.850                   1.02
ti2                    0.868                   1.06
horizontal_divergence  0.874                   0.77
ti1                    0.878                   1.04
rva_magnitude          0.901                   1.00
temperature_gradient   0.913                   1.04
brown2                 0.913                   1.00
ngm1                   0.934                   1.03
wind_speed             0.943                   0.98
nva                    0.973                   1.02
deformation            0.999                   0.99
brown1                 1.087                   0.98
vorticity_squared      1.121                   1.07
ncsu1                  6.072                   9.07
```

The `f2d` and `negative_richardson` trend ratios (87.4 and −5.41) are artefacts
of dividing by a near-zero denominator — Prosser's own changes for those two are
+0.3 % (p = 1.0) and −7.5 % (p = 0.6), i.e. flat and insignificant. The honest
statement is "his is flat and not significant; ours is +26 % / +41 % and also
not significant", not "we are 87× his trend".

### 5.2 Which Phase-2 flags did the empirical comparison NOT flag?

These are formula or implementation findings that **cancel out of a
percentile-calibrated exceedance field**. Real, and correctly low priority.

| Phase-2 finding | Empirical status | Why it cancels |
|---|---|---|
| §3.1 meridional metric double-correction (≤ 0.7 %, latitude-structured) | **Not flagged.** `temperature_gradient` 0.913, `deformation` 0.999, `ngm1` 0.934 — the three most gradient-dependent diagnostics are all clean. | 0.7 % peak-to-peak is far below the ±0.02 read-off precision of Prosser's levels. **This is also a direct falsification of the cos φ hypothesis for the big disagreements** (see §5.4). |
| §4.5 Brown1 ≈ √0.3·f, a latitude field | **Not flagged** — `brown1` 1.087, `brown2` 0.913, both clean. | Both studies compute A13 identically, so the f-dominance is *shared*. It makes Brown1 a poor diagnostic (it carries little flow information) but it is not a replication defect. **§1 item 2 is hereby demoted**: it is a property of the published index, not a discrepancy. |
| §4.6 Brown2's missing length² factor | Not flagged (0.913). | An unknown *constant* L². Cancels exactly. |
| §4.3 CP's un-converted units and constant λ | Level *is* flagged (0.138) but **not by this** — a constant λ² and a constant kt² factor both cancel exactly. The CP disagreement must come from `Sv² − 2N²` itself. This narrows CP rather than explaining it. | constants |
| §4.9 `compute_rojak_diagnostics` returns DEF² | Not flagged — `deformation` 0.999, the single best-replicating diagnostic of the 21. | The production path applies the square root; and squaring is monotone on a non-negative field, so even the uncorrected path would give an identical exceedance. Confirms the §4.9 reading. |
| §3.2 polar amplification of ∂/∂x | Not visible per-diagnostic. | Affects the global threshold, not the box/globe contrast, and cos φ weighting suppresses it. Unfalsified either way by this table. |

### 5.3 Which diagnostics did the empirical comparison flag that Phase 2 found clean?

Seven diagnostics have a level ratio below 0.62. Phase 2 returned MATCHES for
five of them (`colson_panofsky`, `negative_richardson`, `f2d` — algebra —,
`endlich`, `ngm2`) and DIFFERS for one (`magnitude_pv`) and MATCHES-with-two-
ambiguities for one (`ubf`). And one diagnostic, `ncsu1`, is 6× **high** with a
MATCHES verdict.

**The disagreements are physical, not bugs — and the physics is sharper than
"the stencil".** `STATUS.md` §17.3 pre-registered a stencil test and got
0.93 (single-level, 11) against 0.73 (stencil, 10). Splitting the stencil group
by **what is being differentiated vertically** separates it much further:

| group | n | level ratios | median level | median trend |
|---|---|---|---|---|
| **A — no vertical derivative** | 11 | 0.193 … 1.121 | **0.934** | 0.99 |
| **B — vertical derivative of MOMENTUM only** (Sv = \|∂**v**/∂z\|): `vertical_wind_shear`, `ti1`, `ti2`, `brown2` | 4 | 0.850, 0.868, 0.878, 0.913 | **0.873** | 1.03 |
| **C — vertical derivative of a THERMODYNAMIC or DIRECTIONAL field** (N², ∂T/∂z, ∂ψ/∂z, ∂/∂θ): `negative_richardson`, `colson_panofsky`, `ncsu1`, `f2d`, `endlich`, `ngm2` | 6 | 0.138, 0.425, 0.494, 0.511, 0.615, **6.072** | **0.502** (0.494 excluding `ncsu1`) | 1.37 |

**Group B is statistically indistinguishable from group A.** The vertical
stencil, in itself, costs almost nothing. Every diagnostic that differentiates
only `u` and `v` vertically replicates as well as the ones that differentiate
nothing vertically. It is group C — and only group C — that fails.

**Why: products survive a stencil change, ratios and differences do not.**

*An earlier draft of this section attributed the split to the curvature of the
vertical profile — that `T(z)` and `θ(z)` have a tropopause kink at 200 hPa
while `u(z)` is smooth. That mechanism was tested on 2026-09-08 and is NOT
supported; it is withdrawn. What follows replaces it, and is supported by both a
measurement on this project's own data and an independent radiosonde study.*

Measured, by recomputing each vertical derivative on the 175/200/225 levels with
a 50 hPa centred stencil and again with a 25 hPa one-sided stencil (the finest
those three levels allow), and comparing the p3 / median / p97 of each:

```
quantity        stencil        p3      median       p97
Sv^2  (momentum) ratio 25/50   1.612    1.451      1.832
N^2   (thermo)   ratio 25/50   1.210    1.060      1.257
dT/dz (thermo)   ratio 25/50   0.878    2.831      1.915
Ri = N^2/Sv^2    ratio 25/50   0.620    0.854      0.708
```

Halving the stencil moves `Sv²` by roughly a *uniform* factor of 1.45–1.83 and
`N²` by a roughly uniform 1.06–1.26. **Each individual derivative is rescaled;
neither is reshaped much.** And a rescaling — any positive constant — cancels
EXACTLY out of a percentile-calibrated exceedance field (§5.2). That is why
group B survives: every member is a pure PRODUCT of derivatives, `Sv`,
`Sv·DEF`, `Sv(DEF−δ)`, `Φ·Sv²`, so a factor on `Sv` is a factor on the
diagnostic and cancels.

Group C members are not products. They combine two derivatives whose stencil
sensitivities are *different* — 1.45 against 1.06, a factor of seven apart — in
a way that does not cancel:

| # | diagnostic | how the derivatives combine | why the factor survives |
|---|---|---|---|
| 1 | `−Ri` | **ratio** `N²/Sv²` | the two factors divide, leaving 1.06/1.45 ≈ 0.73 |
| 3 | `colson_panofsky` | **difference** `Sv² − 2N²` | two terms rescaled differently, near-total cancellation (2N² ≫ Sv²), so the residual is dominated by the mismatch |
| 4 | `f2d` | **ratio** — `∂u/∂θ = (∂u/∂p)/(∂θ/∂p)` | momentum over thermodynamic, again 1.45 against 1.06 |
| 14 | `endlich` | derivative of an **angle** | `ψ = atan2(−v,−u)` does not rescale at all when `u, v` do |
| 16 | `ngm2` | `∂T/∂z` **crosses zero** near the tropopause | measured median ratio **2.83**, the largest of anything — a sign change is not a rescaling |
| 19 | `ncsu1` | **ratio** `1/Ri` | as `−Ri`, inverted (see §5.6) |

**Published corroboration, from Williams himself.** Williams & Storer (2022)
ran precisely this experiment and published it: they computed seven of these
diagnostics from ERA-Interim at its native resolution *and* from the same
ERA-Interim linearly interpolated onto HadGEM2-ES's much coarser horizontal grid
**and vertical levels** (p. 1427), then measured the spread. Their conclusion,
p. 1431, verbatim:

> "The best agreement between the climate model and reanalysis data evidently
> occurs for the purely **dynamical** diagnostics (i.e., Ellrod's turbulence
> index, wind speed × directional shear, and nonlinear balance equation
> residual), whereas the worst agreement occurs for the diagnostics **explicitly
> involving a thermodynamic component** (i.e., Richardson number,
> Colson–Panofsky index, frontogenesis function, and North Carolina State
> University index)."

Their Table 1 inter-dataset standard deviations (DJF, light CAT, 200 hPa):
Colson–Panofsky **1.00 %**, frontogenesis **0.75 %**, −Ri **0.68 %**,
NCSU1 **0.53 %** — against Ellrod TI1 **0.41 %**, wind speed × directional shear
**0.38 %**, UBF **0.22 %**. **This is the group B / group C split of the table
above, found independently, on different data, by the author of the method being
replicated.** All four of their thermodynamic diagnostics are among this
project's worst-disagreeing: −Ri 0.49, CP 0.14, f2d 0.43, ncsu1 6.07.

**Two of their results cut against §5.3's grouping, usefully.** W&S place
**UBF** among the *most* robust of the seven (0.22 %, best of all) and **wind
speed × directional shear** second-best (0.38 %). This project has those two at
level ratios of **0.19** and **0.51** — the worst and the fourth-worst. A
diagnostic that a published resolution experiment finds robust, and that this
pipeline finds badly off, is **not** suffering from resolution. That is strong
independent support for the §5.5 diagnosis of `ubf` (mixed operator provenance,
fix F2), and it moves `endlich` out of the stencil explanation entirely into a
new open item (§4.14, §6 F10).

**Independent corroboration on the mechanism.** Ko et al. (2023, *Atmos. Chem.
Phys.* **23**, 12589) compare ERA5 against global high-resolution radiosondes and find exactly
this asymmetry: ERA5 **underestimates wind shear by about 50 %** at comparable
vertical resolution, while its **N² is reliable to about 11 %**. They also find
that a sparser vertical grid systematically lowers the occurrence frequency of
subcritical Ri. So the two ingredients of Ri are damaged by very different
amounts — which is harmless for anything that multiplies them and decisive for
anything that divides them. Note this is the OPPOSITE of the withdrawn
mechanism, which blamed N².

**What is NOT established.** The probe above conflates three things — stencil
width, a shift of effective level (187.5 hPa versus 200), and first- versus
second-order accuracy — so it cannot separate group B from group C directly: the
p97 exceedance set flips 81 % for `−Ri` and 80 % for `vertical_wind_shear`
between the two estimators, i.e. *both* families are internally sensitive. The
product/ratio argument explains why that internal sensitivity nevertheless
cancels for one family and not the other, but it has been tested once, on one
small file. **F7 is what settles it, and §6 states the prediction.**

`colson_panofsky` at 0.138 is the extreme case and it should be: with a constant
λ and `N²/0.5 ≫ Sv²`, CP is a **pure monotone function of N²** and nothing else.
It is the cleanest single measurement of the stencil's damage to static
stability that this project has.

### 5.4 Could a latitude-structured error produce these? Testing the named suspects

The prompt's directive: a disagreement in exceedance *level* is a statement
about region-versus-globe contrast, i.e. about spatial structure, and
latitude-dependent terms are the prime suspects. Taking them one at a time:

| suspect | verdict | how it is falsified or supported |
|---|---|---|
| **cos φ metric terms (§3.1)** | **Falsified.** | The defect is real and measured (+0.42 % → −0.27 % over 30–75 °N) but it reaches `temperature_gradient`, `deformation` and `ngm1` most directly, and those are at 0.913, 0.999 and 0.934 — three of the five best. A ≤ 0.7 % error cannot produce 0.14. |
| **the Coriolis parameter f inside ζₐ** | **Falsified.** | `brown1` is ≈ √0.3·f by construction (§4.5) — the most f-saturated diagnostic of the 21 — and it is at 1.087. `nva`, which carries f through ∇(ζ+f), is at 0.973. Both studies compute f identically, so it cancels. |
| **tan φ terms in the spherical Jacobian (§4.17)** | **Live, for UBF only.** | `ubf` is the only one of the 21 containing a Jacobian, it is a residual of near-cancelling terms, and at 60 °N the distributed curvature term `u tan φ/a ≈ 8 × 10⁻⁶ s⁻¹` is comparable to `u_x` itself. Cannot be falsified from this table. Ranked below the mixed-provenance hypothesis below only because that one has a sharper test. |
| **the tropopause as a latitude-structured surface (§5.3)** | **Supported, and it is the main story.** | It predicts the group A / B / C split, the *direction* of the effect (< 1), and its concentration in the five thermodynamic/directional diagnostics. |

### 5.5 The two group-A exceptions, and why Phase 2 flagged both

`ubf` (0.193) and `magnitude_pv` (0.613) are the only two diagnostics with no
vertical derivative that are badly off. Phase 2 independently flagged exactly
these two, and for the same reason: **§3.4, archived versus computed ERA5
fields.**

**The fact that makes this decisive.** Prosser et al. (2023), p. 2, states the
complete input list: "zonal and meridional wind speed, dry bulb temperature, and
geopotential height were extracted … The 21 turbulence diagnostics were then
calculated from the extracted reanalysis fields." **Vorticity, divergence and
potential vorticity are not on that list.** Prosser therefore computed ζ, δ and
PV himself by finite differences from `u, v, T, z`. This project takes all three
from ERA5's archive. Ten of the 21 diagnostics are affected: `brown1`, `brown2`,
`ti2`, `magnitude_pv`, `vorticity_squared`, `ubf`, `horizontal_divergence`,
`ncsu1`, `nva`, `rva_magnitude`.

The archived-versus-computed axis does **not** explain the ensemble — its median
level ratio is 0.907 against 0.850 for the rest, i.e. no separation at all. It
is not a general defect. But it is decisive for the two diagnostics where the
*mechanism* is specific:

* **`ubf` (0.193).** UBF is a residual of an *equation*. Koch & Caracena (2002)
  §2 define it as "a pronounced residual in the **computed sum of the terms**" —
  the terms must be computed with one consistent operator or the residual
  measures the inconsistency between two operators rather than the atmospheric
  imbalance. Here `fζ` uses ERA5's spectrally-derived ζ (which carries
  small-scale power a 0.25° centred difference cannot represent) while `2J` and
  `∇²Φ` are 0.25° finite differences. Prosser's three terms all come from one
  operator. **This is the single most specific, most testable finding in the
  audit**, and Phase 2 reached it from the literature alone (§3.4, §4.17), before
  the level ratio was seen.
* **`magnitude_pv` (0.613).** Phase 2's verdict was **DIFFERS (defensibly)**:
  ERA5's archived `pv` is the full Ertel PV, whereas Sharman A18 (p. 283) —
  which is what Prosser must have implemented, having only `u, v, T, z` —
  defines `PV = −g ζₐ ∂θ/∂p`. These are two different fields, and they differ
  most where the horizontal vorticity and horizontal θ gradient are largest, i.e.
  in the jet, i.e. in the diagnostic's own upper tail. This is the one diagnostic
  where Phase 2 and the empirical comparison agree on **both** the existence and
  the identity of the difference: **case closed, fix specified** (§6, F3).

`horizontal_divergence` at 0.874 is worth one sentence of correction to the
project record: `STATUS.md` §17.4 clears it on the reasoning that "it is a
property of ERA5's archived divergence that both studies inherit". Prosser's own
Methods paragraph says he did not use ERA5's archived divergence. **The
conclusion (cleared) stands on the level ratio and the matching insignificance;
the stated reason does not.**

### 5.6 `ncsu1`, the one that goes the other way — not settled

`ncsu1` is 6.07× **high** in level and 9× high in trend, the only diagnostic
whose error has the opposite sign to everything else. `STATUS.md` §17.5 ranks it
first for exactly that reason and says the tropopause/N² story does not fit it.
**That is correct, and this audit does not settle it either.** Three candidates,
with the test that separates them:

1. **The Ri floor's geography (§4.19).** `1/MAX(Ri, 10⁻⁵)` is 10⁵ wherever
   `N² ≤ 10⁻⁵·Sv² ≈ 0`, against ~10⁻² for a typical `Ri ≈ 100`. Every floored
   cell is amplified by seven orders of magnitude and lands in the extreme tail,
   so NCSU1's upper tail may be *nothing but* the map of `N² ≤ 0`. That map is
   completely different on a 50 hPa stencil than on an 18 hPa one — a shape
   change, which a percentile calibration cannot absorb. Note this is the
   opposite-signed consequence of the same stencil mechanism: where §5.3's
   diagnostics lose their tail because low-N² cells are smoothed away, NCSU1's
   tail is a *threshold* effect on the sign of N², which coarsening can move
   either way. **Test: instrument `ncsu1()` to report the fraction of cells where
   the floor binds, globally and inside the box, and their geography.** If the
   floor binds on ≳ 3 % of global cells the entire global light threshold is set
   by it, which would be conclusive.
2. **`|∇ζ|` from archived vorticity.** NCSU1 is the only one of the 21 that
   contains `|∇ζ|` multiplied by two other flow factors — effectively a second
   spatial derivative of the wind. ERA5's archived ζ carries small-scale power
   that Prosser's `∂v/∂x − ∂u/∂y` at 0.25° cannot, and that power is concentrated
   in the midlatitude storm track (inside the box) rather than in the tropics
   (which dominate the global threshold). That predicts the right *sign*.
   Corroborating but weak: the two highest level ratios in group A,
   `vorticity_squared` 1.121 and `brown1` 1.087, are the two that take archived ζ
   most directly. Against it: `nva` (0.973) and `rva_magnitude` (0.901) also use
   ∇ζ and are clean, so this cannot be the whole story.
3. **The `MAX(u u_x + v v_y, 0)` clip.** It discards roughly half the domain, and
   the surviving half's regional-versus-global distribution is not the full
   field's. Phase 2 confirms the clip is verbatim in A36, p. 284, so this is not
   a bug — but it *is* a source of shape sensitivity, and it interacts with
   candidate 2.

Candidates 1 and 2 predict opposite responses to a stencil-width experiment
(candidate 1 moves, candidate 2 does not), which is what makes F5 in §6 worth
running.

### 5.7 Where both agree — case closed

| diagnostic | Phase-2 verdict | empirical | conclusion |
|---|---|---|---|
| `magnitude_pv` | DIFFERS — full Ertel PV, not A18 | 0.613, and the one that costs the 17th significant trend | **Closed. Fix specified (F3).** |
| `ubf` | MATCHES with two recorded ambiguities, one of which is mixed-provenance ζ in a residual | 0.193, worst in group A | **Closed on mechanism, open on which of the two ambiguities. Fix specified (F2), with a discriminating outcome.** |
| `f2d` | AMBIGUOUS-IN-SOURCE on the absolute value | 0.425, and his trend is exactly flat where ours is +26 % | **Closed as an open scientific decision, not a bug. Sensitivity run specified (F6).** |
| `deformation`, `brown1`, `brown2`, `ti1`, `wind_speed`, `ngm1`, `nva`, `rva_magnitude`, `vorticity_squared`, `temperature_gradient`, `horizontal_divergence` | MATCHES / UNVERIFIABLE-but-consistent | 0.87 – 1.12, all clean | **Closed. Nothing to fix.** |

---

## 6. Phase 4 — propagation, priority, and the evidence that would confirm each fix

Ranked by (diagnostics reached) × (confidence in the finding). Nothing here has
been applied.

### F1 — `spatial_gradient`'s meridional double-correction

* **Reaches:** every horizontal derivative — **13 of 21**, confirmed by running
  both versions rather than by counting call sites: `f2d`,
  `temperature_gradient`, `ubf`, `ncsu1`, `nva`, `rva_magnitude` directly
  through B1, and `brown1`, `brown2`, `ti1`, `ti2`, `deformation`, `ngm1`,
  `ngm2` through B2/B8. The other 8 come back bit-identical.
* **Sign and size:** ∂/∂y too large by `a/M(φ)` — +0.42 % at 30 °N, −0.27 % at
  75 °N. Monotone in latitude. Diagnostics whose gradient is meridionally
  dominated (`temperature_gradient` in the jet, `ngm2`) take the full factor;
  ones that combine both components take roughly half.
* **Confidence:** **highest of anything in this report.** Measured against two
  analytic fields with exactly known gradients, not argued.
* **Type:** **inherited** — one upstream block, thirteen symptoms. Fixing it in
  thirteen places would be the failure mode §4 of the brief warns about.
* **Proposed change (upstream, `rojak/core/derivatives.py:150–170`):** in
  `nominal_grid_spacing`, compute `dy` as the *map* distance `a Δφ` (equatorial
  radius × angular separation), matching the way `dx` is already taken at the
  equator, so that multiplying by `meridional_scale = a/M(φ)` converts it
  correctly. Equivalently, and locally, leave `dy` as the geodesic and do not
  apply `meridional_scale`. Authority: the derivation in §3.1 plus the two
  analytic tests. Local mitigation without touching the pin: a project-owned
  `spatial_gradient` wrapper that divides `dfdy` by `meridional_scale`.
* **Confirming evidence:** re-run the §3.1 analytic test — `∂f/∂y` must return
  `1.000000` at every latitude, as `∂f/∂x` already does. Then re-run the
  per-diagnostic comparison: **no level ratio should move by more than ~1 %**.
  Specifically `temperature_gradient` must stay near 0.913 and `deformation`
  near 0.999. If anything moves by more than a few per cent, the change did
  something other than what was intended.

### F2 — compute UBF's vorticity from the same operator as its other terms

* **Reaches:** 1 (`ubf`), which is the second-worst level ratio in the table.
* **Sign and size:** unknown in sign a priori; the level ratio should move from
  0.193 toward group A's median 0.93 if the mixed-provenance hypothesis is right.
* **Confidence:** high on the *principle* (Koch & Caracena 2002 §2 defines UBF as
  the residual of a *computed sum*, and a residual of near-cancelling terms is
  the worst possible place for two different operators), moderate on the *size*.
* **Type:** local.
* **Proposed change (`2_diagnostics.py::ubf`):** replace
  `zeta = _sel_level(ds, "vorticity", target_level)` with
  `zeta = vd[DV_DX] - vd[DU_DY]` reusing the `vector_derivatives` call already
  made three lines below, so that `∇²Φ`, `2J` and `fζ` share one geometry.
  Keep the archived field available behind a flag for comparison.
  Authority: Sharman A30, p. 284 and Koch & Caracena (2002) §2.
* **Confirming evidence:** `ubf`'s 1979 level ratio should rise from **0.19 to
  ~0.7–0.9** and its trend from +20 % toward Prosser's +31.6 %. **If it does not
  move**, the cause is the other recorded ambiguity — the metpy distribution of
  the spherical curvature term inside `J(u,v)` (§4.17) — and that becomes the
  next test, by recomputing `J` with the curvature term placed symmetrically.

### F3 — compute |PV| from Sharman A18 rather than using ERA5's archived Ertel PV

* **Reaches:** 1 (`magnitude_pv`), and it is the diagnostic that costs the
  replication its 17th significant trend (t = 1.97 against a critical 2.021).
* **Sign and size:** level ratio should rise from 0.613 toward ~0.9; trend from
  +26 % toward Prosser's +50.6 %; and the t-statistic should cross 2.021.
* **Confidence:** high that the two quantities differ (§4.10; ERA5 archives the
  full Ertel PV, A18 is the vertical term only); high that Prosser used A18
  (his Methods list, p. 2, contains no PV); moderate that this is the whole of
  the 0.613.
* **Type:** local.
* **Proposed change:** replace the `MAGNITUDE_PV` rojak diagnostic with
  `abs(rojak.turbulence.calculations.potential_vorticity(vorticity, theta))`,
  which already implements A18 exactly (`calculations.py:281–296`) and is
  currently unused by any diagnostic. Authority: Sharman A18 and the PV
  definition on p. 283; corroborated by Lee et al. (2023) Eq. 6, p. 3.
* **Confirming evidence, including the way it could fail:** this fix **moves
  `magnitude_pv` out of group A and into group C** — A18 contains `∂θ/∂p`, so it
  acquires the 50 hPa stencil. If the tropopause mechanism of §5.3 dominates, the
  level ratio could get *worse* rather than better. That is a genuine two-sided
  prediction and is the reason to run it as an experiment rather than a fix:
  **compute both, compare both against Prosser's panel (j), and keep the one that
  matches his definition** — which is A18, regardless of which scores better.

### F4 — move the deformation square root into `compute_rojak_diagnostics`

* **Reaches:** 0 in production; every non-production caller
  (`tests/test_analytic.py`, `tests/report_errors.py`, `4_verify.py`).
* **Confidence:** certain — read directly from `diagnostic.py:1378` and
  `2_diagnostics.py:999`.
* **Type:** local footgun, not an output defect.
* **Proposed change:** apply `sqrt` inside `compute_rojak_diagnostics` and delete
  the block in `compute_all_21`; or, if the passthrough must stay faithful to
  rojak, rename the key returned by `compute_rojak_diagnostics` to
  `deformation_squared` so that no caller can mistake one for the other.
  Authority: Sharman A17, p. 283 and Williams (2017) Table 2's units
  (10⁻⁶ s⁻¹, not s⁻²).
* **Confirming evidence:** `deformation` in the production output is unchanged
  bit-for-bit; the analytic-suite value changes by exactly its own square root.

### F5 — instrument, do not change, the NCSU1 Richardson floor

* **Reaches:** 1 (`ncsu1`), the largest single disagreement in the table.
* **Confidence:** this is a **measurement**, not a fix. A36 is transcribed
  verbatim (§4.19) and must not be "corrected".
* **IMPLEMENTED 2026-09-08** as `ada/ncsu1_floor_probe.py`. It reports the
  cos φ-weighted fraction of cells where the floor binds, the fraction of those
  that are convectively unstable, a latitude-band breakdown, and — the decisive
  number — **the fraction of NCSU1's own p97 exceedance set that sits on the
  floor**, with a printed verdict against the three thresholds below. Must be
  run on a **global** file, because NCSU1 is calibrated on global percentiles.
  On `era5_validation_subset.nc` (mid-latitude, summer, 2 timesteps) the floor
  never binds at all — which is the expected null and confirms the instrument
  reads zero correctly, nothing more.
* **Confirming evidence and what each outcome means:** if the floor binds on
  ≳ 3 % of *global* cells, it alone sets the light threshold and the 6.07 is
  explained; if it binds on ≪ 1 %, candidate 1 of §5.6 is dead and the `|∇ζ|`
  provenance hypothesis (candidate 2) becomes the leading one, testable by
  recomputing `ncsu1` with `ζ = dv_dx − du_dy` from `vector_derivatives`
  and checking whether the level ratio falls from 6.07 toward 1.

### F6 — revert the frontogenesis default from variant C to variant A

* **Reaches:** 1 (`f2d`), level ratio 0.43 and a +26 % trend where Prosser's
  panel (d) is flat at +0.3 %, p = 1.0.
* **Status change:** this was AMBIGUOUS-IN-SOURCE until 2026-09-08. It is not
  any more. **Williams & Storer (2022) Eq. (3), p. 1427** states the diagnostic
  as `F_θ = D/Dt |∂**u**/∂θ|²` — signed, un-normalised, no absolute value, no
  clip — which is the project's **variant A** up to the constant ½. W&S is
  Williams' own lineage, published one year before Prosser, and Prosser lists
  Williams as a co-author. It outranks every inference the project has made from
  a figure's axis limits.
* **Confidence:** high on what the equation says; the counter-evidence is
  described below and is not dismissed.
* **Proposed change:** set `F2D_DEFAULT_VARIANT = "A"`. Keep C reachable behind
  the existing flag. Record the citation in the constant's comment, replacing
  the 2026-08-30 note.
* **The counter-evidence, which must stay in the record.** Three things pointed
  at a one-sided quantity and none of them is explained by W&S Eq. (3):
  (i) Williams (2017) Fig. 1 plots frontogenesis on **0…300 × 10⁻⁹**, anchored
  at zero, in a figure where −Ri runs −300…0 and CP runs −45…−25, so Williams
  does plot signed diagnostics on negative axes and this one is not among them;
  (ii) Williams (2017) Table 2's p97 over W&J Table 1's median is **13.6**,
  which a material derivative centred near zero in a statistically stationary
  atmosphere cannot produce; (iii) measured p97/median was 754 for A and 22.7
  for C. **This is now an inconsistency between two Williams papers, not an
  ambiguity in this project's reading of one.** A written equation from the
  replicated lineage outranks an inference from an axis, so A is the default —
  but the tension is real and belongs in any write-up.
* **Confirming evidence:** run A and C through `ada/ab_compare.py` on the same
  month. Then check the trend: if A reproduces Prosser's flat, insignificant
  panel (d) where C gives +26 %, the equation and the replication agree and the
  matter is closed. If A gives the "second largest positive trend of the 21"
  that the project previously recorded, the inconsistency above is real and the
  honest output is to publish both with the discrepancy stated.

### F7 — the vertical stencil: not fixable, but the experiment is now sharper

* **Reaches:** 10 of 21 (groups B and C), and it is the *dominant* cause of the
  level deficit in group C.
* **Not fixable:** `STATUS.md` §11.6 establishes that 188/197/206 hPa are ECMWF
  L137 **model levels 73/74/75** (hybrid sigma-pressure, terrain-following),
  available only from the MARS tape archive; 175/200/225 is the finest stencil
  available at 200 hPa in the CDS pressure-level product. Independently
  confirmed here from Prosser (2023) p. 2, which names those three levels.
* **Proposed change:** none. But §5.3 sharpens the pre-registered stencil-width
  experiment (150/175/200/225/250, computing at 200 hPa with both a 50 hPa and a
  100 hPa stencil) into a much stronger prediction than
  "the deficit concentrates in the stencil group":

  > **Pre-registered:** the stencil-width curve must be **steep for group C**
  > (`negative_richardson`, `colson_panofsky`, `f2d`, `endlich`, `ngm2`) and
  > **flat for group B** (`vertical_wind_shear`, `ti1`, `ti2`, `brown2`).
  > Specifically, going from a 50 hPa to a 100 hPa stencil, group C's level
  > ratios should fall further below 0.5 while group B's stay within a few per
  > cent of 0.87. If both groups move together, the curvature explanation of
  > §5.3 is wrong and the stencil is a scale effect after all — which would
  > contradict the fact that a scale effect cannot change an exceedance level.

  This converts "our levels differ from Prosser's" from an apology into a
  two-point extrapolation with a falsifiable shape.

### F8 — latent fragilities, no output change

* `coriolis_parameter`'s degrees-versus-radians heuristic (`latitude.max() > π`,
  `calculations.py:244`) is wrong for any domain inside ±3.14° of the equator.
  Add an explicit units argument or an assertion.
* `Endlich`'s angle wrap-around silently stops working on a non-uniform vertical
  coordinate (§4.14). Add an assertion that the level coordinate is uniformly
  spaced.
* `colson_panofsky`'s λ is constant *only because* it comes from the ICAO
  standard atmosphere. Add a comment, and an assertion if λ is ever sourced from
  geopotential thickness, because at that point §4.3's stencil mismatch stops
  being inert.
* The `REFERENCE_TABLE["num"]` ordering is Williams & Joshi (2013) Table 1
  (ranked by % change), not Williams (2017) Table 2 (the canonical order Prosser
  Figure 4 uses). Both orderings appear in this project. Recording which is which
  in the table itself would remove a standing source of confusion.

### F9 — WITHDRAWN 2026-09-08

Proposed as a test of Kaplan et al. (2004) Eq. 3's full-advection NCSU1 against
Sharman A36's half-contraction. **Williams & Storer (2022) Eq. (7), p. 1428
reproduces A36 exactly, both `max()` operators included.** Williams' lineage
uses A36, the code implements A36, and the variant is therefore not the cause of
`ncsu1`'s 6× disagreement. The candidate is closed. §5.6's remaining two
hypotheses — the Ri floor's geography (F5) and `|∇ζ|` from archived vorticity —
are unaffected and F5 stands.

### F10 — Endlich: find out why a diagnostic that W&S find robust is off by half

* **Reaches:** 1 (`endlich`, level ratio 0.51).
* **Why this is new:** §5.3 originally filed Endlich under the vertical stencil.
  W&S Table 1, p. 1431 rules that out — "wind speed × directional shear" is the
  second-most robust of their seven to a coarsening of both the horizontal grid
  and the vertical levels (0.30–0.40 %). Something else is wrong.
* **Leading suspect:** the wind-direction gradient. Every other diagnostic in
  this pipeline uses a plain centred difference; this one routes through rojak's
  `angles_gradient` and `_WrapAroundAngleArray.__sub__`, which returns
  `min(|Δψ|, 2π − |Δψ|)` — non-negative by construction, so the sign of ∂ψ/∂z is
  destroyed. That is harmless under A25's `|·|` (§4.14) but it is not what a
  plain "centred second-order finite difference" (W&S p. 1428) computes, and the
  two differ wherever the wrap is *not* the shorter arc.
* **MEASURED 2026-09-08, and the lever is large enough.** Both estimators are
  now implemented (`FixSet.endlich_component_shear`). The alternative avoids the
  angle entirely: for the vector angle `α = atan2(v,u)`,
  `∂α/∂z = (u ∂v/∂z − v ∂u/∂z)/(u²+v²)` exactly, and since ψ differs from α by a
  constant offset and a sign, the whole diagnostic collapses to
  `|u ∂v/∂z − v ∂u/∂z| / √(u²+v²)` — no arctangent, no branch cut, nothing to
  wrap. On `era5_validation_subset.nc` the two disagree by **ρ = 0.9636, a
  median relative difference of 14.6 %, and a 67 % flip of the p97 exceedance
  set.** That is the right order of magnitude to move a level ratio of 0.51.
  **Neither is declared wrong**: rojak takes the secant of the *angle* across
  the 50 hPa layer, the closed form takes the angle-rate implied by the secants
  of the *components*, and W&S p. 1428's "centred second-order finite
  differences" describes both.
* **Proposed change:** none to the default. Run both through
  `ada/ab_compare.py --fix endlich_component_shear` on a real month.
* **Confirming evidence:** if the level ratio moves from 0.51 toward the
  no-vertical-derivative median of 0.93, the wrap-around handling is the cause.
  If it does not move, Endlich rejoins the unexplained list alongside `ncsu1`.
  Note that the 67 % flip means this is *not* a tie-breaker that can be left
  undecided: whichever estimator is chosen, roughly two thirds of #14's
  turbulent cells change identity, so the choice has to be made explicitly and
  cited rather than inherited from whichever library was imported.

### Summary ranking

| rank | fix | diagnostics reached | confidence | applied? |
|---|---|---|---|---|
| 1 | **F1** meridional metric | 15 | measured | no |
| 2 | **F2** UBF vorticity provenance | 1 (worst in group A) | high on principle | no |
| 3 | **F3** \|PV\| from A18 | 1 (costs the 17th significant trend) | high | no — run as an experiment |
| 4 | **F7** stencil-width experiment | 10 | n/a — measurement | no |
| 5 | **F5** NCSU1 floor instrumentation | 1 (largest disagreement) | n/a — measurement | no |
| 3= | **F6** revert frontogenesis to variant A | 1 | **W&S Eq. 3 states it** | no |
| 5= | **F10** Endlich angle-gradient choice | 1 | W&S rules out the stencil; 67 % flip measured | no |
| ~~—~~ | ~~**F9** NCSU1 variant~~ | ~~1~~ | **withdrawn** — W&S Eq. 7 is A36 | n/a |
| 7 | **F4** deformation footgun | 0 in production | certain | no |
| 8 | **F8** latent fragilities | 0 | certain | no |

---

## 7. Unresolved

Answers to the five named questions, plus everything the literature does not
settle.

**Q1 — the coefficient in the Brown index.**
It is **0.3**, and it multiplies **ζₐ² = (ζ + f)²**, the square of the vertical
component of *absolute* vorticity, inside a square root:
`Φ = (0.3 ζₐ² + D_SH² + D_ST²)^{1/2}`. What currently justifies it is **Sharman
et al. (2006) A13, p. 283**, which prints it and attributes the index to
**Brown (1973), *Meteor. Mag.* **102**, 347–360**. Brown (1973) is **not on
disk** and I have not read it. So the chain of authority is: code → Sharman A13
(read, verified) → Brown 1973 (unread). The same applies to the **1/24** in A14.
Everything downstream of Brown1 and Brown2 rests on a transcription of a paper
nobody in this project has read. **Cannot be resolved from available sources.**
Getting Brown (1973) would also settle whether A14's `ε` carries the length²
factor that Williams' J kg⁻¹ s⁻¹ units imply and A14 as printed does not.

**Q2 — the frontogenesis leading sign, and which variant Prosser used.**
**ANSWERED 2026-09-08, and the answer changes the code.** Williams & Storer
(2022) Eq. (3), p. 1427 — the only place in the replicated lineage where this
diagnostic is written as an equation — defines it as
`F_θ = D/Dt |∂**u**/∂θ|²`, i.e. `D/Dt[Q]`: **signed, un-normalised, no leading
minus, no absolute value, no clip.** That is the project's variant **A** up to
the constant ½. So: the algebra was already right (established independently
from A10 and Table B1); the *sign treatment* was not, and the executing default
of variant C (`|½ D/Dt Q|`) does not match the published definition.

Sharman A9's printed self-contradiction — its two sides differ by a sign, and
its RHS carries a `|∂**v**/∂θ|⁻¹` that its own Table B1 units exclude — is now
explained rather than merely noted: W&S's form has neither the minus nor the
normalisation, so both are artefacts of A9 as typeset.

**What remains genuinely unresolved** is not the equation but a conflict inside
Williams' own publications: W&S (2022) Eq. 3 is signed, while Williams (2017)
Fig. 1 plots this diagnostic on a 0–300 axis and Williams (2017) / W&J (2013)
give a p97/median ratio of 13.6 that a zero-centred material derivative cannot
produce. This audit takes the written equation over the inferred axis, and
recommends reverting to A (§6 F6) — but the conflict is real, it is in the
source and not in this project, and it should be stated in any write-up rather
than quietly resolved.

**Q3 — is the one-sided clipping on the vorticity-advection and NCSU
diagnostics specified by the source?**
Three different answers for three diagnostics, which is why the question is
worth asking:
* **NVA (#20): specified.** Sharman A37, p. 284 prints `MAX{…, 0}` explicitly.
* **NCSU1 (#19): specified by Sharman, but NOT by the diagnostic's origin.**
  A36, p. 284 prints both `MAX(u u_x + v v_y, 0)` and `MAX(Ri, 10⁻⁵)`, so the
  implementation is a faithful transcription of A36. But Kaplan et al. (2004)
  Eq. 3 — the source A36 itself cites, retrieved after this section was first
  written — has **neither**, and uses the full advection `U·∇U` rather than the
  half-contraction. So the honest answer is: the clipping is an *implementation
  choice made by Sharman*, inherited here in good faith, not a property of the
  index as published by its authors. §6 F9.
* **RVA (#21): not a clip at all, and correctly not.** No primary equation
  exists; the operator is an **absolute value**, authorised by the diagnostic's
  published *name* in Williams (2017) Table 2 and W&J Table 1 ("Magnitude of
  relative vorticity advection"). The code uses `np.abs` here and `.clip(min=0)`
  for NVA — the distinction is correct and is not accidental.

**Q4 — Colson–Panofsky's scale factor and units.**
A4 is dimensionally **m² s⁻²** (length² × s⁻²). Sharman Table B1 tabulates it in
**kt²**, Williams in **10³ kt²**. The code returns native SI and tags it
`m2 s-2`; **no scale factor is applied**, and none is needed for the
percentile-calibrated exceedance field, which is invariant to any positive
constant. The conversion, if a magnitude comparison is ever wanted, is
`CP[10³ kt²] = CP[m² s⁻²] × 3.7785 × 10⁻³` (1 kt = 0.514444 m s⁻¹). The formula
itself **matches A4 exactly** (§4.3). The *λ* in A4 is genuinely ambiguous —
"the local value of vertical grid increment Δz" does not say whether it is the
half-layer or the full centred stencil, and the code's λ (816.7 m, one 25 hPa
gap) is inconsistent with its own Sv²/N² stencil (50 hPa). Since λ here is a
constant, this is inert. **It would stop being inert the moment λ is made
state-dependent** (e.g. from actual geopotential thickness rather than ICAO).

**Q5 — the turbulent tail, for all 21, derived from the source.**
**Upper, for all 21.** Derived in §4.0 from Williams (2017) Table 1 (percentile
ranges 97.0 → 99.9, i.e. successively deeper right-tail cuts) applied to
Table 2's five-column ladders, every one of which is monotonically increasing —
including the two negative-valued diagnostics, −Ri (−15.4 → −5.9) and CP
(−29.3 → −22.2). Independently corroborated for those two by Sharman Table B1,
p. 285 (−Ri: −20 → 0.5; CP: 0 → 30 000 kt²). The `REFERENCE_TABLE` `sign: "+"`
entries agree with this derivation for all 21.

### Other ambiguities and gaps, recorded

1. **Sources not on disk:** Brown (1973); Kaplan et al. (2004) NASA
   CR-2004-213025 (the actual source of A36); Colson & Panofsky (1965);
   Reap (1996) (A28, A29); Endlich (1964) (A25); Bluestein (1992, 1993).
   Everything sourced from these rests on Sharman's transcription alone.
2. **Kaplan et al. (2005), on disk, does not define NCSU1.** It is a synoptic
   case-study paper. The index it names (p. 134) is "Ellrod index ÷ ipv", which
   is not A36. Any project record citing "Kaplan (2005)" as the source of the
   NCSU1 *equation* is citing the wrong paper.
3. **A2 offers two forms of N²** (θ and θₑ). No source says which Williams,
   W&J or Prosser used. The dry form is implemented.
4. **The Jacobian convention on a sphere** (§4.17) is unspecified by A30 and by
   Koch & Caracena, and the metpy distribution used here has an O(1) effect on
   J at 60 °N.
5. **The material-derivative sampling interval** for #4 is unspecified by A9 and
   by every downstream paper; 3-hourly centred differencing is an implementation
   choice with a 1000 km advective displacement per half-step.
6. **Brown2's missing length² factor** (§4.6): A14 gives s⁻³, the published
   tables give m² s⁻³. Unresolvable without Brown (1973).
7. **Williams (2017) and W&J both label #14 "10⁻³ rad s⁻¹"**, which is
   dimensionally inconsistent with A25's `s·∂ψ/∂z`. A published-unit error, not
   a code error.
8. **ERA5's archived PV is not A18's PV** (§4.10). The substitution is
   undocumented.
9. **`compute_rojak_diagnostics()` returns DEF², not DEF** (§4.9). Only
   `compute_all_21` corrects it. Any other caller silently gets the square.
10. **`coriolis_parameter` decides degrees-vs-radians by `latitude.max() > π`**
    (`calculations.py:244`). Correct for every domain this project uses, wrong
    for any domain lying entirely within ±3.14° of the equator. Latent.
11. **`Endlich`'s angle wrap-around depends on NumPy reducing a uniform
    coordinate to a scalar** (§4.14). Correct at 175/200/225; silently
    inoperative on a non-uniform vertical grid.
