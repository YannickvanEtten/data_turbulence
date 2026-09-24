# BATCH RESULTS — job 1124312, 2026-09-21

`jobs/27_final_checks.sbatch`, git `c7fbb11`, node220, **32 min 34 s, exit 0**,
all four steps exit 0. All six file hashes matched the tested versions.

**This document carries the outcomes of the eight pre-registered predictions in
`FINAL_DATASET.md` §4.3 and replaces every cell marked "pending" in its §5.**
Read it alongside that document, not instead of it.

> **Corrected 2026-09-21, later the same day, after the edits were applied to
> `FINAL_DATASET.md`.** Three defects in the first version of this note:
> (i) every number in §2 after the CI table was **SOG only** and was not
> labelled, so "71 % of the level gap" and "+0.4 pp" read as general claims when
> they are the strongest and the weakest severity respectively — all three
> severities are now given; (ii) **P2's outcome was not reported**, only
> referred to the log — it is now stated and resolved; (iii) the ensemble
> robustness figure was quoted for five class-C members after §6 had moved
> `endlich` to make six. The corrected text is below; nothing in §4, §5 or §6
> changed.

---

## 1. THE DATASET IS CERTIFIED

| | result |
|---|---|
| `derived/north_atlantic_audit` | **504 of 504** on disk, `audit_fixes = meridional_metric+four_variable_provenance` on **all 504**, **0 problems** |
| `derived/global_audit` | 12 of 12, same attribute on all 12, **0 problems** |
| `thresholds_2026-09-11.json` | n = **3,039,966,720** for every diagnostic; all signs `+`; 21 ladders finite and strictly increasing; provenance period names `derived/global_audit` |
| logs | 1,084 store writes covering 1,032 stores, no failures |

**P7 holds. There is no repair list.** The dataset is final as it stands.

**P1 holds.** The independent code path reproduces the job-1108148 table to
printed precision at MOG and SOG (largest difference 0.00 × half-last-printed-
digit). Every audit per-diagnostic number is usable.

**P6 holds, and it settles correction C2.** The `f2d` finite fraction at the
first and last step of every month is **(1.0, 1.0)** — the edge steps are
**one-sided, not NaN**. STATUS §15.3 was right; the `tail_thresholds.py`
docstring was wrong. **L4 now reads:** ≤ 0.82 % of steps carry a one-sided
difference, which passes a diurnal signal at 0.974 rather than 0.900; the count
is identical every year, so there is no trend effect.

**P2 holds in substance — this was left open in the first version of this
note.** The five convention-invariant diagnostics are **identical to printed
precision at SOG**, and identical in 34 of 35 cells at MOG. The exception is
`endlich`'s CI upper bound at MOG: 0.522526 against 0.522522, a difference of
**4 × 10⁻⁶, i.e. 0.0004 percentage points**. That is float summation order
between the xarray reduction path and the numpy one, not a difference in the
data. The check reported DIFFERENT because its tolerance was half the last
printed digit, which is tighter than the arithmetic noise between two code
paths; widen it to one printed digit before reusing it. **The restored
`data_prosser/` table is the 09-07 one, so decision D1 stands** — which was the
point of the prediction.

---

## 2. THE RESULT AGAINST PROSSER

**At LOG, all 21 diagnostics contain Prosser's published value inside our 95 %
interval.** Across the three severities, **59 of 63**.

| | inside our 95 % CI | same significance call |
|---|---|---|
| LOG | **21 / 21** | 19 / 21 (differ: `negative_richardson`, `colson_panofsky`) |
| MOG | 20 / 21 (outside: `ncsu1`) | 19 / 21 (differ: `colson_panofsky`, `ncsu1`) |
| SOG | 18 / 21 (outside: `endlich`, `ubf`, `ncsu1`) | 18 / 21 (differ: `colson_panofsky`, `ubf`, `ncsu1`) |

**Only three diagnostics ever fall outside**, and all three are already class C
or move there below.

### The level-gap decomposition — this strengthens the stencil argument

**Every row is labelled with its severity.** The first version of this note
quoted the SOG row alone, which is the strongest of the three.

| severity | stencil group (11) median / contribution | single-level group (10) median / contribution | ensemble gap | stencil share |
|---|---|---|---|---|
| LOG | 0.93 / **−0.059** | 1.00 / −0.039 | −0.098 | **60 %** |
| MOG | 0.85 / **−0.083** | 0.96 / −0.039 | −0.122 | **68 %** |
| SOG | 0.83 / **−0.110** | 0.96 / −0.045 | −0.154 | **71 %** |

Three things follow, and only the first was in the original note:

1. **The stencil group carries 60–71 % of the ensemble level gap, rising with
   severity.** This *revises* `FINAL_DATASET.md` correction C4, which found a
   roughly even split (−0.073 / −0.088, a 45 % stencil share) at **MOG on the
   baseline series**. The reason is the one C4 itself predicted: A18 moves
   `magnitude_pv` out of the single-level group and into the stencil group,
   which is now 11 against 10. **C4 should be rewritten** to say that the split
   is both severity-dependent and convention-dependent, and to quote both.
2. **The single-level contribution is nearly constant** (−0.039, −0.039,
   −0.045) while the stencil contribution nearly doubles. The part of the level
   deficit that grows into the tail is the part the stencil can explain.
3. **A quarter to a third of the gap is in diagnostics the stencil cannot
   touch**, at every severity. The stencil is not the only cause of the level
   deficit, and whatever else is at work is not severity-specific. Say this in
   the paper; a referee who only sees the 71 % will ask about the other 29 %.

Largest contributors, by severity:

| | largest negative | largest positive |
|---|---|---|
| LOG | `ubf` −0.034, `magnitude_pv` −0.020, `endlich` −0.015 | `ncsu1` +0.006 |
| MOG | `ubf` −0.035, `endlich` −0.027, `vertical_wind_shear` −0.015 | `ncsu1` +0.022 |
| SOG | `ubf` −0.030, `endlich` −0.029, `colson_panofsky` −0.027 | `vorticity_squared` +0.005, `ncsu1` +0.013 |

`ubf` is the largest single contributor at all three severities.

### The ERA5 / ERA5.1 bound holds (L2)

| severity | ensemble block anomaly | ensemble change, without → with a 2000–06 dummy |
|---|---|---|
| LOG | −3.6 % (t = −1.23) | +16.5 % → +17.1 % |
| MOG | −4.5 % (t = −1.01) | +37.2 % → +37.9 % |
| SOG | −5.0 % (t = −0.87) | +57.3 % → +58.0 % |

**Under one percentage point at every severity, and no t anywhere near
significance.** Across the 21 diagnostics the individual block anomalies span
−14.3 pp to +10.4 pp at LOG with |t| < 1.8 throughout and 18 of 21 within
±6 pp. The only diagnostic that ever reaches |t| ≥ 2 is `wind_speed`
(−32.4 %, t = −2.12 at MOG; −46.4 %, t = −2.16 at SOG), which is the sparsest
diagnostic in the set (0.05 % of cells at SOG) and is the evidence that blocks
of this size arise from variability.

This is not proof of no ERA5.1 effect. It is a demonstration that the effect is
smaller than this design can see, which is the claim L2 should make.

**New, and a further point for the audit conventions:** `magnitude_pv`'s block
anomaly has **vanished** — +0.7 % (t = +0.09) at SOG under A18, against −17 %
/ −22 % for the archived-PV version. The diagnostic most exposed to the ERA5.1
stratospheric correction is no longer the one showing a 2000–06 anomaly.

### Ensemble robustness

| severity | mean of 21 | without the 6 class-C members | difference |
|---|---|---|---|
| LOG | +16.5 % | +16.3 % | **−0.3 pp** |
| MOG | +37.2 % | +37.3 % | **+0.1 pp** |
| SOG | +57.3 % | +57.7 % | **+0.4 pp** |

Leave-one-out over all 21 spans **15.8–17.6 %** (LOG), **36.3–38.5 %** (MOG)
and **55.1–59.1 %** (SOG). `temperature_gradient` sits at the high end at all
three severities; the low end is `ngm1` at LOG and `brown2` at MOG and SOG. **No single member drives the result**, flagged or not.

(The scorecard's own line prints the five-member class-C list it was written
with; the six-member figures above add `endlich`, which §6 moves.)

---

## 3. THE PREDICTIONS

| | prediction | outcome |
|---|---|---|
| **P1** | new MOG table = job 1108148 | **HOLDS** — identical to printed precision |
| **P2** | five convention-invariant rows identical | **HOLDS in substance** — identical at SOG, 34/35 at MOG, the exception 4 × 10⁻⁶ of float summation order (§1). D1 stands |
| **P3** | `magnitude_pv` MOG level ratio 0.61 → 0.9–1.1 | **PARTLY** — 0.61 → **0.87**. Direction and most of the magnitude, slightly under the band |
| **P4** | `magnitude_pv` SOG trend moves toward his +90 %, inside our CI | **HOLDS** — +70 % [+46, +95] contains +90 % |
| **P5** | `ncsu1` at LOG: ratio in 0.75–1.33 and his +20.7 % inside our CI | **HOLDS** — see §4 |
| **P6** | `f2d` edges one-sided or NaN, nothing else short | **HOLDS** — one-sided |
| **P7** | sweep 0 items | **HOLDS** |
| **P8** | `ubf`'s tail not NH-weighted | **HOLDS, and the real answer is sharper** — see §5 |

**Eight for eight on substance; one, P2, was wrong in its framing.** A
tolerance of half the last printed digit is tighter than the arithmetic noise
between two code paths, so the check reported a difference of 0.0004 pp as
DIFFERENT. The prediction was right; the test was too sharp.

---

## 4. `ncsu1` IS EXPLAINED — P5 confirmed

| `ncsu1` | level ratio | our change | his | inside CI? | same sig? |
|---|---|---|---|---|---|
| **LOG** | **1.23** | **+20 % [+10, +30]** | **+21 %**, p = 9e-04 | **yes** | **yes** |
| MOG | 6.38 | +39 % [+20, +57] | +4 %, p = 0.7 | no | no |
| SOG | 2.57 | +47 % [+21, +73] | +10 %, p = 0.5 | no | no |

**NA share of the global tail, LOG → MOG → SOG:**

| | LOG | MOG | SOG |
|---|---|---|---|
| Prosser | **0.92** | **0.18** | **0.36** |
| ours | 1.13 | 1.12 | 0.93 |

At LOG the two agree almost exactly — +20 % against +21 %, level ratio 1.23,
same significance call. Above LOG his NA share **collapses and partly
recovers** while ours stays flat. That non-monotone signature appears in no
other diagnostic, and it is exactly what `FINAL_DATASET.md` §6 Q3 predicted
before the run.

**`ncsu1` is not wrong.** A36 multiplies by `1/max(Ri, 1e-5)`. On Prosser's
18 hPa stencil, thin statically unstable layers are resolved and the floor
binds often enough that his tail above LOG is floor-dominated — a population of
rare instability events, with a flat trend like his `negative_richardson`. On
our 50 hPa stencil those layers are averaged away, the floor binds on 0.0356 %
of cells, and our tail is dynamics-dominated at every severity.

**Caveat to hand to the econometrician:** *agrees with Prosser at LOG (+20 %
vs +21 %); above LOG his diagnostic is dominated by the Richardson floor on a
finer vertical stencil and ours is not, so the two are measuring different
populations. Ensemble member only. For EVT, note that on a finer stencil the
floor would plant an artificial mass in the extreme tail.*

Class stays **C** — but for a stated, evidenced reason rather than as an
unexplained outlier. This closes STATUS §5.6, open since 2026-09-09 after three
failed candidates.

---

## 5. `ubf` IS THE REMAINING PROBLEM, AND IT IS STRUCTURAL — P8

`ubf` is now the worst-fitting diagnostic, ahead of `ncsu1`.

| | LOG | MOG | SOG |
|---|---|---|---|
| level ratio | 0.39 | 0.19 | **0.08** |
| NA share ours | 0.82 | 0.35 | **0.10** |
| NA share his | 2.10 | 1.80 | **1.35** |
| trend | +13 % vs +18 % (in) | +20 % vs +32 % (in) | **+18 % [−6, +42] vs +44 % (OUT)** |

**The zonal profile of the reference-year tail says why.**

| `ubf` density | 90S–60S | 60S–30S | 30S–0 | 0–30N | 30N–60N | 60N–90N | NH/SH |
|---|---|---|---|---|---|---|---|
| light | 1.49 | 0.93 | 0.47 | 0.56 | 2.01 | 1.57 | 1.59 |
| moderate | 2.93 | 0.81 | 0.34 | 0.39 | 1.75 | 2.27 | 1.33 |
| severe | **4.67** | 0.67 | 0.28 | 0.26 | **1.01** | **3.70** | 0.99 |

**`ubf`'s extreme tail migrates to the poles.** Both caps take over as severity
rises, the 30–60N band drains from 2.01 to 1.01, and NH/SH falls to 1.0.

Compare a well-behaved diagnostic. `vertical_wind_shear` does the **opposite**:
30N–60N density **rises** 1.60 → 2.12 → 2.51, 60S–30S rises 1.68 → 2.18, and
both poles stay low (0.35 / 0.41 at severe). That is the jet-stream signature —
extreme shear concentrates in the midlatitude jets. `endlich` behaves the same
way (30N–60N 1.11 → 1.28, poles low). `ncsu1` is near-flat.

**`ubf` is the only one of the 21 whose extremes are polar rather than
jet-aligned.** That is structural, not a property of our box.

**The mechanism, stated as a hypothesis.** UBF is a residual of near-cancelling
terms, `−∇²Φ + 2J(u,v) + fζ − βu`. Its error concentrates where those terms are
largest and most nearly equal. Approaching the poles the spherical metric terms
diverge as cos φ → 0, β = 2Ω cos φ / R → 0 while f is maximal, and the
cancellation becomes worst-conditioned. A polar-dominated tail is what an
ill-conditioned residual would produce, and nothing else in the ensemble shares
the pattern. Our geometry is verified against Stokes' and the divergence
theorem at mid-latitudes (STATUS §4g); it has never been checked near the
poles.

**One limit on this evidence.** All six densities come from the reference year
alone. Whether the polar tail is a property of the diagnostic or of the year
2000 is not tested, and testing it would need a second calibration year, which
constraint 2 rules out. It does not change the verdict — class C either way.

**Caveat:** *our UBF's extreme tail is dominated by high-latitude cells where
the nonlinear-balance residual is numerically ill-conditioned; Prosser's is
jet-aligned. Ensemble member only. Do not use as a stand-alone indicator and do
not use its magnitudes.*

**This also clears `endlich` of a spatial explanation.** Its distribution is
jet-aligned and normal, so its 0.51 level ratio is consistent with the vertical
stencil, exactly as desk check D2 concluded.

---

## 6. CLASS CHANGES

Applying `FINAL_DATASET.md` §4.3's rule — move only if a LOG or SOG row shows
**both** Prosser outside our 95 % CI **and** a level ratio outside 0.75–1.33:

| | from | to | why |
|---|---|---|---|
| **`endlich`** | B | **C** | SOG: level ratio **0.46** and our +58 % [+37, +78] excludes his +36 %. Both criteria met |

**Nothing else moves.** Checked and not moved:

- `ncsu1`, `ubf` — already C, both criteria met at SOG, confirmed.
- `colson_panofsky`, `negative_richardson`, `f2d` — already C; each has
  Prosser **inside** our CI at every severity, so criterion (i) fails. They are
  C on sparsity and power, not on disagreement. Note `f2d`'s LOG agreement is
  excellent: ratio 0.97, +14 % against his +14 %.
- `magnitude_pv` (LOG ratio 0.72), `ngm2` (0.83 / 0.62 / 0.57) — outside the
  band but inside the CI at every severity, so criterion (i) fails. Stay **B**.

**Final counts: 13 class A, 2 class B (`magnitude_pv`, `ngm2`), 6 class C
(`negative_richardson`, `colson_panofsky`, `f2d`, `ubf`, `ncsu1`, `endlich`).**

Dropping all six class-C members moves the ensemble change by **−0.3 pp (LOG),
+0.1 pp (MOG), +0.4 pp (SOG)** — see the table in §2.

---

## 7. EDITS THIS REQUIRES

In `FINAL_DATASET.md` — **all applied 2026-09-21**:

- §5: every "pending" cell replaced with measured LOG/MOG/SOG values; `endlich`
  moved to class C; the `ncsu1` and `ubf` caveats rewritten from §4 and §5 here.
- §2 C4: rewritten. The split is severity- **and** convention-dependent:
  60 / 68 / 71 % stencil share at LOG / MOG / SOG on the final dataset, against
  45 % at MOG on the baseline series.
- §3 L4: the `f2d` month edges are **one-sided** (factor 0.974), not NaN.
- §3 L7: `ncsu1` is explained — P5 confirmed.
- §3 L10: `ubf` reads **structural**, with the polar-tail evidence, the
  ill-conditioning mechanism and the one-year limit on the evidence.
- §3 L15: new — agreement with Prosser degrades monotonically into the tail,
  which is where the econometrics wants to work.
- §4.3: P1–P8 recorded as resolved, with P2's framing defect noted.
- §6 Q2–Q5 and §7 rewritten; §8 added (what to keep, what to delete).

Still to do, in the other documents:

- `DECISION_RECORD.md`: C1 (the `f2d` factor is 0.900 diurnal, not 0.637, and
  is identical in Prosser), C4 as above, C5, C9, C11 (the trend bound is a
  **combined** bound, not a stencil bound).
- `FORMULAS_AND_DECISIONS.md`: D3 (the damping factor), D10 (`endlich`
  decided, sensitivity not measured on the series), §7 (close the open items).
- `essentials/DECISIONS.md` §7, the same.
- Note: `ada/prosser_published_s4.py` gained the 42 `abs=` values after the
  batch was submitted, so its sha256 no longer matches the table in
  `FINAL_DATASET.md` §4.2. The version marker is unchanged, so the job guard
  still passes — but do not quote the old hash.

---

## 8. WHAT REMAINS

Nothing that affects the dataset.

1. **`ubf`'s polar tail** is a real open question about our implementation, not
   about the data. It would be settled by checking the spherical geometry near
   the poles the way STATUS §4g checked it at mid-latitudes — Stokes' theorem
   on a high-latitude circuit. Cheap, and worth doing before anyone asks. It
   does not change `ubf`'s class either way.
2. **Prosser's Figure 3a**, 42 blue crosses, digitised for a year-by-year
   correlation against our ensemble series. Desk work. The only remaining check
   that could *strengthen* the replication claim rather than qualify it.
3. **The NA peaks-over-threshold archive** (L13) — one pass, ~23 GB, once the
   unit of observation is fixed.
4. **The ITvO backup email** (L14) — before anything is published.

**The data-gathering phase is closed.**
