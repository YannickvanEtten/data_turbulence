# STATUS — CAT replication

**Single source of truth for this repo.** The code cites audit IDs (`Q-CALIB-2`,
`Q-GLOBAL-1`, `Q-AGG-3`, `STATUS_24`, …) from working sessions whose notes are
no longer findable. Those references are not retrievable; this file replaces
them. When a decision gets made, it gets written here, not into a docstring.

Last updated: 2026-08-29 (stage 2 validated against Prosser; the trend check
run at nine seasons; phase 4 unblocked). §13 is the run ledger and the measured
disk inventory — read that first if the question is "what exists and what did
it cost".

---

## 1. What this project is

Replication of **Prosser et al. (2023, GRL)** — clear-air turbulence diagnostics
computed on ERA5 over the North Atlantic — as a correctness check on the
pipeline. The econometric research design is deliberately deferred until the
replication holds up.

Likely econometric direction (not yet fixed): **tail-behaviour modelling** —
whether the shape of the CAT distribution's upper tail changes over 1979–2020,
not merely its mean. See §6 for what this implies for the production run.

---

## 2. Workflow

```
1_download.py       original single-day reference downloader (untouched)
1_download_hpc.py   ADA entrypoint: resumable, atomic, integrity-checked
download_plan.py    pure request planning (LOCKED science constants)

2_diagnostics.py    ALL diagnostic logic — 21 diagnostics + data prep
3_pipeline.py       thin orchestration: load → compute → stats → save
                    also hosts run_layers_2_to_6 (calibration → trend)
4_verify.py         correctness harness for the hand-written diagnostics
5_explore.ipynb     visual exploration

calib_weighted_percentile.py   Layer 2  cos(phi)-weighted percentile
calibration.py                 thresholds as a file: compute/save/load/compare
aggregate.py                   Layers 3-4  exceedance-first, then average
annual_aggregate.py            Layer 5  leap-aware annual normalisation
trend.py                       Layer 6  per-gridpoint OLS 1979–2020
chunk_stitch.py                month-boundary overlap buffer for d/dt

ada/                the ADA drivers — see §11.11
jobs/               the SLURM scripts — see §13.4

old_code/           DEAD. Superseded. Four competing versions of the
                    diagnostics live here. Nothing imports it. Do not
                    consult it to answer a question about current behaviour.
```

---

## 3. State of play

**Done.** All 21 diagnostics implemented: 14 taken from `rojak` after a source
audit, 7 rewritten by hand from Sharman (2006) Appendix A because rojak was
wrong. Root bug found: rojak computes Richardson as N²/Sv rather than N²/Sv²,
which silently broke Colson–Panofsky and NCSU1 as well. Layers 2–6 written and
wired into `3_pipeline.py`. HPC downloader written.

**The stack has now been exercised end to end on real ERA5 and validated
against a published result.** 39 months downloaded and put through all 21
diagnostics; global year-2000 thresholds computed and applied; the resulting
North Atlantic DJF frequencies reproduce Prosser (2023) Table 1 to within 21 %
(§11), and the fitted 1979→2020 trend contains Prosser's published change
inside its 95 % interval at all five severities (§12).

**Not done.** The *production* series does not exist: 39 of 504 months. Phase 4
is now a scheduling problem rather than a correctness one — see §12.7 for what
gates it and §13.2 for the measured cost.

The ADA environment is verified and documented — see §10. The cfgrib concern
its docstring raised is **resolved** (§10.8).

---

## 4. Verification status of the 21 — ALL COVERED as of 2026-08-26

The old three-tier split (5 strong / 14 read-only / 2 unverified) is closed.
Every diagnostic now has numerical evidence behind it, from one of two
independent methods.

| Method | n | What it establishes | Diagnostics |
|---|---|---|---|
| **Analytic (manufactured solution)** — `tests/` | 16 | Computed output matches a **closed-form** expected value on an atmosphere built with known derivatives. Compares against truth, not against another implementation. | `magnitude_pv`, `brown1`, `temperature_gradient`, `horizontal_divergence`, `vertical_wind_shear`, `endlich`, `deformation`, `wind_speed`, `ngm2`, `rva_magnitude`, `nva`, `vorticity_squared`, `ti1`, `ngm1`, `ti2`, `f2d` |
| **Cross-implementation** — `4_verify.py` | 5 | Hand-written version agrees with an in-memory corrected rojak. **All five now match at ρ = 1.0000** (2026-08-26). | `negative_richardson`, `colson_panofsky`, `ubf`, `ncsu1`, `brown2` |

Measured accuracy is in `verification/<date>/analytic_verification.csv`.
Fourteen of the sixteen analytic cases land under 0.2 % relative error; see
§4b for the two that do not and why.

**What this does and does not establish.** It establishes that each diagnostic
correctly computes **its stated formula** — every expected value was derived by
hand from the formula in the relevant docstring, each of which cites its source
paper. It does **not** establish that those formulas are the ones Prosser (2023)
intended. That is a literature question, needs the PDFs in `Articles/`, and is
tracked separately in §7. Conflating the two is how a project convinces itself
it has verified more than it has.

### 4b. The two that are not exact, and why

- **`ngm2` — 2.8 % (accepted).** `NGM2 = |dT/dZ| × DEF`. Temperature is not
  linear in pressure, so rojak's three-level `.differentiate` carries a real
  O(dp²) truncation error at the project's 175/200/225 hPa spacing.
  `test_ngm2_vertical_error_is_truncation_not_formula` proves this is
  discretisation rather than a formula error by halving the level spacing and
  showing the error fall by ~4× each time. Largely absorbed by the percentile
  calibration, but it should be quoted with the caveat.
- **`endlich` — 0.9 % (accepted).** Wind direction is an angle, so its vertical
  derivative inherits error from both wind components; 0.9 % is consistent with
  the metric floor compounding through the arctangent.

### 4c. Two things the analytic suite discovered

1. **Deformation and its whole family depend on curvature corrections.** DEF,
   TI1, TI2, NGM1 and NGM2 initially disagreed with the expectation by 20–40 %.
   The fault was in the expectation: differentiating a *vector component* on a
   curved earth adds metric terms that differentiating a scalar does not, and
   the first draft omitted them. rojak applies them correctly (metpy's
   map-factor correction, reducing to tan(φ)/M on a sphere). This is positive
   evidence about rojak — **and a standing warning: any future hand-written
   diagnostic that differentiates u or v must use `vector_derivatives()`, never
   the scalar `_grad()`.** `ncsu1` already does; `rva` correctly does not need
   to, since it differentiates the scalar vorticity.
2. **`ti2` uses ERA5's own divergence product, not the wind-derived
   divergence.** rojak's `TurbulenceIndex2` accepts `du_dx` and `dv_dy` and
   then never uses them in `_compute`. On real ERA5 the two divergences are
   nearly identical (both descend from the same spectral wind), so this is not
   a defect — but it is invisible in the source and is now pinned by a test.

### 4d. A data caveat the suite exposed for F2D

A centred difference of a signal at angular frequency W sampled at spacing h
returns the true derivative times `sin(Wh)/(Wh)`. For a **diurnal** signal on
**3-hourly** ERA5 that factor is **0.637** — the material derivative in F2D
understates it by 36 %. This is a property of the sampling, not of the code, and
no code change can repair it. F2D should be read as systematically damped for
anything sub-diurnal. Worth an explicit sentence in any write-up.

### 4e. The gradient foundation

`test_gradient_metric_matches_exact_ellipsoid` checks the operator every other
gradient-based diagnostic sits on. rojak reaches ∂/∂x via a nominal equatorial
grid spacing times PROJ's `parallel_scale`; compared against the exact WGS84
metric computed from the defining constants, it agrees to **6 significant
figures** zonally and **0.2 %** meridionally, with cross-terms at 1e-20. That
0.2 % is the noise floor for the whole suite.

### 4g. Two geometry bugs found in `ubf` (#13), fixed 2026-08-26

Found by asking *which derivative operator is physically correct*, and settling
it with theorems rather than with the `_grad()` docstring. Both bugs were
present in **rojak as well**, which is why the cross-check in `4_verify.py`
could never have caught them — the two implementations shared the mistake.
That is the general lesson: agreement between two implementations is only as
good as their common assumptions.

**The arbitration.** Two theorems, neither of which involves choosing a
derivative operator:

| Theorem | Physical value | Flat operator | Curvature-corrected |
|---|---|---|---|
| Stokes (circulation ÷ area = vorticity) | 5.292e-06 s⁻¹ | −2.07e-08 → **100 % wrong, wrong sign** | 5.310e-06 → 0.34 % |
| Divergence (boundary flux ÷ area = ∇²) | 4.4695e-11 s⁻² | 4.473e-13 → **99 % wrong** | 4.4844e-11 → 0.33 % |

Neither correction is a refinement. A zonal flow on a sphere carries curvature
vorticity `u·tan(φ)/M` even with **no shear at all**, and the scalar operator
cannot see it.

**Bug 1 — ∇²Φ was not a Laplacian.** `d2phi = _grad(dphi_dx)[0] + _grad(dphi_dy)[1]`
omits the spherical `−tan(φ)/M · ∂Φ/∂y` term. The Laplacian is the *divergence
of the gradient*: `grad(Φ)` is fine with the scalar operator because Φ is a
scalar, but its **divergence** is a vector divergence and needs the correction.

**Bug 2 — the Jacobian used scalar gradients of vector components.** `J(u,v)`
was built from `_grad(u)`/`_grad(v)`; rojak's own
`CATData.jacobian_horizontal_velocity()` uses `velocity_derivatives()`. The
Jacobian itself was 50 % wrong, though `2J` is only ~0.5 % of the residual.

**Why it matters more here than anywhere else.** UBF is a **residual** — a
near-cancellation of large terms. The `fζ` term uses ERA5's own vorticity,
which is the true spherical vorticity, so computing the other terms in a
flat-earth geometry meant the residual was partly measuring that inconsistency
rather than atmospheric imbalance. And the error scales with `tan(φ)`, so it is
**latitude-structured** — exactly the class of error that does *not* wash out of
a percentile calibration (§5.4) and *does* contaminate spatial patterns and
their trends. **See §12.6: UBF is now the one diagnostic behaving anomalously
on real data, and this is the first place to look.**

**Measured impact on real ERA5** (`era5_validation_subset.nc`):

| | before | after |
|---|---|---|
| median UBF | 3.635e-09 | 3.632e-09 |
| median abs change | — | 3.05e-10 s⁻² (~1.9× the W&J published median) |
| Spearman ρ (before vs after) | — | **0.99109** |
| cells changing > 10 % | — | **41.2 %** |

Rank order changes, and exceedance counting sees ranks. On the *synthetic*
field the change looks like 0.1 % with ρ = 1.00000 — an artifact of that field
not being near geostrophic balance, so UBF there is dominated by ∇²Φ instead of
being a small residual. Real data is the meaningful measurement.

**Confirmation.** `4_verify.py`'s cross-check for UBF went from
ρ = 0.9969 `CLOSE` to **ρ = 1.0000 `MATCH`** once the same correction was
applied to the monkey-patched oracle. All five hand-written diagnostics now
match exactly.

**Guarded by** `tests/test_geometry_theorems.py` — the two theorems plus a
regression test that fails if `ubf()` ever reproduces the old flat formula.

**Standing rule, now with evidence behind it:** any diagnostic taking a spatial
derivative *of* a wind component, or a divergence/Laplacian of any vector,
must use `vector_derivatives()`. The rule is about *what is being
differentiated*, not about which symbols appear: `rva` multiplies by u and v
but differentiates the scalar ζ, and `frontogenesis_isentropic` differentiates
the rotational invariant `(du/dθ)²+(dv/dθ)²` — both correct as written. The
deprecated `frontogenesis_2d` has the same violation as `ubf` did, but nothing
calls it.

### 4f. Running it

```bash
python -m pytest tests/ -v              # 56 tests, ~20 s
python tests/report_errors.py           # writes verification/<date>/analytic_verification.csv
python 4_verify.py                      # the cross-implementation check
```

### 4a. Measured results — run of 2026-08-26

Input `era5_validation_subset.nc`. Full output in `verification/2026-08-26/`.

| Diagnostic | median ratio | Spearman ρ | verdict |
|---|---|---|---|
| `negative_richardson` | 1.000 | **1.0000** | MATCH |
| `colson_panofsky` | 1.000 | **1.0000** | MATCH |
| `ncsu1` | 1.000 | **1.0000** | MATCH |
| `brown2` | 1.000 | **1.0000** | MATCH |
| `ubf` | 0.993 | **0.9969** | CLOSE |

**Two gaps this run exposed, both properties of the validation FILE, not the
code:** `f2d` was computed over an empty array (the subset has only two
timesteps and the material derivative drops the first and last), and the
negative-Ri branches never fired (Ri > 0 everywhere in that file).

**Both are now closed on real data** by the 39 months in §13 — see §10.5 of
`CALIBRATION_REFERENCE.md`. Replacing `era5_validation_subset.nc` with a
≥3-timestep winter subset would still be tidier for the standalone harness, but
it is no longer load-bearing.

---

## 5. What can actually invalidate the results

Ranked by impact on the trend estimates, which are the only numbers that matter.

1. **Rank-order errors** — a wrong term changes *which* cells are flagged, not
   just by how much. The F2D `dv_dy`→`dv_dx` bug was this kind.
2. **Sign / tail errors** — the UBF sign flip was this. Cheap to re-derive all
   21 `REFERENCE_TABLE` sign entries from the source papers; unrecoverable
   downstream if wrong.
3. **Time-varying bias** — anything whose error drifts across 1979–2020 goes
   straight into the slope. Includes ERA5's own observing-system changes as
   satellite instruments enter the assimilation stream. A standard referee
   objection; needs an explicit answer before publication.
4. **Interannual variability swamping the trend in any short sample.** Measured
   in §12.4: across nine DJF seasons the MOG frequency spans a factor of 1.65,
   while the fitted 41-year change is 43 %. This is why Prosser fits 42
   consecutive years and describes his Table 1 as a guide "in the absence of
   interannual variability", and why no subsample can substitute for the full
   series.
5. **Constant magnitude offsets — LOW PRIORITY.** Thresholds are percentiles of
   the data itself, so a uniform scaling `c·D` scales the threshold by the same
   `c` and leaves the exceedance field **identical**:

   ```
   D′ ≥ percentile_p(D′)  ⟺  c·D ≥ c·percentile_p(D)  ⟺  D ≥ percentile_p(D)
   ```

   This holds for any strictly increasing transformation. So the
   `1–2 ORDERS OFF` entries in `cat_outputs/comparison_table.csv` are a smell
   test against Williams & Joshi's published medians, **not** a correctness
   criterion for Prosser's method.

   **UPDATED 2026-08-29 — the proviso has been measured, and it does not
   hold.** §5.4 flagged that the level-substitution argument works only if the
   damping is roughly constant. It is not. The deficit against Prosser runs
   from 13 % at LOG to 31 % at SOG, and the trend excess mirrors it, 3 % to
   43 % (§12.5). The percentile argument holds for the **bulk** and leaks in
   the **tail** — which is exactly where the econometrics wants to live. This
   is now a documented property of the dataset, not an open defect, because
   §11.6 shows it cannot be closed on CDS.

---

## 6. Plan

### Phase 0 — legibility  ✅ COMPLETE (2026-08-26)
- [x] **Fixed a crash that stopped the whole pipeline running.** Both
      `3_pipeline.py` and `4_verify.py` load `2_diagnostics.py` by file path
      with `importlib`, but never registered it in `sys.modules`. Because
      `2_diagnostics.py` uses `from __future__ import annotations`, `@dataclass`
      resolves its field annotations through `sys.modules[cls.__module__]` —
      which was `None` — and the import died on Python 3.12 with a bare
      `'NoneType' object has no attribute '__dict__'`.
- [x] Fix the hardcoded `/home/claude/work/phaseB_out` path in `4_verify.py`.
- [x] Add `.gitignore` so `git add .` stops being dangerous.
- [x] This file.
- [x] Made `summarize()` in `4_verify.py` empty-safe.
- [x] Ran `4_verify.py`; evidence in `verification/2026-08-26/`.

### Phase 1 — close tiers B and C  ✅ MOSTLY COMPLETE (2026-08-26)
- [x] Manufactured-solution suite in `tests/` covering all 16 diagnostics that
      had no numerical evidence. See §4.
- [x] Convergence tests — the stronger evidence, since they show the error
      falls as the grid refines.
- [x] Verified the gradient metric against the exact WGS84 ellipsoid.
- [x] The negative-Ri branch coverage gap is closed on real data (§4a).
- [ ] Build a better validation subset for `4_verify.py` — now cosmetic.
- [ ] **Confirm the formulas against the papers** (§7). Partly done: units are
      closed, `ti2` is closed. The Sharman Appendix A derivations remain.

### Phase 2 — pilot: one month, end to end, on ADA  ✅ COMPLETE
- [x] cfgrib opens all 7 variables as one hypercube — §10.8.
- [x] **Compute-node egress — CONFIRMED (§10.7).**
- [x] `chunk_stitch.py` filename bug — FIXED.
- [x] Job scripts written — §13.4.
- [x] ERA5 pressure-levels Terms of Use accepted on CDS.
- [x] Pixi env built on the project share (§10.10).
- [x] Smoke test, first download, first diagnostics (§10.11–12).

### Phase 3 — decisions  ✅ RESOLVED 2026-08-29
- **Pressure levels — CONSTRAINT, NOT A CHOICE.** 188/197/206 hPa are *model*
  levels 73/74/75, MARS-only (§11.6). 175/200/225 is the finest stencil
  available at 200 hPa on CDS. Accept and document.
- **Calibration domain — GLOBAL.** Demonstrated feasible via day sub-sampling
  (§11.8), which is what made a global year fit in memory at all.
- **Box — NORTH ATLANTIC** for the production series; global only for the
  year-2000 calibration.
- **What the production run persists — ZARR ONLY, float32, 200 hPa.**
  `save_outputs()` writes NetCDF *and* Zarr, which doubles storage for nothing;
  `ada/diagnostics_global.py` writes zarr only and is what production uses.
  Measured cost in §13.2: **~470 GB, 19 % of the 2.5 TB share.** The earlier
  796 GB and 1.32 TB figures assumed dual-format or float64 and are superseded.
  Taking 19 % rather than >50 % also removes the conversation with the Koopman
  group that §10.6 flagged.

### Phase 4 — the 42-year run
CDS-queue bound (504 monthly requests), not compute bound. Measured sizing and
timing in §13.2 (~470 GB, ~36 h of download); gating items in §12.7 — there are
none on the correctness side.

**How to start it.** The array spec and its throttle are declared *inside* the
script, so submit it bare:

```bash
cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
mkdir -p logs
sbatch jobs/01_download.sbatch          # NOT --array=..., see §11.9
```

Passing `--array` on the command line silently discards the `%4` and breaks
roughly two thirds of the downloads. The 27 months already on disk detect their
own output and exit in seconds, so resubmitting after any interruption is free
and needs no list of which months failed.

Then, once the raw files are down, the derived side:

```bash
sbatch --array=0-503%8 jobs/07_diagnostics_na_trend.sbatch   # after editing MONTHS
```

— `07` currently hard-codes 27 DJF months. Production needs a driver over all
504; that is a small edit to the same pattern, not a new script.

### Phase 5 — econometrics

**Design constraint, decided in phase 3, not phase 5.** Layers 3–6 are the
*replication check*, not the research dataset. `aggregate.py` produces **binary**
0/1 exceedance fields by design — magnitude is discarded. That is correct for
Prosser and **fatal for extreme-value work**: a GPD cannot be fitted to an
indicator variable, and peaks-over-threshold needs the size of each excess.

| What is saved | Size | Supports |
|---|---|---|
| Everything raw, 3-hourly, float32 zarr | ~287 GB (measured, §13.2) | anything |
| Monthly per-gridpoint quantiles (p50/p90/p99/p99.9/max) | ~23 GB | trends, coarse tail work |
| Peaks-over-threshold archive (top 1%, magnitudes kept) | ~23 GB | EVT / GPD fitting |
| Annual exceedance grids only (current behaviour) | trivial | replication only |

Keeping the full derived field costs 287 GB and makes the middle two derivable
at any time. That is the plan.

---

## 7. Open: do the formulas match the papers?

The analytic suite verifies code against docstring. It cannot verify docstring
against literature.

**Closed since 2026-08-28:**

- **All 21 `REFERENCE_TABLE` units and `wj_median` values** were checked line by
  line against both papers; all 21 match. Colson–Panofsky's *scale factor* was
  wrong (§11.5) but its recorded unit was right.
- **`ti2` sign convention — CLOSED.** Ellrod & Knapp eq (9)/(10) define
  `CVG = −(∂u/∂x + ∂v/∂y)`, so `DEF + CVG = DEF − δ`, exactly as implemented.
  The DEF and VWS forms are corroborated by the same reading.
- **`f2d`'s units argument — CONFIRMED.** Sharman's own Table B1 units require
  the squared form, which contradicts the printed `^(1/2)` in A9. The table
  wins.

**Still open, in rough order of what a mistake would cost:**

- **`ubf`.** Promoted to the top of this list on 2026-08-29. Not a literature
  question any more but an empirical one: it is the only diagnostic whose trend
  breaks the pattern every other one follows (§12.6). Its geometry was rewritten
  on 26 August (§4g) and it is a residual of near-cancelling terms, so it is the
  most fragile of the 21 by construction.
- **`f2d`'s leading minus sign in A9** is not in the implementation and has not
  been checked. A sign question, i.e. §5's category 2.
- **`brown1` coefficient.** `sqrt(0.3 ζ_a² + D_sh² + D_st²)`. The 0.3 is
  Brown (1973) via rojak's docstring. **Brown (1973), *Meteorological
  Magazine* 102, 347–360 is not in `Articles/`** and must be obtained.
- **`nva` and `ncsu1` clipping.** Both apply `max(·, 0)`. Confirm Sharman (2006)
  A36 and the NVA definition really are one-sided.
- **`brown2`** sits **five orders** from its published value in a way that
  resolution does not explain, and it is absent from Sharman's Table B1 so the
  independent anchor in §11.12 does not cover it. *Partially reassured by
  §12.6*: its level and trend both sit mid-pack among the 21, so its **ranks**
  — which is all the exceedance counting uses — look healthy even though its
  magnitude does not.
- **Every `sign` entry in `REFERENCE_TABLE`.** Two were already found wrong
  (`colson_panofsky`, `negative_richardson` were tagged `"either"` and corrected
  to `"+"`). A cheap, high-value re-derivation.

## 8. Open questions

- Intended econometric unit of observation? (route × month, gridpoint × month, …)
- Supervisor expectations and deadline.
- **Whether `yen230` may use the `unlimited` QOS** (§11.9). Ask
  `itvo.ucit@vu.nl`. Highest-leverage question before phase 4.
- Snapshot retention / backup policy on the Koopman share.
- The `STATUS_24` / `Q-*` notes, if they ever turn up. Not a blocker.

---

## 9. Repo facts worth remembering

- `cat_outputs/` is 421 files / 114 MB. Now git-ignored except its small CSVs.
- `.git/config` has an `[lfs]` section but there is no `.gitattributes`, so LFS
  is half-configured and will not catch large files. Do not rely on it.
- Everything lives in OneDrive, so uncommitted work is backed up and versioned.
- `rojak` requires Python ≥ 3.12.
- **`logs/` is not yet git-ignored.** Should be.
- **`pixi.lock` is deliberately not committed.** The pin that matters is the
  rojak SHA in `pixi.toml` (§10.4) plus the numpy/numba/sparse reference set;
  `pixi install` from `pixi.toml` reproduced the environment on ADA from
  scratch. A lockfile would be stricter but is not what makes this
  reproducible.

---

## 9a. Repository and data layout

### The three-stage science plan this serves
1. **Pilot** — one North Atlantic month, end to end. ✅ §10.11–12
2. **Calibrate** — global year 2000 → severity thresholds → **compare against
   the published tables**. ✅ §11
3. **Analyse** — North Atlantic 1979–2020, thresholds applied, trends computed.
   Partially: 9 DJF seasons as a trend check (§12); the 504-month series is
   phase 4.

Stages 2 and 3 use *different datasets*. `3_pipeline.run_layers_2_to_6` already
took `calibration_fields` separately from `diagnostic_fields`; what was missing
was the plumbing to produce both, and a way to carry thresholds between them —
now `calibration/thresholds_YYYY-MM-DD.json`.

### On-disk layout
```
/scistor/SBE-EDS-ClimateKoopman/yen230/
├── data_turbulence/     the git repo — CODE ONLY, nothing large
├── env/                 the pixi environment (pixi.toml)
├── raw/
│   ├── north_atlantic/  era5_na_YYYY-MM.grib
│   └── global/          era5_glob_2000-MM_d01-09-17-25.grib
├── derived/
│   ├── north_atlantic/  diagnostics_na_YYYY-MM.zarr
│   └── global/          diagnostics_glob_2000-MM.zarr
├── calibration/         thresholds_YYYY-MM-DD.json   <- stage 2 -> 3 handoff
├── results/             annual probabilities, trends
└── logs/                SLURM output
```
**Measured contents as of 2026-08-29 are in §13.1** — that is the live
inventory; this is the schema.

### Changes made 2026-08-27

**1. Domains are named, not ad-hoc boxes.** `download_plan.DOMAINS` defines
`north_atlantic` and `global`; `1_download_hpc.py` takes `--domain`.

**2. FIXED — a silent data-corruption bug.** `month_filename()` did not encode
the domain, so global and North Atlantic files for the same month shared a name
(`era5_2000-01.grib`) and **one would overwrite the other**. The integrity check
could not tell them apart: both have exactly 7 variables, 3 levels and
8×n_days timesteps. Only the grid differs.

Two fixes, both needed:
- filenames now carry a domain code: `era5_na_2000-01.grib`, `era5_glob_2000-01.grib`
- `check_grid()` added, verifying the horizontal grid against `expected_grid(area)`
  (NA 121×301, global 721×{1440,1441}).

**3. NEW `calibration.py` — thresholds are a file now.** Provides
`compute_thresholds`, `save_thresholds` / `load_thresholds` (JSON,
schema-versioned, with provenance: domain, period, levels, sample sizes, rojak
rev), and `compare_to_published`, which reports a ratio against Williams &
Joshi Table 1 and labels it, because there are **three** outcomes and only one
is good:

| ratio vs published | verdict |
|---|---|
| 0.1 – 10, same sign | `consistent` — the calibration works |
| exactly 1.0000 | `SUSPICIOUSLY EXACT` — and `aggregate.assert_thresholds_not_hardcoded_table2` will refuse to run |
| orders apart, or opposite sign | `ORDERS OFF` / `SIGN MISMATCH` |

### Deliberately NOT changed yet
The numbered module names (`2_diagnostics.py`) cannot be imported normally,
which is what caused the `sys.modules` crash. Renaming into a package would
end that class of bug permanently — **but not before phase 4 runs.** The
verification records, job scripts and this document all reference the current
names.

### Open sizing question — ANSWERED
A full global month is 10.6 GB and 5,208 fields in one CDS request, and the
Copernicus per-request limit was unconfirmed. Sidestepped rather than tested:
the calibration uses 4 days per month, **672 fields per request**, comfortably
inside any plausible limit (§11.8).

---

## 10. The ADA environment (verified on-cluster 2026-07-06, extended through 08-29)

Everything here was verified live, not assumed. This section is the answer to
"how do we run this".

### 10.1 Access
- ADA is the **renamed BAZIS**; old `bazis.readthedocs.io` links redirect to
  `rdm.vu.nl/tools/ada/`. Support: `itvo.ucit@vu.nl`.
- Account `yen230` (uid 722622). **Primary group is the Koopman group**
  (`R_SciStor-SBE-EDS-ClimateKoopman_C`), so new files are group-owned by the
  project automatically.
- Route in: **MobaXterm + EduVPN**, jump host `ssh.data.vu.nl` →
  `ada.labs.vu.nl`. **Two different passwords** — VUnetID at the gateway, a
  separate ITVO cluster password at ADA.
- Lands on `login01`/`login02`/`login03`. **Login nodes auto-log-off every 5
  minutes** if used for anything heavy — they are for editing and submitting
  only. Interactive work goes on `inter01`–`inter04`.
- **The repo is at `$BASE/data_turbulence`, the data at `$BASE`.** `git` run
  from `$BASE` fails with *"not a git repository ... Stopping at filesystem
  boundary"*. That is a wrong `cwd`, not a broken repo.

### 10.2 Modules
```bash
module load 2025 && module load pixi     # lowercase 'pixi'
```
The published docs say `module load Pixi` (capital P) — **that fails**.
Modules do not persist across sessions; batch scripts must re-load their own.

### 10.3 Scheduler

| Partition | Time limit | Nodes |
|---|---|---|
| `defq` *(default, use this)* | **7 days** | 17 |
| `defq-thin` | 7 days | 4 |
| `defq-fat` | 7 days | 2 |
| `defq-gpu` | 7 days | 6 |
| `bw` | infinite | 3 (restricted) |

**`defq` is heterogeneous** — this is the fact that governs scheduling:

| Memory | Nodes |
|---|---|
| 59 GB | 003, 004, 220–222 |
| 124 GB | 007, 008, 016, 017 |
| 252 GB | 009–011, 013–015 |
| 1031 GB | 001, 002 |

The fat and GPU nodes are *also* in `defq`. "17 nodes" is true but useless
without knowing that five of them cannot hold a 60 GB job.

**The 7-day cap is a design constraint, not a detail.** 504 monthly CDS
requests cannot run as one job. `1_download_hpc.py --index` exists for exactly
this — one array task per month.

### 10.4 Software — the reference environment
Pixi 0.46.0, **Python 3.12.14**, ecCodes v2.48.0, `import rojak` succeeds.

**rojak is pinned to rev `25b8685`** (`rojak_cat 1.0.2.dev20+g25b8685c6`,
2026-08-17). See `pixi.toml` in the repo root.

> **Pin moved forward from `1a65326` on 2026-08-27, deliberately.** `25b8685`
> is *"FIX: Geospatial laplacian was missing an extra derivative"* (#234). At
> `1a65326`, rojak's `spatial_laplacian` returned `df/dx + df/dy` — not a
> Laplacian at all. Upstream now routes it through `divergence()`, which uses
> `vector_derivatives()` and therefore carries the spherical curvature term.
> **That is the same correction this project derived independently from the
> divergence theorem the day before (§4g).** Two independent routes, one answer.
>
> **Nothing regresses.** The suite was run at both revs: identical error
> magnitudes, all analytic diagnostics unchanged, all 5 cross-checks at
> ρ = 1.0000.
>
> **The cross-check got stronger.** Because upstream's Laplacian is now
> correct, `4_verify.py` no longer monkey-patches it. UBF still matches at
> ρ = 1.0000 — but now against rojak's *own* implementation rather than
> against a second copy of our formula.
>
> **Pin a SHA, never a branch.** The pin is the reproducibility record.

**Six of the seven Phase-A findings still stand at `25b8685`**, so all seven
hand-written diagnostics remain necessary:

| Phase-A finding | Status at `25b8685` |
|---|---|
| `spatial_laplacian` returned gradient components | **FIXED upstream** ✅ |
| Richardson computes `N²/Sv` (must be `N²/Sv²`) | still broken |
| UBF returns `mass + inertial` (must be a residual) | still broken |
| Frontogenesis2D cross-term uses `dv_dy` (must be `dv_dx`) | still broken |
| UBF β via `latitudinal_derivative` | still broken — **docstring contradicts the code**: it states `β = 2Ω cosφ/R_E`, the body returns `coriolis_param / EARTH_AVG_RADIUS`, i.e. `f/R`. Worth reporting upstream. |
| NCSU1 / Colson–Panofsky inherit the broken Ri | still broken (follows from Ri) |
| RVA absent from rojak | still absent |

**numpy is not a free choice.** It is transitively bounded by
`rojak → sparse → numba`. Left to float, numba lags, the resolver falls back to
`sparse 0.3.1` (which needs Python <3.10), and the error message blames Python
rather than the resolver. Reference set: numpy 2.4.6 / numba 0.65.1 /
sparse 0.19.0.

### 10.5 Storage — SciStor, no auto-purge, no `quota` command

| Use | Path | Size | Used |
|---|---|---|---|
| Code, scripts, Pixi env | `/scistor/guest/yen230` (home) | 200 GB | ~1.5 GB |
| **ERA5 data + outputs** | **`/scistor/SBE-EDS-ClimateKoopman/yen230/`** | **2.5 TB** | **~68 GB** (§13.1) |

There is a `.snapshot` directory (point-in-time copies) — **confirm retention
and whether real backups exist** before treating derived output as safe. Raw
ERA5 is always re-downloadable, so it is lower stakes.

### 10.6 Sizing — measured 2026-08-27, superseded by §13.2

| Quantity | Value |
|---|---|
| One full day, 7 vars × 3 levels × 8 steps, NA box | **12.26 MB** (2.00 bytes/gridpoint — 16-bit GRIB packing) |
| Average month | **373 MB** |
| **Raw GRIB, 1979–2020** (15,341 days) | **188 GB** |
| Timesteps in the full run | 122,728 *(matches `chunk_stitch.py`)* |

The derived-side figures in the original version of this table (1.13 TB, 0.38 TB)
assumed dual-format or float64 output. **They are superseded by §13.2**, which
measures 27 real months written by the production driver.

### 10.7 Network
- **Login node egress works.** conda-forge, GitHub and Copernicus all
  reachable; `curl -sSI https://cds.climate.copernicus.eu/api/` → HTTP 308.
- **cdsapi auth works from ADA.** `~/.cdsapirc` present, `chmod 600`, new
  format (`url: https://cds.climate.copernicus.eu/api`, key = Personal Access
  Token).
- **✅ COMPUTE-NODE EGRESS CONFIRMED (2026-08-27).** Tested with `srun` on
  `node222`: CDS → HTTP 202, PyPI → HTTP 200, no proxy vars set. Downloading
  happens in a `defq` **job array**, one month per task.
- **⚠ `--ntasks=1` IS MANDATORY.** Every line of that egress test printed
  **twice** — `srun` allocated two tasks and ran the command twice on one node.
  In a download job that means two processes requesting the same month from CDS
  and racing on the same `<name>.tmp`. The atomic-rename guard in
  `_fetch_verified` protects against a *half-written* file, **not against two
  concurrent writers**. Every job script sets `--ntasks=1` explicitly.

### 10.8 cfgrib multi-variable open — RESOLVED 2026-08-27
```
xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
  -> SUCCESS, single open, no filter_by_keys
  data_vars: ['u','v','t','z','d','vo','pv']
  dims     : {'time': 8, 'isobaricInhPa': 3, 'latitude': 121, 'longitude': 301}
  domain   : lat 30.00..60.00, lon -75.00..0.00, 0.25 deg, levels [225,200,175]
```
Matches `download_plan.NORTH_ATLANTIC_BOX` exactly. `verify_file` passes with
zero problems. Risk closed.

### 10.10 Environment built and verified on ADA — 2026-08-27

Job `1091550` on `node009`, exit 0. `pixi install` from the repo's `pixi.toml`,
1.2 GB, rojak pinned at `25b8685`. **The full verification suite passes on
ADA** — that is the point: the accuracy figures in `verification/` are
reproducible in the environment that will produce the dataset.

**Three bugs found getting there** (jobs 1091491, 1091547, 1091549 in §13.3):

1. **`pixi run python -c "..."` is not safe for multi-line Python.** `pixi run`
   parses the string with its own task-shell parser before Python sees it, and
   that parser rejects English words that happen to be shell reserved words.
   The first run died on the word **"in"**, inside a Python string, inside an
   assertion message — with a perfectly healthy environment. All checks now
   live in `ada/verify_env.py` as a real file.
2. **rojak has undeclared imports.** `dask-geopandas`, `requests` and `rich`
   are imported at module load but not declared. Now named in `pixi.toml`.
3. **`cdsapi.Client()` defaults to `retry_max=500, sleep_max=120`** — over
   **sixteen hours** of silent retrying on an unreachable CDS. Now bounded to 8
   attempts backing off to a minute, then a loud failure — the right shape
   precisely *because* the downloader is resumable: a failed month leaves no
   final file, so re-submitting the array retries exactly the months that
   failed. **Validated in production on 2026-08-29**, §13.3.

Also: `check_authentication()` does not exist across all cdsapi builds. The
probe tries a couple of locations and degrades to a note; a missing or
old-format `~/.cdsapirc` is still fatal.

### 10.11 First real month downloaded — 2026-08-27

Job `1091553_0` on `node009`. `era5_na_1979-01.grib`, 363 MB, integrity check
PASSED, 6 min 55 s.

### 10.12 First month of diagnostics — 2026-08-28

Job `1091555` on `node013`, exit 0. `3_pipeline.py` on `era5_na_1979-01.grib`.
All 21 computed, no failures, no NaN placeholders. 8 min 11 s, 14.8 GB peak,
1.3 GB NetCDF + 1.1 GB Zarr.

> **Three corrections to this section, made 2026-08-29:**
>
> 1. **The 14.8 GB MaxRSS figure is not reproducible by the method that
>    produced it** — job-step accounting is disabled cluster-wide (§11.9), so
>    `sacct`/`sstat` return empty MaxRSS. The global run's independently
>    measured 17.7 GB on a comparable point count corroborates it.
> 2. **This run wrote both NetCDF and Zarr.** Production writes zarr only; the
>    resulting storage plan is §13.2, not the 796 GB quoted here. Its output
>    directory `derived/north_atlantic/1979-01/` was 2.3 GB of dead weight and
>    was **deleted on 2026-08-29**.
> 3. **`colson_panofsky` reading `1-2 ORDERS OFF` was a units bug in
>    `SCALE_TO_TABLE`, not a property of the diagnostic** (§11.5).

**The output is 200 hPa only.** rojak computes on `(121, 301, 248, 3)` and
`compute_all_21(target_level=200)` saves `(121, 301, 248)`. That is correct —
three levels give the vertical-derivative stencil, the middle one is where
diagnostics are evaluated — but it is an implicit decision. **If the
econometrics ever wants 175 or 225 hPa, that is a re-run of all 504 months.**

**Comparison against Williams & Joshi Table 1: read it as season, not error.**
January in the North Atlantic at 200 hPa *is* the jet stream; W&J's medians are
a far broader climatology. Wind speed **29.8 vs 14.9 m/s** is the giveaway —
30 m/s is a textbook winter jet core. Everything downstream inherits that
(VWS 4.67 vs 1.88, DEF 43.8 vs 18.6, TI1 184 vs 31.5). Several land close
regardless: brown1 80.6 vs 77.1, PV 5.19 vs 6.84, NCSU1 8.98 vs 11.1.

### 10.9 Still open on ADA
- Snapshot retention / backup policy on the Koopman share.
- Node-local scratch on compute nodes (`$TMPDIR`) — minor.
- **The `unlimited` QOS question** (§11.9) — the one item that changes phase 4.
- `1_download.py` has **two trial dates in circulation** — the uploaded copy
  uses 2016-01-01, the project-folder copy 2016-07-01. Pick one deliberately.

---

## 11. Stage 2 executed — the calibration check, 2026-08-28/29

**The headline: the chain reproduces Prosser (2023) Table 1 to within 24 %,
and the identity check on real data is exact.** The production dataset does not
exist yet, but the method that will produce it has been validated end to end
against a published result.

### 11.1 What was run

Full ledger with job IDs in §13.3. In summary: **14 CDS requests and about four
hours of compute**, against 504 requests and days for production. All 15 input
files passed integrity checks (`ada/check_calib_downloads.py`); all 15
diagnostics runs produced 21 of 21 with no failures and no all-NaN placeholders.

### 11.2 The result

North Atlantic DJF 1979, in **Prosser's own box (36–60°N / 55–10°W)**:

| | ours | Prosser | ratio |
|---|---|---|---|
| LOG | 5.037 % | 5.968 % | 0.84 |
| LMOG | 1.684 % | 2.111 % | 0.80 |
| **MOG** | **0.785 %** | **1.032 %** | **0.76** |
| MSOG | 0.399 % | 0.560 % | 0.71 |
| SOG | 0.200 % | 0.296 % | 0.68 |

Order of magnitude **PASS**. Ladder shape, worst deviation 11 % — **PASS**.
Absolute agreement 24 % — **PASS**, and this one was not expected to be met.

**Caveat that §12 later supplied:** this compares a single *raw* DJF 1979
against Prosser's *fitted* 1979. Some of the 24 % is 1979 being one draw from a
distribution whose interannual spread is large (§12.4). Against a nine-season
fitted 1979 the MOG ratio is 0.79 rather than 0.76.

**The identity check came back exact**: applying the year-2000 thresholds back
to the year-2000 data returns 3.0000 / 0.9000 / 0.4000 / 0.2000 / 0.1000 %,
0.00 % relative error on all five. That is what licenses reading the table
above; without it, agreement with Prosser would be coincidence.

### 11.3 The one systematic pattern, and its cause

The ratios decline monotonically (0.84 → 0.68): our tail is thinner than
Prosser's, more so the further out you go.

**The vertical stencil predicts exactly this.** We differentiate across
175→225 hPa (50 hPa); Prosser uses 188→206 (18 hPa). A 2.8× wider stencil
smooths across sharp shear layers and under-resolves shear *maxima* — barely
visible in the bulk, progressively damaging in the extreme tail.

This is §5's documented level substitution, finally measured. Per §11.6 it
cannot be closed, so it should be quoted as a known property of the dataset
rather than an open defect. §12.5 shows the same parameter from the other side,
and is the stronger evidence for this reading.

### 11.4 Method facts now pinned

- **The percentile ladder is LOCKED**: LOG p97.0, LMOG p99.1, MOG p99.6,
  MSOG p99.8, SOG p99.9. Williams (2017) Table 1, adopted verbatim by Prosser.
  Origin: a log-normal EDR distribution constrained by an observed 3.0 % LOG
  and 0.4 % MOG probability.
- **Prosser recomputes thresholds from ERA5; he does not reuse Williams'
  numbers.** This was genuinely open and is now settled from the paper's own
  text. Latitude weighting is required by the method, not an embellishment.
- **Prosser's North Atlantic box is 36–60°N / 55–10°W** (his Figure 2 caption).
  It is a strict subset of our 30–60°N / 75–0°W download, so Table 1 can be
  compared exactly. **This coordinate pair was recorded nowhere before.**
  Comparing on our larger box gives systematically different frequencies and
  looks like an error.
- **Prosser publishes no threshold values.** Confirmed by reading his
  Supporting Information directly: Figures S1–S5, no tables. The only published
  per-diagnostic threshold table anywhere is Williams (2017) Table 2, and it is
  GFDL-CM2.1 — magnitudes are not comparable, as Williams states himself.

### 11.5 Bugs found and fixed

1. **Colson–Panofsky was compared in the wrong units.**
   `3_pipeline.SCALE_TO_TABLE["colson_panofsky"]` was `1.0` with a comment
   promising a "see note below" that was never written, so a native value in
   m² s⁻² was compared against a table printed in **10³ kt²**. That single line
   is the entire reason CP has read `1-2 ORDERS OFF` since the comparison table
   was built. Correct factor 1/(1000 × 0.2646526) = **3.778538 × 10⁻³**.

2. **`calibration.compute_thresholds` sorted each array five times.** It asked
   `weighted_percentile` for one severity at a time, and that function sorts its
   whole input on every call — 105 full sorts of 4 × 10⁸ values where 21 would
   do. **The first calibration attempt (job 1092585) hit walltime at 1 h and was
   killed because of this.** Now fixed in `calibration.py` itself: one
   `np.atleast_1d` call with the five percentiles as an array, with a length
   assertion. Verified bit-identical, **5.35× faster**, 31 round-trip tests pass.

3. **A unit error in the new comparison code itself** (mine): Williams prints
   NVA/RVA in 10⁻⁹ s⁻² where `REFERENCE_TABLE` uses 10⁻¹⁰; converting requires
   multiplying by ten and the first version divided, inflating those two ratios
   by a hundred and making them look like wild outliers.

### 11.6 The 188/197/206 levels are MODEL levels — the branch is closed

Checked against ECMWF's L137 definitions: **188.29, 197.37 and 206.81 hPa are
model levels 73, 74 and 75.** Prosser used ERA5 model-level data described by
nominal pressures, not the pressure-level product.

- They **cannot** be requested from `reanalysis-era5-pressure-levels`, whose
  fixed 37-level set offers 150, 175, 200, 225, 250 around 200 hPa and nothing
  between.
- Model-level ERA5 is on the **MARS tape archive** (`reanalysis-era5-complete`)
  — different request format, slow access.
- **Therefore 175/200/225 is not a compromise: it is the finest stencil
  available at 200 hPa on CDS.** Phase 3's "accept 175/200/225, or pursue model
  levels" now has a price attached, and the answer is accept.

Model levels are also **hybrid sigma-pressure**, i.e. terrain-following, so the
actual pressure at level 74 varies with location and time. Even with MARS
access the comparison would not be exact.

**What can be done instead** (§11.1 of `CALIBRATION_REFERENCE.md`): request
150/175/200/225/250 in one download and compute at 200 hPa with both a 50 hPa
and a 100 hPa stencil. Two points on the stencil-width curve let the 18 hPa
case be *extrapolated* rather than measured. ~1 day. Optional but valuable,
because "why are your levels different from Prosser's" is a certain referee
question and this converts the answer from an apology into a number. **The
production configuration cannot change either way** — there is no tighter
pressure-level stencil to move to — so this does not gate phase 4.

### 11.7 ERA5 vs ERA5.1 — also MARS-only

Prosser used **ERA5.1** for 2000–2006, which corrects a lower-stratospheric
cold bias — and year 2000 is our calibration year. ERA5.1 is **not a CDS disk
dataset**; it lives in MARS. ECMWF's own note is that behaviour "in most of the
troposphere is similar to that in ERA5.1". Proceeding with ERA5; revisit at
publication, not before.

### 11.8 The sub-sampled calibration design

Days 1/9/17/25 of each month, all 8 three-hourly steps: **48 days,
3.99 × 10⁸ points per diagnostic, 1.59 × 10⁶ above the MOG threshold.**

- Split-half (days 1+17 vs 9+25): median disagreement **0.81 %**, 90th
  percentile **5.90 %**. **48 days is enough.** The script's "marginal" verdict
  came from a 92 % worst case on `colson_panofsky @ LMOG`, which is an artifact
  — CP's LMOG threshold sits at −1.04, essentially zero, and a relative
  difference across a sign change is meaningless. The metric should divide by
  the diagnostic's spread, not its value.
- It also **made stage 2 possible at all**: a full global year is 3.04 × 10⁹
  points per diagnostic, and `compute_thresholds` ravels the whole field —
  over 100 GB peak. It does not fit any ADA partition.
- And it sidestepped the unconfirmed CDS per-request field limit: 672 fields
  per request instead of 5,208.

### 11.9 ADA facts learned the hard way

- **Memory footprint is the biggest lever on time-to-start.** A 120 GB request
  was scheduled **22 hours out** while 24 GB started in 13 seconds. Do not
  round requests up "to be safe" — it costs a day. Jobs 1092472 (250 GB) and
  1092477/78/82 (120 GB) were all cancelled for this before 1092495 ran at
  24 GB (§13.3).
- **Do not submit to `defq-thin` or `defq-fat` directly.** Their nodes are
  shared with `defq`, which has higher priority, so a job addressed to the
  narrower partition queues behind everything. Ask `defq` for the memory you
  need and let the scheduler pick.
- **Job-step accounting is DISABLED.** `sacct` and `sstat` both return empty
  MaxRSS. Every driver must measure and print its own peak RSS
  (`resource.getrusage`).
- **`export PYTHONUNBUFFERED=1` in every job script.** Python block-buffers
  stdout to a file, so a working job is indistinguishable from a hung one.
  This cost an hour of watching a blank log.
- **`mkdir -p logs` before `sbatch`, not inside the script.** SLURM opens the
  output file before the script runs.
- **Short, honest walltimes get backfilled.** 6 h waits; 1 h runs.
- **HARD CAP: 8 CONCURRENTLY RUNNING JOBS PER USER.** The default `normal`
  QOS sets `MaxJobsPU = 8`; confirmed in `sacctmgr show qos` and in the
  association manager. Array tasks count as jobs, so a 9th sits as
  `QOSMaxJobsPerUserLimit`. **This makes any `%` throttle above 8 a no-op.**

  **There is an `unlimited` QOS on this cluster** with no `MaxJobsPU`, and an
  `ood` QOS capped at 2. Whether `yen230` may use `unlimited` is a question for
  `itvo.ucit@vu.nl`, and it is the single highest-leverage thing to ask before
  phase 4. Check what is currently granted with:

      sacctmgr show assoc user=yen230 format=User,Account,Partition,QOS,DefaultQOS

- **CDS REJECTS AT 8 CONCURRENT REQUESTS PER USER. `%4` IS VALIDATED.**
  Measured, not guessed:

  | Concurrent CDS requests | Outcome |
  |---|---|
  | 4 (global calib at `%4`) | 12/12 succeeded |
  | 6 (that plus 2 NA) | all succeeded |
  | **8** (no throttle, QOS-capped) | **8 of 12 failed** (job 1092924) |
  | 4 (the same 12 resubmitted) | all succeeded (job 1092941) |
  | 4 (12 fresh months) | **12/12 succeeded** (job 1092983) |

  The failures are `400 Client Error: Bad Request` on
  `/retrieve/v1/jobs/<id>/results` — Copernicus refusing, not a transient
  gateway blip like the 502s that recover on attempt 1 of 8.

  **The 504-month production download would have lost roughly two thirds of its
  months this way**, and only discovered it at the integrity check. `%4` is not
  politeness, it is the working limit.

- **⚠ `sbatch --array=...` ON THE COMMAND LINE DISCARDS THE `%N` IN THE
  SCRIPT.** `01_download.sbatch` declares `--array=0-503%4`; overriding with
  `--array=132,133,...` replaced the *entire* spec, throttle included, and the
  array ran at whatever the QOS allowed. That is what caused the failures above.
  Any override must re-append it:

      sbatch --array=72,73,83,192,193,203,312,313,323,432,433,443%4 \
             jobs/01_download.sbatch

  This fails **silently** — you get the array you asked for and no warning that
  the concurrency limit went with it. Confirmed twice: the botched submission
  (8/12 lost) and the corrected one (12/12 landed).

- **THE RESUMABLE DESIGN WORKS UNDER FAILURE, MEASURED.** When job 1092941
  re-ran the full twelve-month index list, the four months that had already
  succeeded detected their own output and exited in **9 seconds each**; the
  eight real re-downloads took 7–23 min. Same on the diagnostics side: job
  1093047 ran 27 tasks, of which 15 skipped in ≤1 s. **Recovery from a partial
  failure is "resubmit the same array", with no hand-maintained list of which
  months failed.** For 504 months, where some CDS failures are certain, this is
  the property that makes phase 4 manageable.

### 11.11 New files, 2026-08-28/29

```
ada/calibration_check.py        stage 2->3: calibrate, verify, compare to Prosser
ada/diagnostics_global.py       21 diagnostics, overlap-buffered time chunking
ada/check_calib_downloads.py    verify all 15 inputs in one pass
ada/trend_check.py              the trend comparison, with slope uncertainty
ada/progress.sh                 one line per job; newest array by default
jobs/02b_download_global_calib.sbatch
jobs/04_diagnostics_global.sbatch
jobs/05_diagnostics_na_djf.sbatch
jobs/06_calibration_check.sbatch
jobs/07_diagnostics_na_trend.sbatch
jobs/08_trend_check.sbatch
tests/test_calibration_roundtrip.py   31 tests, the identity, no ERA5 needed
RUNBOOK_calibration_check.md
calibration/thresholds_2026-08-29.json
```

`--days` added to `1_download_hpc.py`, with a distinct filename
(`era5_glob_2000-01_d01-09-17-25.grib`) so a partial month can never be
mistaken for a complete one.

**`ada/diagnostics_global.py`'s chunking is bit-identical to whole-file
processing** — verified across 1-, 2-, 4-day and whole-file chunk sizes. Each
chunk is computed with one extra timestep on each side and trimmed back, so
F2D loses only the 2 steps at a file's true ends. Naive chunking would cost it
8. Same device as `chunk_stitch.py`, one level down.

**`ada/progress.sh` had a defaults bug, fixed 2026-08-29.** Its default filter
was the substring `diag-`, which matches every diagnostics log ever written.
Bash expands globs lexicographically, so the oldest finished run printed first
and the array you had just submitted printed last, below the fold — it looked
stuck on a completed job. The default now resolves the newest log by mtime,
takes its array id, and shows only that array; `--all` restores the old
behaviour.

### 11.12 Three independent literature anchors

The calibration check compares against Prosser. Three further published tables
were checked so that a single agreement is not carrying the whole argument:

| Source | Data | What it establishes | Result |
|---|---|---|---|
| **Williams (2017) Table 2** | GFDL-CM2.1, ~2° | sign and units only — magnitudes not comparable, as Williams says himself | consistent |
| **Lee et al. (2023) Table 1** | **ERA5** | the closest match in data source | **all five diagnostics within ±23 %** |
| **Sharman (2006) Table B1** | 20-km RUC | the only table covering the hand-written seven | **six of seven within a factor of 3** |

Sharman Table B1 in detail: TI1 1.02, |∇ₕT| 1.03, UBF 1.85, CP 2.20,
NCSU1 2.75, −Ri 0.60, F2D 35.6. **This cleared UBF and NCSU1**, which had
looked alarming against the GCM and are simply resolution-sensitive. F2D's
35.6× is consistent with the §4d diurnal-sampling damping plus a 20-km-versus-
0.25° resolution gap. `brown2` is absent from B1, which is why §7 still lists
it as a suspect.

Full tables and reasoning in `CALIBRATION_REFERENCE.md`.

---

## 12. Stage 2b — the trend check, 2026-08-29

The calibration check verified a **level** in one season. Prosser's actual
claim, and the reason this project exists, is a **change**. Layers 5 and 6
(annual aggregation, per-gridpoint regression) had never touched real data.

`ada/trend_check.py` fits DJF exceedance across several seasons and compares
the fitted 1979 → 2020 change against Prosser Table 1's DJF row. Run twice:
once at five seasons, then at nine after §12.8 showed five could not decide
anything.

### 12.1 Nine seasons — the per-season table

DJF 1979, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, in Prosser's box:

| DJF | LOG | LMOG | MOG | MSOG | SOG |
|---|---|---|---|---|---|
| 1979 | 5.037 % | 1.684 % | 0.785 % | 0.399 % | 0.200 % |
| 1985 | 5.852 % | 2.053 % | 0.991 % | 0.519 % | 0.267 % |
| 1990 | 5.990 % | 2.036 % | 0.965 % | 0.501 % | 0.253 % |
| 1995 | 5.597 % | 1.958 % | 0.951 % | 0.505 % | 0.264 % |
| 2000 | 5.648 % | 1.960 % | 0.945 % | 0.494 % | 0.252 % |
| 2005 | 4.788 % | 1.634 % | 0.801 % | 0.435 % | 0.233 % |
| 2010 | 5.339 % | 1.846 % | 0.878 % | 0.458 % | 0.234 % |
| 2015 | 6.884 % | 2.545 % | 1.293 % | 0.714 % | 0.386 % |
| 2020 | 6.870 % | 2.531 % | 1.285 % | 0.713 % | 0.391 % |

### 12.2 The fitted trend against Prosser

| | our 1979 | our 2020 | ours % | 95 % CI | Prosser % | ratio | R² | t |
|---|---|---|---|---|---|---|---|---|
| LOG | 5.205 % | 6.331 % | 22 % | [−10 %, 53 %] | 21 % | 1.03 | 0.27 | 1.63 |
| LMOG | 1.741 % | 2.304 % | 32 % | [−7 %, 72 %] | 30 % | 1.08 | 0.35 | 1.93 |
| **MOG** | 0.811 % | 1.159 % | **43 %** | **[−3 %, 89 %]** | **37 %** | **1.16** | 0.41 | 2.20 |
| MSOG | 0.411 % | 0.638 % | 55 % | [3 %, 108 %] | 43 % | 1.29 | 0.47 | 2.48 |
| SOG | 0.203 % | 0.345 % | 70 % | [9 %, 131 %] | 49 % | 1.43 | 0.51 | 2.72 |

In hours per season at an average point in the box:

| | ours | Prosser |
|---|---|---|
| LOG | 112.4 → 138.3 | 128.9 → 155.6 |
| LMOG | 37.6 → 50.3 | 45.6 → 59.3 |
| MOG | 17.5 → 25.3 | 22.3 → 30.6 |
| MSOG | 8.9 → 13.9 | 12.1 → 17.2 |
| SOG | 4.4 → 7.5 | 6.4 → 9.6 |

**The four criteria:**

| | result |
|---|---|
| 1. All trends positive | **PASS** |
| 2. Trends distinguishable from zero at 5 % (n=9, df=7, needs \|t\|>2.37) | **2/5 — WEAK** |
| 3. **Prosser inside our 95 % interval** | **5/5 — PASS** |
| 4. Stronger at higher severity | **PASS** (Prosser: 21/30/37/43/49 %) |

### 12.3 How to read that — the honest version

**Criterion 3 is the result**, and it is the one that could have failed:
Prosser's published change lies inside our interval at every severity, with
point estimates tracking his closely (ratios 1.03 → 1.43).

**But criterion 3 passes partly because the intervals are wide, and criterion 2
is the counterweight.** MOG's interval, [−3 %, 89 %], contains Prosser's 37 %
and also contains zero. Only MSOG and SOG clear significance; MOG sits at
t = 2.20, p ≈ 0.064. So what has been established is **consistency with
Prosser, not independent confirmation of a trend.** Anyone reporting this
should say so in those words.

That is not a defect in the pipeline. It is what §12.4 says it is.

### 12.4 The finding that matters: variability exceeds the signal

R² *fell* when seasons were added (MOG 0.59 at n=5 → 0.41 at n=9), because the
new seasons brought scatter rather than confirming a line. The series is not
monotone: 2005 dips to 0.801 %, below several earlier seasons, while 2015 and
2020 are both high at ~1.29 %.

**Across nine winters the MOG frequency spans 0.785 % to 1.293 % — a factor of
1.65 — against a fitted 41-year change of 43 %.** Interannual variability is
larger than the entire trend signal, and the fit leans heavily on two high
winters at the end of the sample. Both 2015 and 2020 were strongly positive-NAO
DJFs, and NAO phase modulates the North Atlantic jet directly, so a nine-season
sample taken every five years can easily land on a biased draw.

This is exactly why Prosser fits **42 consecutive years** and describes his
Table 1 as a guide to the statistics "in the absence of interannual
variability". It is recorded as risk 4 in §5.

**Consequence: no subsample can settle the trend question, and the check was
never going to.** What the trend check was for — showing that Layers 5–6 run on
real data and produce trends of the right sign, magnitude and severity ordering
— it did. Resolving the trend itself is what phase 4 is for.

### 12.5 Level deficit and trend excess are one fact, not two

| | level ratio | trend ratio | product |
|---|---|---|---|
| LOG | 0.87 | 1.03 | 0.90 |
| LMOG | 0.82 | 1.08 | 0.89 |
| MOG | 0.79 | 1.16 | 0.91 |
| MSOG | 0.73 | 1.29 | 0.95 |
| SOG | 0.69 | 1.43 | 0.98 |

Mirror images, both monotone in severity, products between 0.89 and 0.98.

The wide stencil under-resolves shear maxima, so we sit lower on the tail of
the EDR distribution; lower on the tail means the *same* physical
intensification produces a *larger* relative change, because exceedance
probability is convex in the shift. **One parameter seen twice.** A clean
multiplicative offset would instead have left the trend ratio at 1.00 at every
severity; it does not, which is why §5 has been amended.

Two consequences worth stating in any write-up: the cancellation is real but
not exact (products drift from 0.89 to 0.98 rather than sitting flat), and
because the two effects partially offset, **our agreement with Prosser improves
across the record** — the MOG level ratio goes from 0.79 in 1979 to 0.83 in
2020.

### 12.6 Per-diagnostic breakdown — and the one anomaly

MOG exceedance by diagnostic, DJF 1979 → DJF 2020:

| diagnostic | 1979 | 2020 | change |
|---|---|---|---|
| `vertical_wind_shear` | 2.321 % | 3.462 % | 49 % |
| `brown2` | 1.993 % | 3.029 % | 52 % |
| `ti1` | 1.639 % | 2.476 % | 51 % |
| `ti2` | 1.566 % | 2.288 % | 46 % |
| `temperature_gradient` | 1.349 % | 2.270 % | 68 % |
| `ngm1` | 1.234 % | 2.251 % | 82 % |
| `rva_magnitude` | 0.857 % | 1.636 % | 91 % |
| `nva` | 0.827 % | 1.460 % | 77 % |
| `endlich` | 0.785 % | 1.528 % | 95 % |
| `wind_speed` | 0.713 % | 1.462 % | 105 % |
| `deformation` | 0.642 % | 1.014 % | 58 % |
| `ncsu1` | 0.610 % | 0.937 % | 54 % |
| `brown1` | 0.592 % | 0.986 % | 67 % |
| `horizontal_divergence` | 0.572 % | 0.792 % | 38 % |
| **`ubf`** | **0.249 %** | **0.276 %** | **11 %** |
| `vorticity_squared` | 0.236 % | 0.470 % | 99 % |
| `ngm2` | 0.115 % | 0.180 % | 56 % |
| `magnitude_pv` | 0.075 % | 0.175 % | 132 % |
| `colson_panofsky` | 0.056 % | 0.138 % | 148 % |
| `f2d` | 0.026 % | 0.086 % | 226 % |
| `negative_richardson` | 0.022 % | 0.061 % | 176 % |

Spread at DJF 1979: min 0.022 %, median 0.642 %, max 2.321 % — **105×**. That
range is expected, not alarming: the thresholds are global percentiles evaluated
on a regional winter box, and the 21 diagnostics measure different things.

**The table is otherwise strikingly coherent.** Trend magnitude
anti-correlates with exceedance level — the four sparsest diagnostics
(`negative_richardson`, `f2d`, `colson_panofsky`, `magnitude_pv`) carry the four
largest relative trends. That is the same tail-convexity effect as §12.5,
appearing independently across diagnostics rather than across severities. Two
different cuts of the data, one mechanism.

**`ubf` is the exception and is now the top item in §7.** It is 15th of 21 on
level but last on trend, at +11 % where every sibling is +38 % to +226 %. It
breaks the pattern that explains all the others. It is also the diagnostic whose
geometry was rewritten on 26 August (§4g) and a residual of near-cancelling
large terms — the most fragile of the 21 by construction. Prosser's Figure S5 is
the published version of this breakdown and is the right thing to compare
against.

**`brown2` is partially reassured.** §7 flags it as five orders from its
published value, but here its level (2nd of 21) and trend (+52 %, mid-pack) are
both unremarkable. Its *ranks* look healthy even though its magnitude does not,
and ranks are all the exceedance counting uses.

### 12.7 What gates phase 4 — nothing

- The **level deficit cannot be fixed by downloading anything**: 175/200/225 is
  the finest pressure-level stencil that exists (§11.6). It is a bias to
  characterise and report, and §12.5 characterises it.
- The **trend replicates** in sign, magnitude and severity gradient, and
  contains Prosser's value at all five severities.
- The **trend cannot be resolved from a subsample** (§12.4), which is an
  argument for running the full 504 months rather than against it.

Remaining before submitting phase 4: the **`unlimited` QOS question** (§8), and
optionally the `ubf` look (§12.6) — which does not need new data, only the 27
months already on disk.

### 12.8 The five-season run, and a test that was retired

The first run (job 1092975) used 1979, 1990, 2000, 2010, 2020 and printed
**NOT confirmed** on a test framed as: *the trend ratio should be closer to
1.00 than the level ratio, because a constant offset cancels in a change.*
Trend was 26 % off, level 24 % off.

**That verdict carried no information, and the test was a mistake.** With n = 5
and R² = 0.59 on MOG the slope has 3 degrees of freedom: t = 2.1, p ≈ 0.13, and
a 95 % interval running roughly −25 % to +165 %. Prosser's +37 % sat comfortably
inside. An interval containing both the null and the target cannot distinguish
them, so the comparison could not have failed and therefore could not have
passed either.

`ada/trend_check.py` was rewritten the same day to report a standard error, a t
statistic and a 95 % interval on every fitted trend, and to make the verdict
**"does Prosser's value lie inside our interval"** rather than a comparison of
two point estimates. It also gained `--years` so the five-season run can be
reproduced exactly.

**The lesson is general and worth keeping:** a comparison of point estimates
from a small sample is not a test. Every future check in this project that
compares a fitted quantity against a published one should carry an interval, or
it is decoration.

---

## 13. Run ledger and measured disk inventory

### 13.1 What is on the share — measured 2026-08-29

| Path | Size | Contents |
|---|---|---|
| `raw/global/` | 16 GB | 12 GRIB months, year 2000, days 1/9/17/25 (+ cfgrib `.idx` sidecars) |
| `raw/north_atlantic/` | 9.3 GB | 27 GRIB months = 9 DJF seasons × 3 (+ sidecars) |
| `derived/global/` | 26 GB | 12 `.zarr`, ~2.17 GB each |
| `derived/north_atlantic/` | ~15 GB | **27 `.zarr`**, ~0.57 GB each |
| `calibration/` | 7 KB | `thresholds_2026-08-29.json` |
| **total** | **~68 GB** | **~2.7 % of the 2.5 TB share** |

`derived/north_atlantic/1979-01/` — 2.3 GB of NetCDF+Zarr from job 1091555,
the only dual-format artefact on the share — was **deleted on 2026-08-29**.

Directory entry counts exceed file counts (42 vs 27 in `raw/north_atlantic`)
because cfgrib writes `.idx` sidecars. Harmless.

### 13.2 Production sizing, from measured months

| Per month | measured |
|---|---|
| NA raw GRIB, 31-day month | 363 MB (328 MB for February) |
| NA derived zarr, 21 diagnostics, 200 hPa, float32 | **~0.57 GB** |
| Global raw, 4 days sub-sampled | ~1.33 GB |
| Global derived zarr | ~2.17 GB |

Projected to 504 months:

| | size |
|---|---|
| Raw GRIB, 1979–2020 | **~183 GB** (§10.6 predicted 188 GB) |
| Derived zarr, all 21, 200 hPa, float32 | **~287 GB** |
| **Production total** | **~470 GB — 19 % of the share** |

This supersedes §10.12's 796 GB and §10.6's 1.32 TB, both of which assumed
dual-format or float64 output.

**Corrected phase-4 timing** (the §10.12 plan to run diagnostics at `%16` and
finish in 4.3 h is unreachable — the QOS caps concurrency at 8, and shared
filesystem I/O contends: six concurrent tasks stretched a 30-minute job to 36,
mostly in GRIB decode):

| Stage | concurrency | 504 months |
|---|---|---|
| Download | `%4` (CDS limit, §11.9) | **~36 h** |
| Diagnostics | `%8` (QOS limit) | **~7 h** |

**The download figure was revised upward on 2026-08-29 and the earlier one was
wrong.** An estimate of 7 min/month was carried forward from job 1091553_0,
which ran at midnight against an idle CDS. Twenty months downloaded during
European daytime (jobs 1092941, 1092983) averaged **18.4 min**, range 5–27.
At `%4` that is 477 remaining months × 18.4 min ÷ 4 ≈ **36 hours**, i.e. a day
and a half, possibly two if CDS is busy.

Nothing about that is a problem — each task is one month, far inside the 7-day
per-task cap, and the array is resumable (§11.9), so it survives interruption.
But it is a day and a half, not the "8–16 h" this table previously claimed.

If `unlimited` is granted the diagnostics stage shrinks toward the I/O ceiling;
**the download stage will not**, because it is bounded by Copernicus, not by
ADA. The `%4` throttle is a hard external limit, so phase 4's wall clock is
essentially fixed.

### 13.3 Run ledger

Every job since the environment build. Elapsed is `sacct` wall clock. `ReqMem`
is what was *requested*, not used — step accounting is off (§11.9).

| Job | Name | Result | Elapsed | Mem | What it did |
|---|---|---|---|---|---|
| 1091491 | cat-setup | FAILED | 3:40 | 8G | `pixi run python -c` shell-parsing bug (§10.10 #1) |
| 1091547 | cat-setup | FAILED | 0:27 | 8G | rojak undeclared imports (§10.10 #2) |
| 1091549 | cat-setup | FAILED | 0:11 | 8G | same |
| **1091550** | cat-setup | **COMPLETED** | 1:04 | 8G | **environment built and verified on ADA** |
| 1091551 | cat-smoke | COMPLETED | 0:59 | 8G | 7 pre-flight checks, one real CDS request |
| 1091553_0 | era5-dl | COMPLETED | 7:10 | 4G | first real month, `era5_na_1979-01.grib` |
| 1091555 | cat-diag | COMPLETED | 8:11 | 64G | first diagnostics, via `3_pipeline.py` (dual-format) |
| 1092324_0 | era5-glob-calib | COMPLETED | 1:49 | 8G | single-task trial of the global download |
| 1092393_0–11 | era5-glob-calib | COMPLETED | 0:17–22:54 | 8G | **12 global months at `%4` — 12/12** |
| 1092394_1, _11 | era5-dl | COMPLETED | 12:15, 14:30 | 4G | NA Feb + Dec 1979 |
| 1092472 | cat-diag-glob | CANCELLED | — | **250G** | never scheduled — memory request too large |
| 1092477/78/82 | cat-diag-glob | CANCELLED | — | **120G** | queued 22 h out (§11.9) |
| 1092490_0 | cat-diag-glob | CANCELLED | 3:34 | 24G | superseded |
| 1092495_0 | cat-diag-glob | COMPLETED | 30:19 | 24G | global month 0 — the 24 GB request started in 13 s |
| 1092536_1–11 | cat-diag-glob | COMPLETED | 15:29–36:11 | 24G | global months 1–11 |
| 1092537_0–2 | cat-diag-na | COMPLETED | 5:08–8:22 | 24G | NA DJF 1979 |
| 1092585 | cat-calib | CANCELLED | 1:01:55 | 32G | **hit walltime — the 5× sort bug** (§11.5 #2) |
| **1092630** | cat-calib | **COMPLETED** | **1:57:28** | 32G | **the calibration check — §11.2** |
| 1092924_* | era5-dl | **8 FAILED**, 4 ok | 0:27–29:27 | 4G | **the dropped-`%4` incident** (§11.9) |
| 1092941_* | era5-dl | COMPLETED ×12 | 0:08–23:15 | 4G | re-run at `%4`: 8 downloaded, **4 skipped in 9 s** |
| 1092956_0–14 | cat-diag-trend | COMPLETED ×15 | 0:00–8:31 | 24G | NA diagnostics, DJF 1990/2000/2010/2020 (12 real, 3 skipped) |
| 1092975 | cat-trend | COMPLETED | 5:49 | 16G | the five-season trend check — §12.8 |
| 1092983_* | era5-dl | COMPLETED ×12 | 5:10–26:54 | 4G | **1985/1995/2005/2015 — 12/12 at `%4`** |
| 1093047_0–26 | cat-diag-trend | COMPLETED ×27 | 0:00–8:30 | 24G | **all 9 seasons; 12 real, 15 skipped in ≤1 s** |
| *(see `logs/trend-*.out`)* | cat-trend | COMPLETED | 6:08 | 16G | **the nine-season trend check — §12.1–12.6** |

**Aggregate: 39 months downloaded, 39 months of diagnostics, one global
calibration, two trend fits — roughly 9 hours of compute and 41 CDS requests.**
Against phase 4's 504 requests and ~15–23 h (§13.2).

### 13.4 Job scripts

```
jobs/00_smoke_test.sbatch              7 pre-flight checks
jobs/01_download.sbatch                NA months, array 0-503%4
jobs/02_download_global.sbatch         global 2000, full months (unused)
jobs/02b_download_global_calib.sbatch  global 2000, days 1/9/17/25, array 0-11%4
jobs/03_diagnostics.sbatch             NA diagnostics via 3_pipeline.py (superseded)
jobs/04_diagnostics_global.sbatch      global diagnostics, --chunk-days 1
jobs/05_diagnostics_na_djf.sbatch      NA DJF 1979
jobs/06_calibration_check.sbatch       the stage-2 gate
jobs/07_diagnostics_na_trend.sbatch    NA diagnostics, 9 DJF seasons, array 0-26
jobs/08_trend_check.sbatch             the trend comparison
```

Production uses `01` for raw and `07`'s driver (`ada/diagnostics_global.py`)
for derived — **not** `03`, which writes both formats.
