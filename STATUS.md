# STATUS — CAT replication

**Single source of truth for this repo.** The code cites audit IDs (`Q-CALIB-2`,
`Q-GLOBAL-1`, `Q-AGG-3`, `STATUS_24`, …) from working sessions whose notes are
no longer findable. Those references are not retrievable; this file replaces
them. When a decision gets made, it gets written here, not into a docstring.

Last updated: 2026-08-27 (ADA environment verified).

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
aggregate.py                   Layers 3-4  exceedance-first, then average
annual_aggregate.py            Layer 5  leap-aware annual normalisation
trend.py                       Layer 6  per-gridpoint OLS 1979–2020
chunk_stitch.py                month-boundary overlap buffer for d/dt

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

**Not done.** No month of real ERA5 has been through the full stack, and the
downloader has never run on ADA. **The dataset does not exist yet.**

The ADA environment itself, however, is **verified and documented** — see §10.
The cfgrib concern its docstring raised is now **resolved** (§10.8).

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
  calibration (a near-constant factor moves the threshold by the same factor),
  but it should be quoted with the caveat.
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
no code change can repair it. It matters because turbulence-relevant features
evolve fast; F2D should be read as systematically damped for anything
sub-diurnal. Worth an explicit sentence in any write-up.

### 4e. The gradient foundation

`test_gradient_metric_matches_exact_ellipsoid` checks the operator every other
gradient-based diagnostic sits on. rojak reaches ∂/∂x via a nominal equatorial
grid spacing times PROJ's `parallel_scale`; compared against the exact WGS84
metric computed from the defining constants, it agrees to **6 significant
figures** zonally and **0.2 %** meridionally, with cross-terms at 1e-20. The
meridional 0.2 % is a slight redundancy between rojak's geodesic `dy` and its
`meridional_scale`; harmless at this magnitude. That 0.2 % is the noise floor
for the whole suite.

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
**latitude-structured** — which is exactly the class of error that does *not*
wash out of a percentile calibration (§5.4) and *does* contaminate spatial
patterns and their trends.

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
match exactly. The residual disagreement had been precisely this geometry.

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
pip install pytest                      # rojak requires Python >= 3.12
python -m pytest tests/ -v              # 25 tests, ~4 s
python tests/report_errors.py           # writes verification/<date>/analytic_verification.csv
python 4_verify.py                      # the cross-implementation check
```

### 4a. Measured results — run of 2026-08-26

Input `era5_validation_subset.nc`. Full output in
`verification/2026-08-26/` (two CSVs plus the console log).

Check (b), hand-written vs corrected rojak:

| Diagnostic | median ratio | Spearman ρ | verdict |
|---|---|---|---|
| `negative_richardson` | 1.000 | **1.0000** | MATCH |
| `colson_panofsky` | 1.000 | **1.0000** | MATCH |
| `ncsu1` | 1.000 | **1.0000** | MATCH |
| `brown2` | 1.000 | **1.0000** | MATCH |
| `ubf` | 0.993 | **0.9969** | CLOSE |

**The README's "ρ ≥ 0.99, two of them exact" claim holds, and understates the
result** — four of the five are exact to four decimal places, not two. This is
now evidenced on disk rather than asserted.

**Two gaps this run exposed, both properties of the validation FILE, not the code:**

1. **`f2d` was computed on nothing.** It came out with shape `(61, 61, 0)` —
   zero timesteps. `frontogenesis_isentropic` drops the first and last timestep
   for its centred material derivative, and the validation subset only has
   **two** timesteps, so nothing survives. F2D is currently unverified by any
   check at all, and the old dimensional "PLAUSIBLE" verdict for it was
   computed over an empty array.
2. **The negative-Ri check (c) never fired.** Zero negative-Ri cells at any of
   the three levels; max of −Ri is −0.324, i.e. Ri > 0 everywhere in this file.
   So the negative-Ri branches of `negative_richardson`, `colson_panofsky` and
   `ncsu1` are **not** covered by the ρ = 1.0000 result above — only their
   positive-Ri branches are.

**Action (phase 1):** replace `era5_validation_subset.nc` with a subset that has
**at least 3 timesteps** and comes from a **winter day** likely to contain
genuine Ri < 0 cells. Both gaps close with a better input file; neither needs a
code change.

---

## 5. What can actually invalidate the results

Ranked by impact on the trend estimates, which are the only numbers that matter.

1. **Rank-order errors** — a wrong term changes *which* cells are flagged, not
   just by how much. The F2D `dv_dy`→`dv_dx` bug was this kind. Tier B has no
   defence against it.
2. **Sign / tail errors** — the UBF sign flip was this. Cheap to re-derive all
   21 `REFERENCE_TABLE` sign entries from the source papers; unrecoverable
   downstream if wrong.
3. **Time-varying bias** — anything whose error drifts across 1979–2020 goes
   straight into the slope. Includes ERA5's own observing-system changes as
   satellite instruments enter the assimilation stream. A standard referee
   objection; needs an explicit answer before publication.
4. **Constant magnitude offsets — LOW PRIORITY.** Thresholds are percentiles of
   the data itself, so a uniform scaling `c·D` scales the threshold by the same
   `c` and leaves the exceedance field **identical**:

   ```
   D′ ≥ percentile_p(D′)  ⟺  c·D ≥ c·percentile_p(D)  ⟺  D ≥ percentile_p(D)
   ```

   This holds for any strictly increasing transformation. So the
   `1–2 ORDERS OFF` and `FLAG: magnitude mismatch` entries in
   `cat_outputs/comparison_table.csv` are a smell test against Williams &
   Joshi's published medians, **not** a correctness criterion for Prosser's
   method. The same argument largely defuses the 175/200/225 vs 188/197/206 hPa
   substitution — *provided* the ~2.8× damping is roughly constant rather than
   varying with latitude, level or time. That proviso is worth checking; it is
   the one way a magnitude problem becomes a real problem.

---

## 6. Plan

### Phase 0 — legibility  ✅ COMPLETE (2026-08-26)
- [x] **Fixed a crash that stopped the whole pipeline running.** Both
      `3_pipeline.py` and `4_verify.py` load `2_diagnostics.py` by file path
      with `importlib`, but never registered it in `sys.modules`. Because
      `2_diagnostics.py` uses `from __future__ import annotations`, `@dataclass`
      resolves its field annotations through `sys.modules[cls.__module__]` —
      which was `None` — and the import died on Python 3.12 with a bare
      `'NoneType' object has no attribute '__dict__'`. One line per loader
      (three in total). Confirmed: both files now import cleanly with all 21
      diagnostics and all six layers wired.
- [x] Fix the hardcoded `/home/claude/work/phaseB_out` output path in
      `4_verify.py`; results now land in `verification/<date>/`, dated so a
      re-run never overwrites earlier evidence, and with the console log
      captured (checks b and c printed their rho values to stdout and nowhere
      else).
- [x] Add `.gitignore` so `git add .` stops being dangerous.
- [x] This file.
- [x] Made `summarize()` in `4_verify.py` empty-safe — it crashed in
      `np.quantile` when a diagnostic had no finite values, which is exactly
      what happens to `f2d` on a 2-timestep file. It now reports `n=0` instead
      of killing the run halfway through.
- [x] **Ran `python 4_verify.py`; evidence committed to
      `verification/2026-08-26/`.** Results and the two gaps it exposed are in
      §4a above.

### Phase 1 — close tiers B and C  ✅ MOSTLY COMPLETE (2026-08-26)
- [x] Built the manufactured-solution suite in `tests/` — 22 tests covering all
      16 diagnostics that had no numerical evidence. See §4.
- [x] Convergence tests, which are the stronger evidence: they show the error
      falls as the grid refines, proving a correct second-order scheme rather
      than a wrong formula that happens to land nearby on one grid.
- [x] Verified the gradient metric itself against the exact WGS84 ellipsoid.
- [ ] **Build a better validation subset** for `4_verify.py`: ≥3 timesteps and a
      winter day, so `f2d` and the negative-Ri branches get exercised there too
      (see §4a). The analytic suite covers `f2d` now, so this is no longer
      urgent — but the negative-Ri branches of `colson_panofsky`,
      `negative_richardson` and `ncsu1` remain untested on real data.
- [ ] **Confirm the formulas against the papers** (§7). The suite proves the code
      computes what the docstrings say; it cannot prove the docstrings match
      Prosser. Needs `Articles/Sharman et al (2006).pdf`,
      `Williams and Joshi (2013).pdf`, `Ellrod Knapp (1992).pdf`.

### Phase 2 — pilot: one month, end to end, on ADA
- [x] Confirm cfgrib opens all 7 variables as one hypercube — **done, §10.8**.
      The integrity check in `1_download_hpc.py` also passes on a real GRIB.
- [x] **Compute-node egress — CONFIRMED (§10.7).** CDS reachable from
      `node222`, no proxy. Batch downloading is viable; job array it is.
- [x] **`chunk_stitch.py` filename bug — FIXED.** Its `filename_pattern`
      defaulted to `.nc` while the downloader writes `.grib`, so it could never
      have found a file the pipeline produces. Now defaults to `.grib` and
      dispatches on extension via a new `_open()` that uses the cfgrib engine
      explicitly and suppresses the `.idx` sidecar.
- [x] **Job scripts written**: `jobs/00_smoke_test.sbatch` (7 checks, one real
      CDS request) and `jobs/01_download.sbatch` (array `0-503%4`).
- [ ] Accept the ERA5 pressure-levels Terms of Use on the CDS website.
- [ ] Build the real Pixi env (rojak pinned to `1a65326`) in a dedicated folder
      on the project share, not the throwaway `~/env-test`.
- [ ] Run the smoke test, then launch the download array.
- [ ] Run Layers 2–6 with `validate_complete=False` and
      `calibration_domain='regional-plumbing-test'`.

### Phase 3 — decisions (needs supervisor)
- **What the production run persists.** *Largely resolved by §10.6:* raw
  (188 GB) plus full-resolution derived (1.13 TB) is 1.32 TB against a 2.5 TB
  share, so everything can be kept at 3-hourly resolution and the tail work
  gets its magnitudes without a second 42-year run. What remains is a
  conversation with the Koopman group about taking >50 % of a shared volume.
- Calibration domain: Prosser's global year-2000 pull, or regional with a
  stated caveat.
- Pressure levels: accept 175/200/225 as documented approximation, or pursue
  model levels.
- Box: North Atlantic or global (~28× cost).

### Phase 4 — the 42-year run
CDS-queue bound (504 monthly requests), not compute bound.

### Phase 5 — econometrics

**Design constraint, decided in phase 3, not phase 5.** Layers 3–6 are the
*replication check*, not the research dataset. `aggregate.py` produces **binary**
0/1 exceedance fields by design — magnitude is discarded. That is correct for
Prosser and **fatal for extreme-value work**: a GPD cannot be fitted to an
indicator variable, and peaks-over-threshold needs the size of each excess.

Rough storage for the North Atlantic box (121 × 301 grid, 3 levels, 21
diagnostics, ~122,700 timesteps, float32):

| What is saved | Size | Supports |
|---|---|---|
| Everything raw, 3-hourly | ~1 TB (less in practice; several diagnostics collapse to one level) | anything |
| Monthly per-gridpoint quantiles (p50/p90/p99/p99.9/max) | ~28 GB | trends, coarse tail work |
| Peaks-over-threshold archive (top 1%, magnitudes kept) | ~11 GB | EVT / GPD fitting |
| Annual exceedance grids only (current behaviour) | trivial | replication only |

Doing both of the middle two costs ~40 GB and covers nearly everything either
research direction needs. Storing annual grids alone means re-running 42 years
to get anything finer.

---

## 7. Open: do the formulas match the papers?

The analytic suite verifies code against docstring. It cannot verify docstring
against literature. These are the specific claims that still rest on a reading
rather than on a check, in rough order of how much a mistake would cost:

- **`ti2` sign convention.** `Sv × (DEF − δ)`. Ellrod & Knapp write it with
  convergence, i.e. `DEF + CVG` where `CVG = −δ`. rojak's implementation matches
  that. Worth one look at Ellrod & Knapp (1992) to be certain the sign is right,
  since a flipped sign here would be invisible in magnitude and fatal in rank.
- **`brown1` coefficient.** `sqrt(0.3 ζ_a² + D_sh² + D_st²)`. The 0.3 is
  Brown (1973) via rojak's docstring. Not independently checked.
- **`nva` and `ncsu1` clipping.** Both apply `max(·, 0)`. Confirm Sharman (2006)
  A36 and the NVA definition really are one-sided.
- **`f2d`** is the isentropic A9 form, chosen over the Miller physical-space form
  on a citation-chain argument (Prosser → Williams 2017 → Sharman 2006) plus a
  units argument. The units argument is strong; the citation chain deserves one
  confirmation pass against the actual Sharman appendix.
- **Every `sign` entry in `REFERENCE_TABLE`**, which decides which tail counts as
  turbulent. Two were already found wrong (`colson_panofsky`, `negative_richardson`
  were tagged `"either"` and corrected to `"+"`). A cheap, high-value re-derivation.

## 8. Open questions

- Intended econometric unit of observation? (route × month, gridpoint × month, …)
- ADA: storage quota on `/scistor/SBE-EDS-ClimateKoopman/yen230`, partitions,
  walltime, and whether compute nodes have outbound internet — usually not,
  which pins downloading to the interactive nodes.
- CDS: credentials current, licence accepted, realistic queue throughput?
- Supervisor expectations and deadline.
- The `STATUS_24` / `Q-*` notes, if they ever turn up. Not a blocker — every
  decision in them is re-derivable from the papers and the code.

---

## 9. Repo facts worth remembering

- Git tracks 12 files / 58 KB, and they are the **old** versions now sitting in
  `old_code/`. Last commit 2026-05-26; files edited through August 2026.
- `cat_outputs/` is 421 files / 114 MB. Now git-ignored except its small CSVs.
- `.git/config` has an `[lfs]` section but there is no `.gitattributes`, so LFS
  is half-configured and will not catch large files. Do not rely on it.
- Everything lives in OneDrive, so uncommitted work is backed up and versioned.
- `rojak` requires Python ≥ 3.12.

---

## 9a. Repository and data layout

### The three-stage science plan this serves
1. **Pilot** — one North Atlantic month, end to end. Proves the chain on real data.
2. **Calibrate** — global year 2000 → severity thresholds → **compare against
   the published tables**. This is the step that validates the whole chain.
3. **Analyse** — North Atlantic 1979–2020, thresholds applied, trends computed.

Stages 2 and 3 use *different datasets*. `3_pipeline.run_layers_2_to_6` already
took `calibration_fields` separately from `diagnostic_fields`; what was missing
was the plumbing to produce both, and a way to carry thresholds between them.

### On-disk layout
```
/scistor/SBE-EDS-ClimateKoopman/yen230/
├── data_turbulence/     the git repo — CODE ONLY, nothing large
├── env/                 the pixi environment (pixi.toml, pixi.lock)
├── raw/
│   ├── north_atlantic/  era5_na_YYYY-MM.grib      504 files, 188 GB
│   └── global/          era5_glob_2000-MM.grib     12 files, 128 GB
├── derived/
│   ├── north_atlantic/  diagnostics_na_YYYY-MM.zarr
│   └── global/          diagnostics_glob_2000-MM.zarr
├── calibration/         thresholds_YYYY-MM-DD.json   <- stage 2 -> 3 handoff
├── results/             annual probabilities, trends
└── logs/                SLURM output
```
Total raw: **316 GB** of the 2.5 TB share.

### Changes made 2026-08-27

**1. Domains are named, not ad-hoc boxes.** `download_plan.DOMAINS` defines
`north_atlantic` and `global`; `1_download_hpc.py` takes `--domain`. The two
downloads are now symmetric commands rather than "one is the default and the
other is an `--area` flag you must remember to get right". `--area` survives as
an explicit override for one-off experiments only.

**2. FIXED — a silent data-corruption bug.** `month_filename()` did not encode
the domain, so global and North Atlantic files for the same month shared a name
(`era5_2000-01.grib`) and **one would overwrite the other**. Worse, the
integrity check could not tell them apart: both have exactly 7 variables, 3
levels and 8×n_days timesteps. Only the grid differs. Demonstrated on the real
2016-07-02 GRIB — `check_counts` reports *nothing wrong* when an NA file is
checked against a global request.

Two fixes, both needed:
- filenames now carry a domain code: `era5_na_2000-01.grib`, `era5_glob_2000-01.grib`
- `check_grid()` added to the integrity check, verifying the horizontal grid
  against `expected_grid(area)` (NA 121×301, global 721×{1440,1441} — the pair
  allows for the wrap-around meridian being repeated or not).

**3. NEW `calibration.py` — thresholds are a file now.** Previously the
thresholds existed only as a local variable inside `run_layers_2_to_6`,
computed and consumed in one call. Stage 2 of the plan requires them to be
inspectable, comparable to the literature, and reusable weeks later without
recomputing a global year. The module provides `compute_thresholds`,
`save_thresholds` / `load_thresholds` (JSON, schema-versioned, with provenance:
domain, period, levels, sample sizes, rojak rev), and `compare_to_published`.

`compare_to_published` reports a ratio against Williams & Joshi Table 1 and
labels it, because there are **three** outcomes and only one is good:

| ratio vs published | verdict |
|---|---|
| 0.1 – 10, same sign | `consistent` — the calibration works |
| exactly 1.0000 | `SUSPICIOUSLY EXACT` — and `aggregate.assert_thresholds_not_hardcoded_table2` will refuse to run |
| orders apart, or opposite sign | `ORDERS OFF` / `SIGN MISMATCH` |

**4. Job scripts** `jobs/00_smoke_test.sbatch`, `01_download.sbatch` (NA,
array 0-503%4), `02_download_global.sbatch` (global 2000, array 0-11%2,
16 GB memory, and a standing instruction to run `--array=0` alone first).

### Deliberately NOT changed yet
The numbered module names (`2_diagnostics.py`) cannot be imported normally,
which is what caused the `sys.modules` crash. Renaming into a package would
end that class of bug permanently — **but not before the pilot runs.** The
verification records, job scripts and this document all reference the current
names, and renaming mid-flight before a single month has been processed is how
a project loses its thread. Revisit after stage 1.

### Open sizing question
A global month is **10.6 GB and 5,208 fields** in one CDS request. Copernicus
enforces a per-request limit whose current value is unconfirmed. If a whole
month is rejected, split into daily requests (168 fields) —
`download_plan.build_request` already accepts a `days` list, so that is a
change to the job script, not to the library.

---

## 10. The ADA environment (verified on-cluster 2026-07-06)

Recovered from a prior on-cluster investigation. Everything here was verified
live, not assumed. This section is the answer to "how do we run this".

### 10.1 Access
- ADA is the **renamed BAZIS**; old `bazis.readthedocs.io` links redirect to
  `rdm.vu.nl/tools/ada/`. Support: `itvo.ucit@vu.nl`.
- Account `yen230` (uid 722622). **Primary group is the Koopman group**
  (`R_SciStor-SBE-EDS-ClimateKoopman_C`), so new files are group-owned by the
  project automatically.
- Route in: **MobaXterm + EduVPN**, jump host `ssh.data.vu.nl` →
  `ada.labs.vu.nl`. **Two different passwords** — VUnetID at the gateway, a
  separate ITVO cluster password at ADA.
- Lands on `login01`/`login02`. **Login nodes auto-log-off every 5 minutes** if
  used for anything heavy — they are for editing and submitting only.
  Interactive work goes on `inter01`–`inter04`.

### 10.2 Modules
```bash
module load 2025 && module load pixi     # lowercase 'pixi'
```
The published docs say `module load Pixi` (capital P) — **that fails**.
Modules do not persist across sessions; batch scripts must re-load their own.

### 10.3 Scheduler — verified end-to-end (hello-world ran on `node220`)

| Partition | Time limit | Nodes |
|---|---|---|
| `defq` *(default, use this)* | **7 days** | 17 |
| `defq-thin` | 7 days | 4 |
| `defq-fat` | 7 days | 2 |
| `defq-gpu` | 7 days | 6 |
| `bw` | infinite | 3 (restricted) |

**The 7-day cap is a design constraint, not a detail.** 504 monthly CDS
requests cannot run as one job. `1_download_hpc.py --index` already exists for
exactly this — one array task per month.

### 10.4 Software — the reference environment
Pixi 0.46.0, **Python 3.12.13**. `python -m cfgrib selfcheck` → *"Found:
ecCodes v2.47.0. Your system is ready."* `import rojak` succeeds.

**rojak is pinned to rev `25b8685`** (`rojak_cat 1.0.2.dev20+g25b8685c6`,
2026-08-17 — HEAD at the time of writing). See `pixi.toml` in the repo root for
the buildable environment file.

> **Pin moved forward from `1a65326` on 2026-08-27, deliberately.** `25b8685`
> is *"FIX: Geospatial laplacian was missing an extra derivative"* (#234). At
> `1a65326`, rojak's `spatial_laplacian` returned `df/dx + df/dy` — not a
> Laplacian at all. Upstream now routes it through `divergence()`, which uses
> `vector_derivatives()` and therefore carries the spherical curvature term.
> **That is the same correction this project derived independently from the
> divergence theorem the day before (§4g).** Two independent routes, one answer.
>
> **Nothing regresses.** The suite was run at both revs: identical error
> magnitudes, 25/25 tests pass, all 16 analytic diagnostics unchanged, all 5
> cross-checks at ρ = 1.0000.
>
> **The cross-check got stronger.** Because upstream's Laplacian is now
> correct, `4_verify.py` no longer monkey-patches it. UBF still matches at
> ρ = 1.0000 — but now against rojak's *own* implementation rather than
> against a second copy of our formula. Every monkey-patch removed from the
> oracle makes that comparison more independent, not less.
>
> **Pin a SHA, never a branch.** The pin is the reproducibility record.

**Six of the seven Phase-A findings still stand at `25b8685`**, so all seven
hand-written diagnostics remain necessary. Re-confirmed against the source:

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
sparse 0.19.0. (Python is **3.12+**, not ≤3.9 — an earlier note claiming
otherwise was wrong; `1a65326` declares `requires-python >=3.12`.)

### 10.5 Storage — SciStor, no auto-purge, no `quota` command

| Use | Path | Size | Used |
|---|---|---|---|
| Code, scripts, Pixi env | `/scistor/guest/yen230` (home) | 200 GB | ~1.5 GB |
| **ERA5 data + outputs** | **`/scistor/SBE-EDS-ClimateKoopman/yen230/`** | **2.5 TB** | ~0 |

The Koopman share is already one-folder-per-member; slot into `yen230/`, no
request needed. There is a `.snapshot` directory (point-in-time copies) — but
**confirm retention and whether real backups exist** before treating derived
output as safe. Raw ERA5 is always re-downloadable, so it is lower stakes.

### 10.6 Sizing — measured, not estimated (2026-08-27)

Settles the open 250 GB vs 500 GB question. Measured from
`climate_data_02_07_2016.grib` (a real file of exactly the locked
configuration) times the true calendar:

| Quantity | Value |
|---|---|
| One full day, 7 vars × 3 levels × 8 steps, NA box | **12.26 MB** (2.00 bytes/gridpoint — 16-bit GRIB packing) |
| Average month | **373 MB** |
| **Raw GRIB, 1979–2020** (15,341 days) | **188 GB** |
| Timesteps in the full run | 122,728 *(matches `chunk_stitch.py`)* |
| Derived, all 21 diagnostics × 3 levels, 3-hourly, float32 | **1.13 TB** |
| Derived, 21 diagnostics at 200 hPa only | 0.38 TB |
| Monthly per-gridpoint quantiles (5 stats) | 23 GB |
| Peaks-over-threshold, top 1 % with magnitudes | 23 GB |

**This changes the phase-3 "what to persist" decision.** Raw + full-resolution
derived is **1.32 TB against a 2.5 TB share** — it fits, with headroom. The
painful tradeoff between replication output and econometric output does not
have to be made: everything can be kept at full 3-hourly resolution, and the
tail-behaviour work gets the magnitudes it needs without a second 42-year run.
**Caveat:** the 2.5 TB is *shared* with other Koopman members (`info/`,
`skn370/`, `wqt200/`). Taking >50 % of it is a conversation with the group, not
a unilateral choice. The 46 GB compact option (quantiles + POT) remains the
polite fallback.

### 10.7 Network — and the one thing still unverified
- **Login node egress works.** conda-forge, GitHub and Copernicus all
  reachable; `curl -sSI https://cds.climate.copernicus.eu/api/` → HTTP 308.
- **cdsapi auth works from ADA.** `~/.cdsapirc` present, `chmod 600`, new
  format (`url: https://cds.climate.copernicus.eu/api`, key = Personal Access
  Token); `check_authentication()` returns the account profile.
- **✅ COMPUTE-NODE EGRESS CONFIRMED (2026-08-27).** Tested with an `srun`
  on `node222`:

  ```
  CDS  -> HTTP 202
  PyPI -> HTTP 200
  proxy vars: http_proxy=unset  https_proxy=unset
  ```

  **Compute nodes reach the internet directly, no proxy.** The architecture
  question is settled: downloading happens in a `defq` **job array**, one month
  per task, and there is no need for a long-lived process on the interactive
  nodes. See `jobs/01_download.sbatch`.

- **⚠ `--ntasks=1` IS MANDATORY, discovered from that same test.** Every line of
  the egress output printed **twice** — `srun` allocated two tasks and ran the
  command twice on one node. In a download job that means two processes
  requesting the same month from CDS and racing on the same `<name>.tmp`. The
  atomic-rename guard in `_fetch_verified` protects against a *half-written*
  file, **not against two concurrent writers**. Both job scripts set
  `--ntasks=1` explicitly, and `jobs/00_smoke_test.sbatch` §1 asserts it took
  effect.

### 10.8 cfgrib multi-variable open — RESOLVED 2026-08-27
`1_download_hpc.py`'s docstring flagged as untested whether cfgrib can open all
7 variables as a single hypercube, or whether a `filter_by_keys` split would be
needed. Tested against the real `climate_data_02_07_2016.grib`:

```
xr.open_dataset(path, engine="cfgrib", backend_kwargs={"indexpath": ""})
  -> SUCCESS, single open, no filter_by_keys
  data_vars: ['u','v','t','z','d','vo','pv']
  dims     : {'time': 8, 'isobaricInhPa': 3, 'latitude': 121, 'longitude': 301}
  domain   : lat 30.00..60.00, lon -75.00..0.00, 0.25 deg, levels [225,200,175]
```
The domain matches `download_plan.NORTH_ATLANTIC_BOX` exactly. The integrity
check in `1_download_hpc.py` (`verify_file`) was run against this file and
**passes with zero problems**. That risk is closed.

### 10.10 Environment built and verified on ADA — 2026-08-27

Job `1091550` on `node009`, exit status 0. `pixi install` from the repo's
`pixi.toml`, 1.2 GB, rojak pinned at `25b8685`.

```
python 3.12.14                                    [ok]
rojak-cat 1.0.2.dev20+g25b8685c6 -- matches pin   [ok]
cfgrib + ecCodes v2.48.0                          [ok]
all 21 diagnostics import and register            [ok]
~/.cdsapirc, new-format url                       [ok]
CDS authentication from a compute node            [ok]
cfgrib selfcheck: "Your system is ready."
25 passed in 6.34s          <-- the full verification suite, ON ADA
```

**The last line is the point.** The accuracy figures in `verification/` are now
reproducible in the environment that will produce the dataset, at the same
rojak rev, rather than only in the environment where they were first measured.

**Three bugs found getting there**, all mine, all now fixed:

1. **`pixi run python -c "..."` is not safe for multi-line Python.** `pixi run`
   parses the string with its own task-shell parser before Python sees it, and
   that parser rejects English words that happen to be shell reserved words.
   The first run died on the word **"in"**, inside a Python string, inside an
   assertion message — with a perfectly healthy environment. All checks now
   live in `ada/verify_env.py` as a real file. The same pattern was waiting in
   `jobs/00_smoke_test.sbatch` in four more places.
2. **rojak has undeclared imports.** `dask-geopandas`, `requests` and `rich`
   are imported at module load but not declared, so a resolver working from its
   metadata builds a clean environment that then dies on
   `import rojak.core.data`. Now named explicitly in `pixi.toml`.
3. **`cdsapi.Client()` defaults to `retry_max=500, sleep_max=120`** — over
   **sixteen hours** of silent retrying on an unreachable CDS. In
   `1_download_hpc.py` that meant any month CDS could not serve would hold a
   compute node until SLURM killed it at walltime, printing nothing useful.
   Now bounded to 8 attempts backing off to a minute, then a loud failure —
   which is the right shape precisely *because* the downloader is resumable: a
   failed month leaves no final file, so re-submitting the array retries
   exactly the months that failed. Found by accident when a pre-flight probe
   hung against a blocked endpoint.

Also worth recording: `check_authentication()` does not exist across all cdsapi
builds (conda-forge returns a `LegacyClient` without it; PyPI 0.7.7 has no such
class). The probe now tries a couple of locations and degrades to a note rather
than failing setup — a missing or old-format `~/.cdsapirc` is still fatal.

### 10.11 First real month downloaded — 2026-08-27 23:58

Job `1091553_0`, `--array=0`, on `node009`. Smoke test passed first
(`--ntasks=1` confirmed, 25 tests passed on ADA, one trial-day download).

```
era5_na_1979-01.grib   363 MB   integrity check PASSED
request -> successful  6 min 55 s
```

**Measured timings for the full run.** 7 minutes per month at `--array=...%4`
is 126 batches, roughly **15 hours** of wall clock for all 504 — likely one to
two days once CDS is busier in European working hours than it was at midnight.
Size: 363 MB for a 31-day month is ~11.7 MB/day, giving ~183 GB for the full
raw archive, against the 188 GB predicted in 10.6. The storage plan holds.

**Not yet done:** run the diagnostics on this month and look at the output
before launching the full array. A single month cannot exercise
`chunk_stitch.py`'s month-boundary buffer, but it can show whether the
exceedance field is jet-stream-shaped rather than noise.

### 10.12 First month of diagnostics — 2026-08-28 00:15

Job `1091555` on `node013`, exit 0. `3_pipeline.py` on `era5_na_1979-01.grib`.
All 21 diagnostics computed, no failures, no NaN placeholders.

| | measured |
|---|---|
| wall time | **8 min 11 s** (vs 6 min 55 s to download the same month) |
| peak memory (`MaxRSS`) | **14.8 GB** — the job asked for 64 GB |
| output, NetCDF | 1.3 GB |
| output, Zarr | 1.1 GB |

**Download and diagnostics cost roughly the same per month, for opposite
reasons.** The download is idle waiting on Copernicus; the diagnostics are real
arithmetic on the node. They therefore want opposite job shapes: downloads
narrow (`%4`, CDS-limited, 4 GB), diagnostics wide (memory-limited, 24 GB).
504 months of diagnostics: 17 h at `%4`, **4.3 h at `%16`**.

`jobs/03_diagnostics.sbatch` now requests **24 GB** (1.6× the measured peak)
rather than 64 GB. On ADA memory is scheduled explicitly, so the request size
directly sets how many array tasks fit on a node — this is what makes `%16`
schedulable at all.

**Storage, settled.** One format × 504 months = **613 GB**; plus ~183 GB of raw
GRIB that is **796 GB, about 32 % of the 2.5 TB share.** Comfortable. Writing
both formats would be 1226 GB — `save_outputs()` writes the same data as
NetCDF *and* Zarr, and production must pick one.

**The output is 200 hPa only.** Visible in the log: rojak computes on
`(121, 301, 248, 3)` and `compute_all_21(target_level=200)` saves
`(121, 301, 248)`. That is correct — three levels give the vertical-derivative
stencil, the middle one is where diagnostics are evaluated — but it is an
implicit decision. **If the econometrics ever wants 175 or 225 hPa, that is a
re-run of all 504 months.** Worth deciding deliberately before phase 4.

**Comparison against Williams & Joshi Table 1: read it as season, not error.**
January in the North Atlantic at 200 hPa *is* the jet stream; W&J's medians are
a far broader climatology. Wind speed **29.8 vs 14.9 m/s** is the giveaway —
30 m/s is a textbook winter jet core. Everything downstream of the wind
inherits that (VWS 4.67 vs 1.88, DEF 43.8 vs 18.6, TI1 184 vs 31.5). Several
land close regardless: brown1 80.6 vs 77.1, PV 5.19 vs 6.84, NCSU1 8.98 vs 11.1.

**F2D moved from `FLAG: magnitude mismatch` to `PLAUSIBLE` (91 vs 56.6).** On
the 2-timestep validation file it was computed over an *empty array* (§4a);
with 248 timesteps its material derivative is real for the first time.

Still flagged `1-2 ORDERS OFF`: `colson_panofsky`, `rva_magnitude`, `ubf`,
`nva`. Worth examining — `ubf` especially, since its geometry changed on
2026-08-26 (§4g). But recall §5.4: a *constant* offset washes out of the
percentile calibration entirely. Only a structured error matters.

### 10.9 Still open on ADA
- **Compute-node egress** (§10.7) — the blocker.
- Snapshot retention / backup policy on the Koopman share.
- Node-local scratch on compute nodes (`$TMPDIR`) — for job I/O; minor.
- **Accept the ERA5 pressure-levels Terms of Use on the CDS website** —
  one-time, per-dataset, and the first real download fails without it.
- Build the real project env in a dedicated folder (the verified `~/env-test`
  was a throwaway).
- `1_download.py` has **two trial dates in circulation** — the uploaded copy
  uses 2016-01-01, the project-folder copy 2016-07-01. Pick one deliberately.
