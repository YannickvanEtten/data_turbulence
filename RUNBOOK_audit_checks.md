# RUNBOOK — the FORMULA_AUDIT checks

What changed in the code on 2026-08-29, what to run, in what order, and what
each result licenses. Companion to `FORMULA_AUDIT.md`; same shape as
`RUNBOOK_calibration_check.md`.

---

## 0. The one-paragraph version

Push, pull on ADA, and submit `jobs/09_audit_checks.sbatch` — it needs nothing
that is not already on disk and takes minutes. **Submit the production download
first and independently**: the raw ERA5 fetch is diagnostic-independent, it is
the only 36-hour item, and nothing in this runbook can invalidate a byte of it.

```bash
cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
git pull
mkdir -p logs

sbatch jobs/01_download.sbatch          # bare — NOT --array=..., see STATUS §11.9
sbatch jobs/09_audit_checks.sbatch      # while it runs
```

---

## 1. Code changes in this push

Two, and only two. Everything else is new files that read data and print.

### 1.1 `deformation` is no longer persisted squared

**Files:** `2_diagnostics.py` (`compute_all_21`), `3_pipeline.py`
(`PRETRANSFORM`, `SCALE_TO_TABLE` comment).

rojak's DEF diagnostic returns DEF²; Sharman A17 and every published table
define DEF = (D_SH² + D_ST²)^½. The square root used to be applied only in
`3_pipeline.py`, on the path that builds the comparison table — so every zarr
written by `ada/diagnostics_global.py` held DEF² under the name `deformation`.

- **Nothing in STATUS §11 or §12 changes.** Squaring is strictly increasing on
  a non-negative field, so the exceedance field is identical (STATUS §5.5).
- **It matters for phase 5.** A GPD fitted to DEF² is not a GPD fitted to DEF.
- `TI1`, `TI2`, `NGM1`, `NGM2` never had the problem — rojak hands those
  `CATData.total_deformation()`, already un-squared.

**The two halves must move together.** Un-squaring in `compute_all_21` while
leaving the √ in `3_pipeline` would take the fourth root and the comparison
table would drift silently. `tests/test_audit_fixes.py::TestDeformationUnsquared`
pins both.

**Existing zarr on disk is the old convention.** `ada/check_composition.py`
detects which convention a file uses and says so, so nothing depends on
remembering.

### 1.2 `f2d` has a variant flag; the default is unchanged

**Files:** `2_diagnostics.py` (`F2D_VARIANTS`, `frontogenesis_isentropic`,
`compute_all_21`), `ada/diagnostics_global.py` (`--f2d-variant`).

| variant | expression | |
|---|---|---|
| **A** | `+0.5 D/Dt[Q]` | current behaviour, **still the default** |
| B | `-0.5 D/Dt[Q]` | A9's leading minus applied |
| C | `\|0.5 D/Dt[Q]\|` | magnitude |
| D | `-D/Dt[sqrt(Q)]` | literal A9, normalisation included |

with `Q = (∂u/∂θ)² + (∂v/∂θ)²`. **The default reproduces every existing number
exactly**, so this push cannot move a published result by itself. The variant is
recorded in the output zarr's attributes, so a file can always be traced to the
reading that produced it.

Why the question is open at all: `FORMULA_AUDIT.md` §4.

---

## 2. What to run, in order

### Step 1 — `jobs/09_audit_checks.sbatch`  (minutes, 32 GB)

Runs the composition identities on one NA month, then the shape-ratio test
twice (daily-mean and instantaneous).

**What a pass licenses.** The composition identities verify nine of the 21
against each other cell by cell with no external reference — they catch a
mis-wired argument, a level mismatch or a chunk-boundary error, the class of
fault that produces a plausible field and survives every distributional check.
`brown2` in particular becomes fully closed: A14 confirmed against the paper
(`FORMULA_AUDIT.md` §2), magnitude explained (§6), composition confirmed
pointwise here.

**What to look at.** Median relative errors should sit at the float32 floor,
1e-7 to 1e-5. `1e-3` or above is a real disagreement. Read the p99 column too —
an identity that holds in the bulk and fails in the tail is exactly what a
percentile-calibrated pipeline would carry into the econometrics.

**Shape ratios: three predictions, written down in advance so the test can
fail.**

1. Most diagnostics land within ~2× of the published ratio.
2. Shear-based ones come in systematically **low** — a second, independent
   measurement of the stencil damping, from distribution shape rather than
   exceedance frequency.
3. `f2d` under variant A comes in **far too high** (hundreds, not 13.6). If it
   does, variant A is refuted without needing to read A9 correctly at all.

Exit code 2 from the composition step means an identity failed; that is a stop
sign, not a warning.

### Step 2 — `jobs/10_f2d_variants.sbatch`  (~30 min, 48 GB)

Measures which of the four readings is distributed like the published
diagnostic, on the global calibration month, in both boxes.

**What it cannot do:** separate A from B. `B = -A`, so every symmetric statistic
agrees; what differs is which tail a threshold selects, and the tail-overlap
matrix in the output shows those sets are essentially disjoint. Separating them
needs the sign of the *trend* — step 4.

### Step 3 — `jobs/11_ubf_diagnosis.sbatch`  (~15 min, 32 GB)

UBF at two months forty-one years apart. Reports term magnitudes and the
cancellation ratio, float32 vs float64 by rank, and ERA5's `vo` against the
vorticity derived from u and v.

**The number that matters** is the last one: `median |f·(ζ_ERA5 − ζ_derived)|`
against `median |UBF residual|`. Approaching 1 means UBF is substantially
measuring a product-versus-derivative inconsistency rather than atmospheric
imbalance — which would explain both the flat trend in STATUS §12.6 and the
disagreement with Prosser's Figure S5.

Compare the two months. A discrepancy that is the same size at both ends of the
record cannot by itself explain a flat *trend*; one that grows across the record
can, and is also a candidate for STATUS §5's risk 3.

### Step 4 — only if steps 2 or 3 point somewhere

A full re-run under the candidate configuration, then the trend comparison:

```bash
sbatch --array=0-11%8  jobs/04_diagnostics_global.sbatch    # + --f2d-variant X
sbatch                 jobs/06_calibration_check.sbatch
sbatch --array=0-26%8  jobs/07_diagnostics_na_trend.sbatch  # + --f2d-variant X
sbatch                 jobs/08_trend_check.sbatch --per-diagnostic
```

Then compare the per-diagnostic table against Prosser Figure S5
(`Articles/Prosser2023_FigureS5_per_diagnostic.png`), converting the relative
changes to **absolute percentage points** first — S5 is an absolute-change map.
It is a sign-and-ranking test, not a magnitude test: S5 is annual and global,
§12.6 is DJF in Prosser's box, and S5's colour scale saturates at ±0.5 pp while
our DJF `vertical_wind_shear` change alone is +1.14 pp.

`jobs/04` and `jobs/07` do not yet forward `--f2d-variant`; add it to the
`ada/diagnostics_global.py` line in each when step 4 is actually reached, rather
than now, so the scripts that produced the current results stay untouched.

### Step 5 — the gate before the diagnostics array

Only these change what reaches disk, so only these must be settled before
`jobs/07`'s production successor is submitted:

- [ ] the `f2d` variant (steps 2 and 4)
- [ ] whether UBF should compute `f·ζ` from the derived vorticity (step 3)
- [ ] the phase-5 persistence question — monthly quantiles and the
      peaks-over-threshold archive, STATUS §6

The `deformation` fix is already in and needs no decision.

---

## 3. Tests

```bash
python -m pytest tests/ -v
```

`tests/test_audit_fixes.py` adds 11 tests and needs no ERA5. Existing tests are
unchanged except one docstring in `tests/test_analytic.py`, which used to say
`3_pipeline.py` takes the square root — no longer true, and the pairing it
guarded now lives in the new file.

---

## 4. What was and was not exercised before this push

Run on synthetic data in a scratch container, end to end, including a negative
control with a mis-wired `ngm1` and a wrong `brown2` constant — both were caught
and the unaffected identities stayed green:

- `ada/check_composition.py` — both DEF conventions, plus the negative control
- `ada/check_shape_ratios.py` — both sampling modes, box subsetting on
  descending ERA5 latitudes, cos φ weighting under a transposed dimension order

**Not executed anywhere yet**, because both need rojak and real GRIB:

- `ada/check_f2d_variants.py`
- `ada/check_ubf.py`
- the code changes in `2_diagnostics.py` and `3_pipeline.py`

All are syntax-checked and follow the patterns of the scripts that were
exercised, but the first run of steps 2 and 3 on ADA is also their first real
run. `check_ubf.py` guards itself: it cross-checks its own term decomposition
against `2_diagnostics.ubf()` by rank correlation and refuses to print anything
if they disagree.

Run `python -m pytest tests/ -v` before the sbatch queue, since that is the
cheapest place for the two code changes to fail.
