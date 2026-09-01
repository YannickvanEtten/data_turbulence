# RUNBOOK — the FORMULA_AUDIT checks

What changed in the code on 2026-08-29, what to run, in what order, and what
each result licenses. Companion to `FORMULA_AUDIT.md`; same shape as
`RUNBOOK_calibration_check.md`.

---

## 0. The one-paragraph version

Push, pull on ADA, run the tests, then submit the checks. Nothing here needs a
new CDS request and nothing here can invalidate a byte of a running download —
the raw ERA5 fetch is diagnostic-independent, and it is the only 36-hour item.

```bash
export BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
cd $BASE/data_turbulence          # NOT $BASE — git fails there, STATUS §10.1
module load 2025 && module load pixi
git pull
mkdir -p logs                     # before any sbatch, STATUS §11.9

sbatch jobs/12_tests.sbatch       # the gate — wait for this to be green
sbatch jobs/09_audit_checks.sbatch
sbatch jobs/10_f2d_variants.sbatch
sbatch jobs/11_ubf_diagnosis.sbatch
```

**Concurrency.** The default QOS caps you at **8 running jobs** (STATUS §11.9).
A download array at `%4` holds four, leaving four — exactly enough for `12`,
`09`, `10` and `11`. None of them needs to wait for the download.

**Two ADA facts these scripts encode, learned the hard way.** Modules do not
persist onto compute nodes, so a bare `srun ... pixi run ...` dies with
`execve(): pixi: No such file or directory`; every job script re-loads its own
(STATUS §10.2). And `mkdir -p logs` belongs on the command line, not inside a
script — SLURM opens the `--output` file before the body runs.

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

### 1.2 `f2d` has a variant flag

> **Superseded 2026-08-30: the default is now C, not A.** The evidence is in
> `FORMULA_AUDIT.md` §10.4 — decisively, Williams (2017) Figure 1 plots this
> diagnostic anchored at zero on 0..300 in the same figure where Negative
> Richardson runs −300..0. It is a magnitude. The paragraph below describes the
> state on 2026-08-29, when the flag was added and A was still the default;
> §5 is the current instruction.

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

### Step 0 — `jobs/12_tests.sbatch`  (under a minute, 8 GB)

The cheapest possible place for the two code changes to fail. **Do not queue
09–11 until this is green.** Expect **69 tests** as of 2026-08-30: the original
56, 12 in `tests/test_audit_fixes.py`, and one added to `tests/test_analytic.py`
when the f2d default moved to C. None of them need ERA5.

```bash
mkdir -p logs && sbatch jobs/12_tests.sbatch
tail -f logs/tests-*.out
```

If `TestDeformationUnsquared` fails, the `2_diagnostics`/`3_pipeline` pairing is
wrong. If `TestF2dVariants` fails, the A9 variant flag is. Anything else means
the change broke something it should not have touched.

### Step 1 — `jobs/09_audit_checks.sbatch`  (minutes, 16 GB)

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

### Step 2 — `jobs/10_f2d_variants.sbatch`  (~30 min, 24 GB)

Measures which of the four readings is distributed like the published
diagnostic, on the global calibration month, in both boxes.

**What it cannot do:** separate A from B. `B = -A`, so every symmetric statistic
agrees; what differs is which tail a threshold selects, and the tail-overlap
matrix in the output shows those sets are essentially disjoint. Separating them
needs the sign of the *trend* — step 4.

### Step 3 — `jobs/11_ubf_diagnosis.sbatch`  (~15 min, 16 GB)

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

## 3. Running during the download

All four jobs are safe to run while `jobs/01_download.sbatch` is fetching the
production series.

- **Nothing they read is being written.** `09` reads `derived/`, which the
  download does not touch. `11` reads `raw/north_atlantic/*.grib`, which the
  download *is* writing into — but it writes `<name>.tmp` and atomically
  renames (STATUS §10.7), so a complete file is never seen half-written, and
  months already on disk are skipped rather than rewritten.
- **Nothing they change is being read.** The download runs `1_download_hpc.py`
  and `download_plan.py`; neither was modified. After any `git pull`, confirm
  with `pixi run python -m py_compile 1_download_hpc.py download_plan.py` —
  with hundreds of array tasks still queued, each reading those files as it
  starts, a working tree left in a conflicted state would take them all down.
- **The one thing not to do** is pull midway through a *diagnostics* array.
  Tasks that have not started yet would pick up new code, and the dataset would
  be built half one way and half the other.

---

## 4. What was and was not exercised before this push

Exercised on synthetic data in a scratch container, end to end:

- `ada/check_composition.py` — both DEF conventions, **plus a negative control**
  with a mis-wired `ngm1` (2 % perturbation) and a doubled `brown2` constant.
  Both faults were caught, the three unaffected identities stayed green, and
  the job exited 2. A check that cannot fail is not a check — STATUS §12.8
  records that lesson being paid for once already.
- `ada/check_shape_ratios.py` — both sampling modes, box subsetting on
  descending ERA5 latitudes, cos φ weighting verified cell-by-cell under a
  transposed dimension order, and the graceful fallback when `2_diagnostics`
  cannot be imported.
- `check_f2d_variants.py`'s **chunk-boundary arithmetic**, exhaustively: for
  every timestep count and chunk size, the retained absolute indices are
  exactly `{1 … n-2}`, with no gaps and no duplicates. The chunked result is
  therefore identical to an unchunked one, which is what lets `--chunk-days 1`
  hold peak memory at 24 GB instead of asking for a request that would sit
  hours in the queue.

**Not executed anywhere yet**, because all need rojak and real GRIB:

- `ada/check_f2d_variants.py` end to end (only its chunking maths was proved)
- `ada/check_ubf.py`
- the code changes in `2_diagnostics.py` and `3_pipeline.py`

They are syntax-checked and follow the patterns of the scripts that were
exercised, but the first run of steps 2 and 3 on ADA is also their first real
run. `check_ubf.py` guards itself: it cross-checks its own term decomposition
against `2_diagnostics.ubf()` by rank correlation and refuses to print anything
if they disagree.

Which is why **step 0 exists** and why the four job scripts discover their
inputs with a glob rather than hard-coding a filename — STATUS §13.1 records a
derived directory being deleted on 2026-08-29, and a hard-coded path is exactly
how a script starts failing for a reason that has nothing to do with what it
tests.

---

## 5. The production diagnostics run — `jobs/14`

Added 2026-08-30, once `f2d` was resolved (`FORMULA_AUDIT.md` §10.4) and the
default switched from variant A to C.

### 5.1 Run it alongside the download, not after it

The raw fetch is CDS-bound and takes days; the diagnostics are CPU-bound and
take hours. There is no reason to serialise them, and `jobs/14` is built to
interleave: a task whose GRIB has not landed yet prints `NOT-YET-DOWNLOADED`
and exits 0, so you submit the same array again when the download finishes and
it picks up only what is left.

```bash
sbatch jobs/12_tests.sbatch                              # 1. the gate
sbatch --array=0-11%4 jobs/14_diagnostics_production.sbatch   # 2. canary
sbatch jobs/14_diagnostics_production.sbatch             # 3. the full 504
#   ... and again after jobs/01 completes
```

**Do step 2.** Twelve months costs half an hour and catches a configuration
mistake before five hundred months of it. Check one output before releasing
the rest:

```bash
pixi run python ada/check_composition.py \
    $BASE/derived/north_atlantic/diagnostics_na_1979-01.zarr
```

Expect all five identities at ~1e-8, and — this is the new part — the
deformation line should now read **DEF, un-squared**, not DEF². If it still
says DEF², the run picked up stale output and `--skip-if-matching` is not doing
its job.

### 5.2 Why "output already exists" is not a safe skip any more

`jobs/07` skips whenever the output directory is present. That was correct when
every output had been built identically. It is wrong now, because two
conventions changed underneath:

| written before | holds |
|---|---|
| 2026-08-29 | `deformation` = DEF², not DEF (§5) |
| 2026-08-30 | `f2d` = variant A, not C (§10.4) |

The 27 DJF months already in `derived/north_atlantic/` are **both**. A
directory check would preserve them and hand back a 504-month series that is 27
months of one definition and 477 of another — a real, plausible-looking dataset
quietly built from two conventions, which is the same failure shape as
Q-INTEG-3's silent level substitution.

So `ada/diagnostics_global.py --skip-if-matching` reads the output's recorded
provenance instead of its existence, and distinguishes three states:

```
no output                     -> compute
output from an older config   -> RECOMPUTE
output from this config       -> skip in ~1 s
incomplete store              -> RECOMPUTE
```

Unit-tested against all five cases before the first submission.

### 5.3 Concurrency

The QOS caps you at 8 running jobs and the download array holds 4. `jobs/14`
asks for `%4`, which fills the cap exactly and leaves nothing spare. If you
need a slot for a check job:

```bash
scontrol update JobId=<jobid> ArrayTaskThrottle=3
```

Expect some stretching from filesystem contention — §13.2 measured six
concurrent tasks turning a 30-minute job into 36, mostly in GRIB decode.

### 5.4 What still has to be re-run afterwards

`jobs/14` covers the North Atlantic series. Two things are still on variant A
and need re-running before the replication numbers can be quoted under the new
configuration:

```bash
sbatch --array=0-11%4 jobs/04_diagnostics_global.sbatch   # global calibration
sbatch jobs/06_calibration_check.sbatch                   # thresholds
sbatch jobs/08_trend_check.sbatch --per-diagnostic        # the trend
```

`jobs/04` overwrites by default, so no flag is needed — the default variant is
now C. `jobs/08` reads whatever `jobs/14` wrote.

**The prediction to check when `08` reports** (registered in `FORMULA_AUDIT.md`
§10.8, before any of this was run): under variant A, `f2d`'s DJF 1979 MOG
exceedance was **0.026 %**, second-lowest of the 21 and ~25× below the median
of 0.642 %. Under C it should land in the **0.1–1 %** band where its siblings
sit. If it does not, C is wrong too and the clipped variant
`max(±½ D/Dt[Q], 0)` is next.

Also worth watching: whether the §11 calibration agreement (24 %) and the §12
trend result move at all. One diagnostic out of 21 should shift the ensemble by
roughly 5 %, no more. A larger move would mean `f2d` had been doing more work
in the ensemble than its share, which would itself be worth knowing.
