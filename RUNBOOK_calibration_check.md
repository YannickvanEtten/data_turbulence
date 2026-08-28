# RUNBOOK — the calibration check on ADA

The 14-request experiment that tests the whole chain against Prosser (2023)
Table 1 before committing to the 504-month download.

Target: **North Atlantic DJF 1979, in Prosser's box (36–60°N, 55–10°W),
diagnostic-mean MOG exceedance ≈ 1.03 %.** Full acceptance table in §6.

---

## 0. Before you start — get the code onto ADA

The new pieces (`--days` in `1_download_hpc.py`,
`jobs/02b_download_global_calib.sbatch`, `tests/test_calibration_roundtrip.py`)
are on the laptop, not on the cluster.

```bash
# laptop
cd ~/OneDrive/Documenten/GitHub/data_turbulence
git add -A
git commit -m "calibration check: --days sub-sampling, global calib job, roundtrip test"
git push
```

```bash
# ADA, after logging in (see §1)
cd /scistor/SBE-EDS-ClimateKoopman/yen230/data_turbulence
git pull
git log --oneline -1        # confirm the commit arrived
```

---

## 1. Logging in

Per STATUS §10.1 — two different passwords, do not mix them up.

1. **EduVPN** connected.
2. **MobaXterm** → jump host `ssh.data.vu.nl` (VUnetID password)
   → `ada.labs.vu.nl` (separate ITVO cluster password).
3. You land on `login01` or `login02`.

**Login nodes auto-log-off after ~5 minutes of anything heavy.** They are for
editing and submitting only. Everything below is `sbatch`, which is fine —
nothing here runs interactively.

---

## 2. Set up the environment and directories

```bash
module load 2025
module load pixi                     # lowercase; 'module load Pixi' FAILS

BASE=/scistor/SBE-EDS-ClimateKoopman/yen230
cd $BASE/data_turbulence             # pixi.toml and .pixi/ live here

mkdir -p logs                        # <-- MUST exist BEFORE sbatch
mkdir -p $BASE/raw/global $BASE/raw/north_atlantic
mkdir -p $BASE/derived/global $BASE/derived/north_atlantic
mkdir -p $BASE/calibration $BASE/results
```

> **`mkdir -p logs` is not optional.** `#SBATCH --output=logs/...` is resolved
> relative to the directory you submit from, and SLURM opens that file *before*
> the job script runs — so the `mkdir` inside the script is too late. A missing
> `logs/` makes the job fail instantly with no output to tell you why.
>
> Modules do not persist across sessions; every job script re-loads its own.
> The two `module load` lines above are only for *your* shell.

Sanity check the environment is intact after the `git pull`:

```bash
pixi run python ada/verify_env.py
pixi run python -m pytest tests/ -q       # 25 existing + 31 new = 56
```

The 31 new tests in `tests/test_calibration_roundtrip.py` need no ERA5 and
take ~15 s. They check that a p-th percentile threshold, applied back to the
data it was calibrated on, is exceeded by exactly (100−p)% of it.

---

## 3. Submit the global calibration download

**One month first.** `--days` is new code — the request building and filenames
are unit-tested, but it has never touched CDS.

```bash
cd $BASE/data_turbulence
sbatch --array=0 jobs/02b_download_global_calib.sbatch
squeue -u $USER
```

Watch it:

```bash
tail -f logs/global-calib-*_0.out
```

What a good log looks like:

```
=== task 0 on nodeNNN at ... ===
SLURM_NTASKS=1   (must be 1)
global 2000-01 days=[1, 9, 17, 25]: download
=== done at ... ===
era5_glob_2000-01_d01-09-17-25.grib   ~1.4 GB
```

Then confirm the file is what it claims to be:

```bash
ls -lh $BASE/raw/global/
pixi run python -c "
import xarray as xr
ds = xr.open_dataset('$BASE/raw/global/era5_glob_2000-01_d01-09-17-25.grib',
                     engine='cfgrib', backend_kwargs={'indexpath': ''})
print(dict(ds.sizes))
print(sorted(ds.data_vars))
"
```

Expect `time: 32` (4 days × 8 steps), `isobaricInhPa: 3`,
`latitude: 721`, `longitude: 1440 or 1441`, and all seven of
`d, pv, t, u, v, vo, z`. The downloader's own integrity check already verified
this — if it had not, there would be no final file, only a `.tmp`.

**Then release the other eleven:**

```bash
sbatch jobs/02b_download_global_calib.sbatch
```

The array is `0-11%4`. Task 0 will find its file already present, re-verify it
and print `skip` — that is the resumability working, not an error.

---

## 4. Submit the two missing North Atlantic months

`01_download.sbatch` spans 1979-01 … 2020-12, so index 0 = Jan 1979 (already
downloaded), **1 = Feb 1979**, **11 = Dec 1979**.

```bash
sbatch --array=1,11 jobs/01_download.sbatch
```

DJF 1979 = Jan + Feb + Dec 1979 (Prosser's record starts 1 Jan 1979, so there
is no Dec 1978 to use). 90 days = 2160 hours — the denominator in §6.

---

## 5. Check on it

```bash
squeue -u $USER
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem
ls -lh $BASE/raw/global/ $BASE/raw/north_atlantic/
du -sh $BASE/raw/global/            # expect ~17 GB when all 12 are in
grep -il "error\|Traceback\|FAILED" logs/global-calib-*.out
```

**Expected cost.** Reference point: a full NA month (5,208 fields, 363 MB) took
6 min 55 s. Each of these is 672 fields but ~1.4 GB — fewer fields to retrieve,
more bytes to move. Budget 5–15 minutes each; at `%4`, all twelve in roughly
20–45 minutes if CDS is quiet, a few hours if it is busy. The 6m55s benchmark
was measured at midnight.

**If a month fails**, it leaves no final file — only a `.tmp`. Re-submitting the
whole array retries exactly the months that failed and skips the ones that
succeeded. Do not use `--force` unless you actually want to re-download.

---

## 6. The acceptance criterion (for when the data is in)

Subset to **36–60°N, 55–10°W** — Prosser's box, from his Figure 2 caption, and
a strict subset of our 30–60°N / 75–0°W download. Apply the year-2000 global
thresholds. Compute the diagnostic-mean exceedance frequency over DJF 1979.

| Severity | Prosser Table 1, DJF 1979 | over 2160 h |
|---|---|---|
| LOG | 128.9 h | **5.97 %** |
| LMOG | 45.6 h | **2.11 %** |
| MOG | 22.3 h | **1.03 %** |
| MSOG | 12.1 h | **0.56 %** |
| SOG | 6.4 h | **0.30 %** |

Read it at three strengths:

1. **Order of magnitude** — MOG near 1 %, not 10 % or 0.01 %. Catches a gross
   calibration or units failure. Failing this means stop.
2. **The ladder shape** — the five frequencies in the ratios
   5.97 : 2.11 : 1.03 : 0.56 : 0.30. This is the strong test: five constraints
   at once on the tail shape of the North Atlantic distribution relative to the
   global one. No accidental agreement fakes it.
3. **Absolute agreement within ~30 %** — a strong result, given the sub-sampled
   calibration and the 175/200/225 vs 188/197/206 hPa level substitution.

**Also run the split-half check**, which costs nothing extra: compute
thresholds on days 1+17 and on days 9+25 separately. Agreement to a few percent
means 48 days was enough sampling. Tens of percent means add days 5/13/21/29.
Raw point counts (4.0 × 10⁸) overstate the effective sample size because
gridpoints within a timestep are spatially correlated — so measure it rather
than assume it.

---

## 7. Known, recorded, not blocking

- **ERA5, not ERA5.1.** Prosser used ERA5.1 for 2000–2006, which corrects a
  lower-stratospheric cold bias. ERA5.1 is a **MARS tape product, not a CDS
  disk dataset**, so switching is not a one-constant change — it means the
  `reanalysis-era5-complete` entry with a different request format and slow tape
  access. ECMWF's own note says behaviour "in most of the troposphere is similar
  to that in ERA5.1". Proceed with ERA5; revisit at publication, not now.
- **Levels.** We evaluate at 200 hPa from a 175/225 stencil; Prosser uses
  197 hPa from 188/206. A narrower stencil gives slightly larger vertical
  derivatives — a near-constant factor, which largely cancels in the percentile
  calibration (STATUS §5.4).
- **F2D loses timesteps at every file boundary.** Its centred material
  derivative drops the first and last step of each file. On a 4-day (32-step)
  calibration file that is 6 %; if the diagnostics are run day-by-day (see §8)
  it becomes 25 %. Acceptable for a percentile, but record which was used.

---

## 8. NOT YET BUILT — the next two pieces

**a) Diagnostics on the global calibration files.** `03_diagnostics.sbatch` is
hardcoded to `raw/north_atlantic/era5_na_${MONTH}.grib` and needs a global
variant.

There is a sizing problem to settle first. One NA month is
121 × 301 × 248 = 9.0 × 10⁶ points and peaked at **MaxRSS 14.8 GB**. One global
4-day file is 721 × 1440 × 32 = **3.3 × 10⁷ points — 3.7× larger**, so expect
roughly **55 GB peak**. That will not fit the current `--mem=24G`, and may not
fit `defq` at all. Two options:

  - process **one day at a time** (721 × 1440 × 8 = 8.3 × 10⁶ points, about the
    size of an NA month, ~15 GB) → 48 array tasks on `defq`, at the cost of F2D
    losing 2 of every 8 steps; or
  - keep the 4-day files and run on `defq-fat`.

**b) The calibration + comparison script.** Loads the 12 calibration files,
computes thresholds with the split-half check, writes them to
`$BASE/calibration/thresholds_<date>.json` via `calibration.save_thresholds`,
applies them to NA DJF 1979, subsets to Prosser's box, and prints §6's table.

Neither blocks the download. Start §3 and §4 now; these can be written while
CDS works.
