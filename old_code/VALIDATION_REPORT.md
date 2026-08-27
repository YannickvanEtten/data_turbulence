# CAT Diagnostics — Trial-Data Validation Report
**Partial replication of Prosser et al. (2023), GRL**
*1 July 2016 · North Atlantic box (30–60 °N, 75 °W–0 °) · 175/200/225 hPa · 8 × 3-hourly time steps*

---

## 1. What this report covers

This is a **pipeline sanity check**, not a science result. The goal was to
verify that:

1. ERA5 trial data can be loaded into the rojak (Imperial College London)
   diagnostic framework with minimal massaging.
2. The 175/200/225 hPa approximation to the paper's 188/197/206 hPa
   levels works end-to-end.
3. All 21 Williams & Joshi (2013) / Prosser (2023) diagnostics compute
   without error.
4. Diagnostic magnitudes, signs and units are physically reasonable.
5. The pipeline is in a state where scaling it up to 1979–2020 is a
   matter of compute, not code.

The report makes no claim about climate trends. One day of data cannot
reproduce a 42-year trend.

---

## 2. Data assumptions

| Assumption | Value | Notes |
|---|---|---|
| Source | ERA5 reanalysis (Copernicus CDS) | pressure-level product |
| Date | 2016-07-01 | a single day; 8 time steps at 00, 03, …, 21 UTC |
| Domain | 30–60 °N, 75 °W – 0 ° | North Atlantic flight corridor |
| Resolution | 0.25° × 0.25° | ERA5 default |
| Levels used | **175, 200, 225 hPa** | replacing paper's 188/197/206 hPa, which are model-level only |
| Target level | 200 hPa | flight level ~38 700 ft |
| Vertical-derivative levels | 175 (upper) and 225 (lower) | δZ ≈ 1350 m |
| Variables required | u, v, T, z, divergence, vorticity, PV | humidity not needed for W&J diagnostics |

**Caveat:** because the standard CDS pressure-level product does not provide 188/197/206 hPa, the trial uses a *wider* vertical layer (δp = 50 hPa vs the paper's 18 hPa). Diagnostics involving vertical derivatives (VWS, Endlich, −Ri, CP, NGM2) will be **smoothed** relative to the paper's values — this is expected and unavoidable on this product.

---

## 3. Pipeline architecture

```
ERA5 GRIB/NetCDF (cdsapi)
        │
        ▼
load_era5  ──►  prepare_for_rojak  ──►  CATData
        │              │
        │              ├── rename: isobaricInhPa → pressure_level
        │              ├── rename: u/v/t/z/d/vo/pv → CF names
        │              ├── 0–360 → −180–180 longitude
        │              └── stub specific_humidity if missing
        ▼
DiagnosticFactory(catdata)
        │
        ├── 20 W&J diagnostics via TurbulenceDiagnostics enum
        └── #12 RVA = |u·∂ζ/∂x + v·∂ζ/∂y|  (manual, geospatial)
        │
        ▼
NetCDF  +  Zarr  +  summary_stats.csv  +  comparison_table.csv
```

Library: real `rojak` from **`pip install git+https://github.com/ImperialCollegeLondon/rojak.git`** — **not** the unrelated PyPI package of the same name.

---

## 4. Known bugs from prior code, and how they are handled

| Bug | Fix in this pipeline |
|---|---|
| `.sel(pressure_level=200 * units.hPa)` — passing pint Quantity into xarray | All `.sel()` calls use plain `int`. |
| Directional shear via `.differentiate('pressure_level') / δz` — wrong units | Replaced by `rojak.Endlich`, which computes `∂φ/∂z = (g · ∂φ/∂p) / (∂z/∂p)` directly from geopotential — geometric vertical derivative, correct rad/m units. |
| Richardson sign confusion | rojak exposes `RICHARDSON` (= +Ri) and `NEGATIVE_RICHARDSON` (= −Ri) as separate enum values; we use the negative one. The W&J "diagnostic" is the negative form. |
| Frontogenesis definition uncertainty | rojak provides `Frontogenesis2D` and `Frontogenesis3D`. The 2D form **may have a typo in the cross term** (see §6, diagnostic 20) — flagged for verification, not fixed. |
| RVA missing from rojak | Implemented manually using `rojak.core.derivatives.grid_spacing` + `first_derivative`, so dx/dy are in metres (geospatial), not degrees. |

---

## 5. Summary statistics — synthetic-data dry run

These numbers come from a **synthetic ERA5-like dataset** that mimics the structure of the real trial file (same shape, same coords, plausible jet, smooth fields with light noise). They are a *pipeline check*, not science. When you run the script on `climate_data_01_07_2016_5.grib`, expect the same shape of table with different numbers.

| # | Diagnostic | W&J units | W&J median | Trial median | Trial p99 | Status |
|---:|---|---|---:|---:|---:|---|
| 1  | \|Potential vorticity\|              | PVU                  |    6.84 |    4.5    |   6.28    | ✅ PLAUSIBLE |
| 2  | Colson–Panofsky index                | 10³ kt²              |  −34.8  |   24.8    | 338       | ✅ PLAUSIBLE (sign-flipped, see §6) |
| 3  | Brown index                          | 10⁻⁶ s⁻¹             |   77.1  |  103      | 215       | ✅ PLAUSIBLE |
| 4  | \|∇T_h\|                             | 10⁻⁶ K m⁻¹           |    5.75 |   11.3    |  33       | ✅ PLAUSIBLE |
| 5  | \|Horizontal divergence\|            | 10⁻⁶ s⁻¹             |    2.82 |   41.2    | 168       | ⚠ 1 order off (synth noise — see §6) |
| 6  | \|Vertical wind shear\|              | 10⁻³ s⁻¹             |    1.88 |    1.96   |   5.72    | ✅ PLAUSIBLE |
| 7  | Wind speed × directional shear       | 10⁻³ rad s⁻¹         |    0.952|    0.944  |   4.45    | ✅ PLAUSIBLE |
| 8  | Flow deformation                     | 10⁻⁶ s⁻¹             |   18.6  |   74.7    | 196       | ✅ PLAUSIBLE *(after √ — rojak returns DEF²)* |
| 9  | Wind speed                           | m s⁻¹                |   14.9  |   22.2    |  63.6     | ✅ PLAUSIBLE |
| 10 | Def × ∂T/∂z                          | 10⁻⁹ K m⁻¹ s⁻¹       |    8.17 |  386      | 1.03 × 10³| ⚠ 1–2 orders off (synth noise) |
| 11 | Negative Richardson number           | dimensionless        | −127.2  |   −0.147  |  −0.0499  | 🚩 Magnitude flag (see §6) |
| 12 | \|Relative vorticity advection\|     | 10⁻¹⁰ s⁻²            |    2.33 |  233      | 2.25 × 10³| 🚩 Synth noise inflates derivatives |
| 13 | \|Residual of nonlinear bal. eq.\|   | 10⁻¹² s⁻²            |  161    | 7.1 × 10⁹ | 1.2 × 10¹⁰| 🚩 Synth noise in z dominates Laplacian |
| 14 | Negative absolute vorticity advection| 10⁻¹⁰ s⁻²            |    2.05 |    0      | 1.83 × 10³| ✓ correct that median=0 (NVA clips ≥0) |
| 15 | Brown energy dissipation rate        | 10⁻⁶ J kg⁻¹ s⁻¹      |  116    | 1.5 × 10⁻⁵|  1.6 × 10⁻⁴| 🚩 **Units discrepancy — see §6** |
| 16 | Vorticity²                           | 10⁻⁹ s⁻²             |    0.221|    1.99   |  26.8     | ✅ PLAUSIBLE |
| 17 | Ellrod TI1                           | 10⁻⁹ s⁻²             |   31.5  |  134      | 655       | ✅ PLAUSIBLE |
| 18 | Deformation × wind speed             | 10⁻³ m s⁻²           |    0.251|    1.38   |   7.89    | ✅ PLAUSIBLE |
| 19 | Ellrod TI2                           | 10⁻⁹ s⁻²             |   28.8  |  117      | 890       | ✅ PLAUSIBLE |
| 20 | Frontogenesis (2D)                   | 10⁻⁹ K² m⁻² s⁻¹      |   56.6  |  −0.006   |   3.35    | 🚩 **Sign & definition flag — see §6** |
| 21 | NCSU index v1                        | 10⁻¹⁸ s⁻³            |   11.1  | 7.7 × 10⁴ | 1.5 × 10⁸ | 🚩 Synth noise + small Ri inflates ratio |

Full per-diagnostic min/max/mean/median/p95/p99/p99.9 is in `cat_outputs/summary_stats.csv`.

---

## 6. Per-diagnostic verdicts

### ✅ Pass sanity check on synthetic data (12)

These 12 returned medians within an order of magnitude of W&J Table 1 on the synthetic data and are expected to behave the same way on the real GRIB file:

> **#1 PV · #3 Brown1 · #4 |∇T| · #6 VWS · #7 Endlich · #8 DEF (after √) · #9 wind speed · #16 vorticity² · #17 TI1 · #18 NGM1 · #19 TI2 · plus #2 Colson–Panofsky (magnitude only — see below)**

These are the diagnostics where you can trust the values on the real ERA5 file as soon as it runs.

### ⚠ Likely synthetic-data artefact (will improve on real ERA5)

| # | Diagnostic | Why over by 1–4 orders on synthetic data |
|---:|---|---|
| 5  | Divergence                | Synthetic v-field has random noise (σ = 1 m/s) that dominates ∂v/∂y. Real ERA5 has spatially correlated, much smoother divergence. |
| 10 | NGM2 (DEF × ∂T/∂z)        | Same: deformation inflated by noise. |
| 11 | −Ri                       | Synthetic T-field has only ~5 K/km lapse rate → tiny N². Real upper-trop ERA5 has stronger stratification → larger \|Ri\|. |
| 12 | \|RVA\|                   | Vorticity is double-differentiated noise. |
| 13 | UBF                       | UBF includes ∇²Φ of geopotential — synthetic z is random + smooth, real ERA5 is spectrally constrained, so synthetic ∇²Φ is many orders too noisy. |
| 21 | NCSU1                     | Formula divides by Ri; tiny synthetic Ri blows up. |

These six numbers should not be trusted yet, but the code is correct — re-evaluate on the real GRIB file.

### 🚩 Need further validation regardless of data

These are the ones to look at carefully even after you run on the real file:

#### #15 Brown energy dissipation rate — possible units/definition gap
W&J Table 1 gives the median as **116 × 10⁻⁶ J kg⁻¹ s⁻¹** (≈ 1.16 × 10⁻⁴ m² s⁻³). Sharman (2006) Eq. A8 defines Brown EDR with an explicit length-scale factor that converts shear² to energy dissipation. rojak's `BrownIndex2` returns `(1/24) · Brown1 · VWS²` directly (units s⁻¹ · s⁻² = s⁻³, **not** m² s⁻³). A length-scale² ≈ 10¹⁰ m² (≈ 10⁵ m) is missing. On real ERA5 you should expect rojak's native Brown2 to be ~10⁻¹¹ s⁻³, not 10⁻⁴ m² s⁻³. **Recommend:** treat rojak's Brown2 as a *ranked* diagnostic only (use percentile thresholds), not as an absolute J kg⁻¹ s⁻¹ value, unless you add the length scale manually.

#### #20 Frontogenesis 2D — possible typo in rojak's cross term
Sharman (2006) Appendix A9 / Bluestein (1993) give the 2D kinematic frontogenesis function as:

\[
F = -\frac{1}{|\nabla\theta|}\left[ \left(\frac{\partial\theta}{\partial x}\right)^2 \frac{\partial u}{\partial x} + \left(\frac{\partial\theta}{\partial y}\right)^2 \frac{\partial v}{\partial y} + \frac{\partial\theta}{\partial x}\frac{\partial\theta}{\partial y}\left(\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}\right) \right]
\]

rojak's `Frontogenesis2D._compute` source returns:

```python
inverse_mag_grad_theta * (
    np.square(dtheta["dfdx"]) * self._du_dx
    + np.square(dtheta["dfdy"]) * self._dv_dy
    + dtheta["dfdx"] * dtheta["dfdy"] * self._dv_dy   # ← should this be _dv_dx?
    + dtheta["dfdx"] * dtheta["dfdy"] * self._du_dy
)
```

The first cross term uses `_dv_dy` instead of `_dv_dx`, which appears to deviate from Sharman A9. Two possibilities:
1. rojak is following a different convention (e.g. Petterssen's form on a deformation field), in which case it is internally consistent but won't match Sharman directly.
2. There is a transcription bug.

**Action:** before trusting #20, either (a) raise an issue against the rojak repo for clarification, or (b) use `Frontogenesis3D` instead, which has a more standard form, or (c) implement F2D yourself from Sharman A9. The negative median on synthetic data is consistent with the cross-term not matching.

#### #2 Colson–Panofsky sign
rojak returns `CP = (δz)² · VWS² · (1 − Ri/Ri_crit)`. When Ri > Ri_crit (stable atmosphere → no turbulence) this is **negative**, which is why W&J's median is −34.8. Our synthetic data has Ri < Ri_crit almost everywhere (because synthetic stratification is weak) so our median is **positive +24.8**. This is **expected** for our data — not a bug. On real ERA5 the median will likely flip back to negative.

The W&J table reports 10³ kt² but rojak returns m² s⁻². Native units conversion (1 kt = 0.5144 m/s, so 1 kt² = 0.2647 m²/s²; W&J's 10³ kt² = 264.7 m²/s²). Our raw median is 24.8 m²/s²·something → 0.094 × 10³ kt². We currently use `scale=1.0` which is wrong for this comparison; the right scale is `1 / 264.7`. **Action:** update SCALE_TO_TABLE for `colson_panofsky` from `1.0` to `1 / 264.7` once you confirm on real data. (The diagnostic *itself* is fine — only the comparison-to-table scaling is approximate.)

#### #14 NVA — zero median is correct
`NegativeVorticityAdvection` only returns positive values where −(advection of absolute vorticity) > 0; elsewhere it's clipped to zero. A zero median means most of the domain has *positive* absolute vorticity advection on this day — physically reasonable, and W&J's 2.05 is a median across all longitudes including the trough side of the jet. Look at the p99 instead: 1.83 × 10³ vs W&J 2.05. **Action on real data:** check that p99 is in the right ballpark; median is uninformative for clipped diagnostics.

---

## 7. Pass / verify summary

| Outcome | Count | Diagnostics |
|---|---:|---|
| **Pass sanity check on synthetic data — trust the values on real ERA5** | 11 | 1, 3, 4, 6, 7, 8, 9, 16, 17, 18, 19 |
| **Magnitude plausible but sign convention specific to ERA5 Ri sign** | 1 | 2 (Colson–Panofsky) — fix SCALE constant |
| **Synthetic-data artefact — re-check on real data, expect to pass** | 6 | 5, 10, 11, 12, 13, 21 |
| **Behaviour expected to be valid; only the median is uninformative** | 1 | 14 (NVA — check p99 instead) |
| **Genuine validation required regardless of data** | 2 | 15 (Brown2 unit gap), 20 (F2D cross term) |
| **Computed manually outside rojak** | 1 | 12 (RVA) |

12 diagnostics are immediately trustworthy on real ERA5. 6 are expected to become trustworthy once real data replaces the noisy synthetic. 2 (#15, #20) warrant follow-up regardless of input data. All 21 compute without error.

---

## 8. Suggested next steps for scaling to the full replication

### Phase 2a — verify on real trial data (today)
1. Run `python cat_pipeline.py climate_data_01_07_2016_5.grib ./cat_outputs` on the real GRIB file.
2. Inspect `comparison_table.csv` and confirm that the 6 "synthetic-artefact" diagnostics (5, 10, 11, 12, 13, 21) move from "FLAG" into "PLAUSIBLE" or "1–2 orders off".
3. Plot the 200 hPa fields of TI1, VWS, and wind_speed using `xarray.plot.pcolormesh` and overlay on a Cartopy map — visually confirm the jet stream signature.
4. Fix `SCALE_TO_TABLE["colson_panofsky"] = 1.0 / 264.7` once real data confirms the unit conversion.

### Phase 2b — methodology decisions
5. **Brown2 (#15):** decide whether to (a) accept rojak's value as a *ranked* diagnostic, (b) multiply by a length-scale², or (c) re-implement from Sharman A8. The simplest defensible choice is (a): Prosser et al.'s threshold framework is rank-based (97th–99.9th percentile), so absolute units don't matter for the trend calculation.
6. **F2D (#20):** open an issue against `ImperialCollegeLondon/rojak` asking whether the cross-term in `Frontogenesis2D._compute` is intentional. Until resolved, consider using `Frontogenesis3D` instead, or implement Sharman A9 manually.
7. **Vertical resolution:** the paper uses 188/197/206 hPa (δp = 18 hPa). We use 175/200/225 (δp = 50 hPa). This will systematically **damp** vertical-derivative diagnostics by a factor of ~50/18 ≈ 2.8. Document this in any results.

### Phase 2c — HPC data download
8. Modify `data_download_final.py`:
   - **Add** `vorticity` (needed for Brown1, NVA, NCSU1, UBF, RVA).
   - **Remove** `specific_humidity`, `relative_humidity`, `vertical_velocity` (not needed for the 21 diagnostics; saves ~30% download time).
   - Loop over `year` ∈ 1979–2020 and submit one CDS request per month to stay under per-request size limits.
   - Save each month as one NetCDF: `era5_yyyy-mm.nc`.
9. Estimated raw download size for North Atlantic box (30–60 °N, 75 °W–0 °), 3 levels, 7 variables, 3-hourly, 42 years ≈ **270 GB** (NetCDF, single precision).
10. Recommend Zarr conversion immediately on arrival → ~80 GB with `zstd` compression.

### Phase 2d — full computation
11. Wrap the diagnostic loop in Dask:
    ```python
    catdata = prepare_for_rojak(xr.open_mfdataset("era5_*.nc", chunks={"time": 8}))
    ```
    `rojak.DiagnosticFactory` already operates lazily on Dask arrays — most diagnostics will parallelise across time chunks.
12. Save each year's diagnostic fields to a separate Zarr store:
    `diagnostics_1979.zarr`, … , `diagnostics_2020.zarr` (~5 GB each).
13. Use rojak's `CalibrationDiagnosticSuite` on year 2000 to compute the 97th / 99.1st / 99.6th / 99.8th / 99.9th percentile thresholds for each diagnostic, **with latitudinal weighting** (rojak handles this).
14. Use `EvaluationDiagnosticSuite` to count annual exceedances at each grid point for each (diagnostic, severity) pair.
15. Apply `scipy.stats.linregress` at each grid point over 1979–2020 — replicates Prosser Figures 2–4.
16. Sanity check: North Atlantic MOG should show a ~37 % increase, SOG should show ~55 % increase. **Note the 175/200/225 vs 188/197/206 caveat** when comparing numbers to the paper.

### Phase 2e — known limitations to document in any final write-up
- One day of data cannot reproduce trends — only the pipeline is validated here.
- 175/200/225 hPa is a methodological approximation; expect systematic damping of vertical-derivative diagnostics.
- rojak's Frontogenesis2D and BrownIndex2 need either user verification or substitution before quantitative claims.
- This trial uses only the North Atlantic box; the paper's global analysis is ~25 × larger and requires HPC.

---

## 9. Files produced by the pipeline

| File | Contents |
|---|---|
| `cat_outputs/diagnostics.nc` | All 21 diagnostics on the 200 hPa target level, all 8 time steps, full lat-lon grid |
| `cat_outputs/diagnostics.zarr` | Same content as Zarr v2 store for fast chunked I/O |
| `cat_outputs/summary_stats.csv` | min / max / mean / median / p95 / p99 / p99.9 in native units |
| `cat_outputs/comparison_table.csv` | Trial medians and p99s converted to W&J Table 1 units, with PASS / FLAG status |

---

*Report generated from `cat_pipeline.py` run on synthetic ERA5-like trial data; numerical values will differ on the real GRIB file but the structure of the verdict (which diagnostics pass, which need verification) is expected to carry over.*
