"""
ada/_checkutil.py
=================
Shared helpers for the FORMULA_AUDIT.md verification scripts
(`check_composition.py`, `check_shape_ratios.py`, `check_f2d_variants.py`,
`check_ubf.py`).

Deliberately small, and deliberately NOT a refactor of the equivalent helpers
inside `ada/trend_check.py`. Those are load-bearing for a result already
recorded in STATUS.md §12; duplicating four short functions costs less than
risking a silent change to the numbers in that section.
"""
from __future__ import annotations

import importlib.util
import resource
import sys
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[1]
BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")

# Williams (2017) Table 2 / Williams & Joshi (2013) Table 1 box.
# NOT Prosser's box -- see check_shape_ratios.py.
WILLIAMS_BOX = dict(lat=(50.0, 75.0), lon=(-60.0, -10.0))
PROSSER_BOX = dict(lat=(36.0, 60.0), lon=(-55.0, -10.0))


def load_module(name: str, filename: str):
    """Load a module whose filename is not a legal Python identifier.

    Registering in sys.modules BEFORE exec_module is mandatory: 2_diagnostics.py
    uses `from __future__ import annotations`, so @dataclass resolves its field
    annotations through sys.modules[cls.__module__]. Leaving that None is the
    Phase-0 crash in STATUS.md §6.
    """
    spec = importlib.util.spec_from_file_location(name, REPO / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def peak_rss_gb() -> float:
    """Peak RSS so far, GB. Linux reports ru_maxrss in kB.

    Every driver prints its own, because SLURM job-step accounting is disabled
    on this cluster and sacct/sstat return empty MaxRSS (STATUS.md §11.9).
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def subset_box(ds, lat: tuple, lon: tuple):
    """Subset by bounds regardless of coordinate order.

    ERA5 latitudes descend, so a naive slice(50, 75) returns an EMPTY selection
    with NO error -- which would make every number downstream meaningless
    rather than obviously wrong. Raise instead.
    """
    la, lo = ds["latitude"].values, ds["longitude"].values
    lat_sl = slice(lat[1], lat[0]) if la[0] > la[-1] else slice(lat[0], lat[1])
    lon_sl = slice(lon[1], lon[0]) if lo[0] > lo[-1] else slice(lon[0], lon[1])
    out = ds.sel(latitude=lat_sl, longitude=lon_sl)
    if out.sizes["latitude"] == 0 or out.sizes["longitude"] == 0:
        raise ValueError(f"box lat={lat} lon={lon} selected nothing from "
                         f"lat[{la[0]}..{la[-1]}] lon[{lo[0]}..{lo[-1]}]")
    return out


def cos_phi_weights(da: xr.DataArray) -> np.ndarray:
    """cos(latitude) weights broadcast to `da`'s shape, flattened.

    The same weighting Prosser applies and calib_weighted_percentile
    implements; a regular lat-lon grid over-samples the poles and an
    unweighted percentile is biased toward them.
    """
    w = xr.DataArray(np.cos(np.deg2rad(da["latitude"].values)),
                     coords={"latitude": da["latitude"]}, dims=("latitude",))
    return np.asarray(w.broadcast_like(da).values, dtype=np.float64).ravel()


def rel_error(got: np.ndarray, want: np.ndarray, floor: float | None = None):
    """Relative error |got-want| / max(|want|, floor), NaNs dropped.

    `floor` guards the case that makes naive relative error useless: a
    denominator passing through zero. CALIBRATION_REFERENCE.md §10.4 records
    the split-half check printing a 92 % "failure" that was entirely this --
    colson_panofsky's LMOG threshold sits essentially at zero, so a relative
    difference across it is meaningless. Default floor is a small fraction of
    the field's own robust scale rather than an absolute constant.
    """
    got = np.asarray(got, dtype=np.float64).ravel()
    want = np.asarray(want, dtype=np.float64).ravel()
    ok = np.isfinite(got) & np.isfinite(want)
    got, want = got[ok], want[ok]
    if got.size == 0:
        return np.array([])
    if floor is None:
        scale = np.median(np.abs(want))
        floor = 1e-6 * scale if scale > 0 else 1e-30
    return np.abs(got - want) / np.maximum(np.abs(want), floor)


def summarise_error(err: np.ndarray) -> dict:
    if err.size == 0:
        return dict(n=0, median=float("nan"), p99=float("nan"), max=float("nan"))
    return dict(n=int(err.size),
                median=float(np.median(err)),
                p99=float(np.quantile(err, 0.99)),
                max=float(err.max()))


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rank correlation without scipy.

    scipy is not in pixi.toml and adding a dependency to run one correlation
    would change the pinned environment that STATUS.md §10.4 treats as the
    reproducibility record. Ranking with argsort costs three lines.
    """
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if a.size < 3:
        return float("nan")

    def ranks(x):
        order = np.argsort(x, kind="stable")
        r = np.empty(x.size, dtype=np.float64)
        r[order] = np.arange(x.size, dtype=np.float64)
        return r

    ra, rb = ranks(a), ranks(b)
    ra -= ra.mean()
    rb -= rb.mean()
    denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / denom) if denom > 0 else float("nan")


def daily_mean(ds, time_dim: str = "time"):
    """Average to daily means.

    Williams (2017) §2 states the diagnostics there are "calculated from daily
    mean temperature and wind fields", and Williams & Joshi (2013) Table 1 is
    from the same data. Our fields are 3-hourly instantaneous, and averaging
    thins a tail, so comparing an instantaneous p97 against a daily-mean p97
    would bias the comparison in a known direction. This removes that
    difference rather than arguing about its size.

    NOTE this averages the DIAGNOSTIC, whereas Williams averages the WIND AND
    TEMPERATURE and then computes the diagnostic. For a non-linear diagnostic
    those are not the same operation. It is the closer of the two available
    approximations and the residual difference is a caveat, not a correction --
    see check_shape_ratios.py's output note.
    """
    return ds.resample({time_dim: "1D"}).mean()
