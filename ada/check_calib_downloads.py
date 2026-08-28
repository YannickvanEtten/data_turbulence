#!/usr/bin/env python
"""
ada/check_calib_downloads.py
============================
Verify that every file the calibration check needs is present and correct,
BEFORE spending compute on diagnostics.

    pixi run python ada/check_calib_downloads.py

Checks, for each expected file:
  * it exists, and its size is sane
  * it opens with cfgrib as a single hypercube
  * 7 variables, 3 pressure levels, the right number of timesteps
  * the horizontal grid matches the domain it claims in its filename

The per-file check is `1_download_hpc.verify_file` -- the SAME function the
downloader ran before renaming the .tmp. Re-running it here is not redundant:
it catches a file that was truncated or corrupted by the filesystem AFTER the
download passed, and it puts a single "everything is ready" statement on the
record before the next stage starts.

Exit status is 0 only if every expected file passes.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")
REPO = Path(__file__).resolve().parents[1]

sys.path.insert(0, str(REPO))
import download_plan  # noqa: E402


def _load_downloader():
    """Load 1_download_hpc.py, whose name is not a legal Python identifier.

    Registering it in sys.modules is not optional: 2_diagnostics.py uses
    `from __future__ import annotations`, so @dataclass resolves its field
    annotations through sys.modules[cls.__module__]. Leaving that unset is what
    caused the bare "'NoneType' object has no attribute '__dict__'" crash
    recorded in STATUS.md Phase 0.
    """
    path = REPO / "1_download_hpc.py"
    spec = importlib.util.spec_from_file_location("download_hpc", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["download_hpc"] = module
    spec.loader.exec_module(module)
    return module


CALIB_DAYS = [1, 9, 17, 25]

# What the calibration check needs, and why.
EXPECTED = []

# 12 global months, sub-sampled -- the year-2000 calibration domain.
for month in range(1, 13):
    days = [f"{d:02d}" for d in CALIB_DAYS]
    EXPECTED.append({
        "path": BASE / "raw/global" /
                f"era5_glob_2000-{month:02d}_d{'-'.join(days)}.grib",
        "request": download_plan.build_request(2000, month, days=days,
                                               domain="global"),
        "why": f"calibration 2000-{month:02d}",
    })

# 3 North Atlantic months -- DJF 1979, the analysis period for the check.
for (year, month) in [(1979, 1), (1979, 2), (1979, 12)]:
    EXPECTED.append({
        "path": BASE / "raw/north_atlantic" /
                download_plan.month_filename(year, month, "north_atlantic"),
        "request": download_plan.build_request(year, month,
                                               domain="north_atlantic"),
        "why": f"DJF 1979 ({year}-{month:02d})",
    })


def main() -> int:
    dl = _load_downloader()

    print(f"{'file':<44}{'why':<22}{'size':>9}  status")
    print("-" * 92)

    ok = bad = missing = 0
    total_bytes = 0

    for item in EXPECTED:
        path, why = item["path"], item["why"]
        name = path.name

        if not path.exists():
            tmp = path.with_name(path.name + ".tmp")
            note = "MISSING (.tmp present -- download failed)" if tmp.exists() \
                   else "MISSING"
            print(f"{name:<44}{why:<22}{'-':>9}  {note}")
            missing += 1
            continue

        size = path.stat().st_size
        total_bytes += size
        size_gb = size / 1e9

        try:
            problems = dl.verify_file(path, item["request"])
        except Exception as exc:                      # noqa: BLE001
            print(f"{name:<44}{why:<22}{size_gb:>8.2f}G  "
                  f"UNREADABLE: {type(exc).__name__}: {exc}")
            bad += 1
            continue

        if problems:
            print(f"{name:<44}{why:<22}{size_gb:>8.2f}G  FAILED: {problems}")
            bad += 1
        else:
            print(f"{name:<44}{why:<22}{size_gb:>8.2f}G  ok")
            ok += 1

    print("-" * 92)
    print(f"{ok} ok, {bad} failed, {missing} missing "
          f"(of {len(EXPECTED)} expected)")
    print(f"total on disk: {total_bytes / 1e9:.1f} GB")

    if missing:
        print("\nMissing months: re-submit the array. The downloader is "
              "resumable -- months with a passing final file are skipped, so "
              "re-submitting retries exactly the ones that failed:")
        print("    sbatch jobs/02b_download_global_calib.sbatch")
        print("    sbatch --array=1,11 jobs/01_download.sbatch")
    if bad:
        print("\nFailed files are NOT safe to use. Delete each one and "
              "re-download it; do not pass them to the diagnostics.")

    if ok == len(EXPECTED):
        print("\nAll inputs present and verified. Ready for the diagnostics "
              "stage.")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
