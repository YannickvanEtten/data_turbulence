#!/usr/bin/env python
"""
ada/check_production_complete.py
=================================
Is the 42-year production series actually complete, and is it ONE consistent
convention throughout -- including the thresholds it will be compared against?

    pixi run python ada/check_production_complete.py
    pixi run python ada/check_production_complete.py --domain both

WHY THIS EXISTS
---------------
`ada/diagnostics_global.py --skip-if-matching` (RUNBOOK_audit_checks.md §5.2)
already refuses to silently reuse output written under an older convention for
a SINGLE file. This script asks the question at the level the project actually
cares about: across all 504 North Atlantic months plus the 12 global
calibration months, is every store (a) present, (b) a complete zarr write, and
(c) built with the SAME f2d_variant and deformation_convention as every other
store?

A directory listing that is 100% present says nothing about (c). The failure
mode this guards against is explicit in diagnostics_global.py's own docstring:
27 DJF months were written before the 2026-08-29/08-30 convention fixes, and a
naive "does the file exist" check would hand back a 504-month series that is
27 months of one convention and 477 of another -- a real, plausible-looking
dataset quietly built from two definitions of the same variable.

THE SECOND, MORE EASILY MISSED PROBLEM
---------------------------------------
`deformation`'s un-squaring is a monotone transform of a non-negative field,
so by STATUS.md §5.5 it cannot silently corrupt an EXCEEDANCE field computed
against thresholds that were derived under the same convention -- but it can
absolutely corrupt one derived under a DIFFERENT convention. A threshold
calibrated as a percentile of DEF^2 is not the square root away from being a
valid threshold for DEF; applying it to un-squared production data compares
the wrong quantity outright. And `f2d`'s variant change (A -> C) is not even a
monotone transform of the same field -- FORMULA_AUDIT.md §10.4 measured
Spearman rho(A, D) = -0.94 and an exact-zero tail overlap between A and B, so a
threshold calibrated under one variant has no defined relationship to data
computed under another.

So this script also checks whether `derived/global/` (the calibration input)
carries the SAME provenance as `derived/north_atlantic/`, and whether the
newest `calibration/thresholds_*.json` was written on or after the date the
conventions were fixed (2026-08-30, RUNBOOK_audit_checks.md §1.2). The JSON
schema (calibration.py) does not currently record f2d_variant or
deformation_convention in its own provenance block -- so a date check is what
is available; if this fires, the honest fix is to add those two fields to
`calibration.save_thresholds`'s provenance going forward, not to trust the
date alone indefinitely.

READING THE OUTPUT
-------------------
Exit 0 only if every expected NA and global-calibration store is present,
complete, and carries matching provenance, AND the thresholds file postdates
the convention fix. Anything else is a stop sign, not a warning: every check
downstream of this one (calibration comparison, trend fit, per-diagnostic
breakdown) assumes the answer here is yes.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

BASE = Path("/scistor/SBE-EDS-ClimateKoopman/yen230")

# The date the deformation and f2d fixes landed (RUNBOOK_audit_checks.md §1,
# FORMULA_AUDIT.md §10.4's "superseded 2026-08-30" note). A thresholds file
# created before this cannot have been calibrated under the current
# conventions, whatever its provenance block does or does not say.
CONVENTION_FIX_DATE = date(2026, 8, 30)

EXPECTED_F2D_VARIANT = "C"
EXPECTED_DEF_CONVENTION = "DEF (un-squared, Sharman A17)"


def read_zattrs(store: Path) -> dict | None:
    """Read a zarr v2 store's top-level attrs without opening it with xarray.

    Same technique as ada/diagnostics_global.py's existing_output_matches() --
    deliberately duplicated rather than imported, so this checker has no
    dependency on the production driver and cannot be silently broken by an
    unrelated change to it.
    """
    zattrs = store / ".zattrs"
    if not zattrs.exists():
        return None
    try:
        return json.loads(zattrs.read_text())
    except (OSError, ValueError):
        return None


def store_status(store: Path, expected_variant: str, expected_def: str) -> str:
    """One of: MISSING, INCOMPLETE, NO_PROVENANCE, WRONG_VARIANT, WRONG_DEF, OK."""
    if not store.exists():
        return "MISSING"
    complete = (store / ".zmetadata").exists() or (store / "zarr.json").exists()
    if not complete:
        return "INCOMPLETE"
    attrs = read_zattrs(store)
    if attrs is None:
        return "NO_PROVENANCE"
    variant = attrs.get("f2d_variant")
    defo = attrs.get("deformation_convention")
    if variant is None or defo is None:
        return "NO_PROVENANCE"
    if variant != expected_variant:
        return f"WRONG_VARIANT({variant})"
    if defo != expected_def:
        return f"WRONG_DEF({defo})"
    return "OK"


def sweep_north_atlantic(base: Path, start_year: int, end_year: int,
                          variant: str, defo: str) -> dict[str, list[str]]:
    by_status: dict[str, list[str]] = {}
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            name = f"diagnostics_na_{year}-{month:02d}.zarr"
            store = base / "derived" / "north_atlantic" / name
            status = store_status(store, variant, defo)
            by_status.setdefault(status, []).append(name)
    return by_status


def sweep_global_calib(base: Path, variant: str, defo: str) -> dict[str, list[str]]:
    by_status: dict[str, list[str]] = {}
    for month in range(1, 13):
        name = f"diagnostics_glob_2000-{month:02d}.zarr"
        store = base / "derived" / "global" / name
        status = store_status(store, variant, defo)
        by_status.setdefault(status, []).append(name)
    return by_status


def check_thresholds(base: Path) -> tuple[bool, str]:
    found = sorted((base / "calibration").glob("thresholds_*.json"))
    if not found:
        return False, "no thresholds_*.json found in calibration/"
    newest = found[-1]
    try:
        payload = json.loads(newest.read_text())
    except (OSError, ValueError) as exc:
        return False, f"{newest.name}: unreadable ({type(exc).__name__})"

    created_str = payload.get("created", "")
    try:
        created = date.fromisoformat(created_str[:10])
    except ValueError:
        return False, f"{newest.name}: unparseable created={created_str!r}"

    notes = payload.get("provenance", {}).get("notes", "")
    if created < CONVENTION_FIX_DATE:
        return False, (
            f"{newest.name} was created {created.isoformat()}, BEFORE the "
            f"convention fix on {CONVENTION_FIX_DATE.isoformat()}. This "
            f"thresholds file was almost certainly calibrated against "
            f"deformation-squared and/or f2d variant A data. Re-run "
            f"jobs/04_diagnostics_global.sbatch (with the corrected "
            f"defaults) and jobs/06_calibration_check.sbatch before trusting "
            f"any exceedance or trend number computed against it."
        )
    return True, f"{newest.name}, created {created.isoformat()} (notes: {notes!r})"


def report(title: str, by_status: dict[str, list[str]], expected_total: int) -> bool:
    print("=" * 78)
    print(title)
    print("=" * 78)
    ok = len(by_status.get("OK", []))
    print(f"   expected {expected_total}, OK {ok}, "
          f"problems {expected_total - ok}")
    for status, names in sorted(by_status.items()):
        if status == "OK":
            continue
        preview = ", ".join(names[:8])
        more = f" (+{len(names) - 8} more)" if len(names) > 8 else ""
        print(f"   {status:<20} {len(names):>4}   {preview}{more}")
    return ok == expected_total


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=str(BASE))
    ap.add_argument("--start-year", type=int, default=1979)
    ap.add_argument("--end-year", type=int, default=2020)
    ap.add_argument("--domain", choices=["na", "global", "both"], default="both")
    ap.add_argument("--expect-f2d-variant", default=EXPECTED_F2D_VARIANT)
    ap.add_argument("--expect-deformation", default=EXPECTED_DEF_CONVENTION)
    args = ap.parse_args()

    base = Path(args.base)
    all_ok = True

    if args.domain in ("na", "both"):
        n_expected = (args.end_year - args.start_year + 1) * 12
        by_status = sweep_north_atlantic(base, args.start_year, args.end_year,
                                         args.expect_f2d_variant,
                                         args.expect_deformation)
        all_ok &= report(
            f"NORTH ATLANTIC — derived/north_atlantic/ "
            f"({args.start_year}-{args.end_year})", by_status, n_expected)
        print()

    if args.domain in ("global", "both"):
        by_status = sweep_global_calib(base, args.expect_f2d_variant,
                                       args.expect_deformation)
        all_ok &= report("GLOBAL CALIBRATION — derived/global/ (year 2000)",
                         by_status, 12)
        print()

    print("=" * 78)
    print("THRESHOLDS FILE — calibration/thresholds_*.json")
    print("=" * 78)
    thresh_ok, thresh_msg = check_thresholds(base)
    print(f"   {'OK' if thresh_ok else '!! PROBLEM'}: {thresh_msg}")
    all_ok &= thresh_ok

    print("\n" + "=" * 78)
    if all_ok:
        print("ALL CLEAR — production series and thresholds are complete and")
        print("consistent. Safe to run ada/full_trend_check.py and treat its")
        print("output as the real 42-year result, not a provisional one.")
    else:
        print("NOT CLEAR — see the problems listed above. Nothing computed from")
        print("this data should be reported as a validated result until every")
        print("section above reads OK. Likely fix: resubmit jobs/14 for any")
        print("NA gaps, and jobs/04 + jobs/06 if the global/thresholds section")
        print("or the thresholds-file check failed.")
    print("=" * 78)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
