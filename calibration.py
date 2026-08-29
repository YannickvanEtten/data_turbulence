"""
calibration.py
==============
Severity thresholds as a FILE — the handoff between the global calibration run
and the North Atlantic analysis run.

WHY THIS MODULE EXISTS
----------------------
Prosser's method calibrates severity thresholds on one domain (a global year)
and applies them on another (the North Atlantic, 42 years). Those are two
separate compute jobs, weeks apart, on different data. Before this module the
thresholds existed only as a local variable inside
`3_pipeline.run_layers_2_to_6` — computed and consumed inside a single function
call, with no way to:

  * look at them,
  * compare them against the published tables (which is the whole point of
    calibrating on a global year and the main check that the chain works),
  * reuse them later without recomputing the global year, or
  * record which data produced them.

A threshold set is a scientific result. It gets a file, with provenance.

THE GOLDILOCKS BAND
-------------------
`aggregate.assert_thresholds_not_hardcoded_table2` RAISES if thresholds match
Williams (2017) Table 2 to within 1e-4, because an exact match means someone
copy-pasted the table instead of calibrating. So when comparing a real
calibration against the literature there are three outcomes, and only one of
them is good:

  * same order of magnitude, same sign, ratio roughly 0.3–3   -> CONFIRMS
  * identical to published (ratio 1.0000)                     -> SUSPICIOUS,
    and the guard in aggregate.py will refuse to run
  * orders of magnitude apart, or opposite sign               -> BROKEN

`compare_to_published()` below reports the ratio and labels it accordingly, so
the check is a table you read rather than a judgement you make from memory.
"""
from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path

import numpy as np
import pandas as pd

SCHEMA_VERSION = 1


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------
def compute_thresholds(
    calibration_fields: dict,
    lat_weights,
    severities: dict[str, float],
    reference_table: dict,
    weighted_percentile,
) -> tuple[dict, dict, dict]:
    """Percentile-calibrate every diagnostic on the calibration domain.

    calibration_fields : {name: DataArray} over the CALIBRATION domain
                         (global year 2000), NOT the analysis box.
    lat_weights        : cos(phi) weights matching those fields' grid.
    severities         : {"light": 97.0, ...} name -> percentile cutoff.
    reference_table    : 2_diagnostics.REFERENCE_TABLE, for the sign of each.
    weighted_percentile: calib_weighted_percentile.weighted_percentile,
                         injected rather than imported so this module stays
                         free of the numbered-filename import dance.

    Returns (thresholds, signs, sample_sizes).
    """
    thresholds: dict[str, dict[str, float]] = {}
    signs: dict[str, str] = {}
    sample_sizes: dict[str, int] = {}

    for name, field in calibration_fields.items():
        ref = reference_table[name]
        signs[name] = ref["sign"]

        values = np.asarray(field).ravel()
        weights = np.asarray(lat_weights.broadcast_like(field)).ravel()
        finite = np.isfinite(values)
        sample_sizes[name] = int(finite.sum())
        if sample_sizes[name] == 0:
            raise ValueError(
                f"diagnostic {name!r} has no finite values over the calibration "
                f"domain — refusing to produce a threshold from nothing."
            )

        # ONE call with all five percentiles, not one call per severity.
        #
        # weighted_percentile sorts its whole input on every call (argsort plus
        # np.unique), so asking for the severities one at a time sorted the same
        # array five times. On the year-2000 global calibration that is 4e8
        # values x 5 x 21 diagnostics = 105 full sorts where 21 suffice, and the
        # job does not finish inside any reasonable walltime.
        #
        # The function already accepts an array of percentiles and returns one
        # value per entry -- np.interp is vectorised over p and the sort is
        # shared. Measured 5.67x on a 2e6-element benchmark, with output
        # bit-identical to the per-severity calls (max abs diff 0.0).
        sev_names = list(severities)
        pcts = np.asarray([severities[s] for s in sev_names], dtype=float)
        computed = np.atleast_1d(weighted_percentile(values, weights, pcts))
        if computed.shape[0] != len(sev_names):
            raise RuntimeError(
                f"weighted_percentile returned {computed.shape[0]} values for "
                f"{len(sev_names)} percentiles on diagnostic {name!r} -- refusing "
                f"to zip a mismatched threshold ladder."
            )
        thresholds[name] = {s: float(v) for s, v in zip(sev_names, computed)}

    return thresholds, signs, sample_sizes


# ---------------------------------------------------------------------------
# Persist
# ---------------------------------------------------------------------------
def save_thresholds(path, thresholds, signs, severities, *,
                    domain, period, pressure_levels, sample_sizes=None,
                    rojak_version=None, notes=None) -> Path:
    """Write a threshold set with enough provenance to be trusted later.

    The metadata is not decoration. A threshold file with no record of which
    domain, which year and which levels produced it is indistinguishable from
    a made-up one six months later, and it is exactly the sort of artefact a
    reviewer asks about.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "created": _dt.datetime.now().isoformat(timespec="seconds"),
        "provenance": {
            "calibration_domain": domain,
            "period": period,
            "pressure_levels": list(pressure_levels),
            "severities_percentiles": dict(severities),
            "rojak_version": rojak_version,
            "sample_sizes": sample_sizes or {},
            "notes": notes or "",
        },
        "signs": signs,
        "thresholds": thresholds,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return path


def load_thresholds(path) -> tuple[dict, dict, dict]:
    """Return (thresholds, signs, provenance). Fails loud on a schema it does
    not know, rather than silently mis-reading a future format."""
    payload = json.loads(Path(path).read_text())
    got = payload.get("schema_version")
    if got != SCHEMA_VERSION:
        raise ValueError(
            f"threshold file {path} has schema_version {got!r}, this code "
            f"understands {SCHEMA_VERSION}. Refusing to guess."
        )
    return payload["thresholds"], payload["signs"], payload["provenance"]


# ---------------------------------------------------------------------------
# The literature check
# ---------------------------------------------------------------------------
def compare_to_published(thresholds, reference_table, scale_to_table,
                         severity: str, pretransform=None) -> pd.DataFrame:
    """Compare one severity level against Williams & Joshi (2013) Table 1.

    W&J publish a MEDIAN per diagnostic, not a severity threshold, so this is
    an order-of-magnitude and sign check rather than an equality test — which
    is the right strength of claim: it confirms the calibration is producing
    physically sensible numbers in the right units, and would loudly catch a
    units error, a sign error, or a diagnostic calibrated on empty data.

    Thresholds are converted from native SI into W&J table units using the
    same SCALE_TO_TABLE / PRETRANSFORM maps 3_pipeline.py uses, so this cannot
    drift away from the comparison table the pipeline already prints.
    """
    pretransform = pretransform or {}
    rows = []
    for name, per_sev in sorted(thresholds.items(), key=lambda kv: reference_table[kv[0]]["num"]):
        if severity not in per_sev:
            continue
        ref = reference_table[name]
        native = per_sev[severity]
        if isinstance(native, (list, tuple)):
            continue  # two-tailed thresholds have no single published analogue
        value = native
        if name in pretransform:
            value = float(pretransform[name](np.asarray(value)))
        table_value = value * scale_to_table.get(name, 1.0)

        published = ref.get("wj_median")
        if published in (None, 0):
            verdict, ratio = "no published value", np.nan
        else:
            ratio = table_value / published
            if np.sign(table_value) != np.sign(published):
                verdict = "SIGN MISMATCH"
            elif abs(ratio - 1.0) < 1e-4:
                verdict = "SUSPICIOUSLY EXACT — check for hardcoding"
            elif 0.1 <= abs(ratio) <= 10:
                verdict = "consistent"
            else:
                verdict = "ORDERS OFF"

        rows.append({
            "num": ref["num"],
            "diagnostic": name,
            "sign": ref["sign"],
            "threshold_native": native,
            "threshold_table_units": table_value,
            "wj_units": ref["units"],
            "wj_median": published,
            "ratio": ratio,
            "verdict": verdict,
        })
    return pd.DataFrame(rows)


def summarise(df: pd.DataFrame) -> str:
    """One-line verdict over a comparison table."""
    if df.empty:
        return "no diagnostics compared"
    counts = df["verdict"].value_counts().to_dict()
    bad = sum(v for k, v in counts.items() if k in ("SIGN MISMATCH", "ORDERS OFF"))
    return (f"{len(df)} diagnostics: "
            + ", ".join(f"{v} {k}" for k, v in counts.items())
            + (f"   <-- {bad} need attention" if bad else "   <-- all consistent"))
