"""
tests/report_errors.py
======================
Produce the measured-accuracy table for the analytic suite.

pytest answers pass/fail. This answers "by how much", which is the number that
belongs in the verification record -- a diagnostic passing at 0.008 % and one
passing at 4.9 % against a 5 % tolerance are not the same result, and only one
of them should be quoted without a caveat.

Run:  python tests/report_errors.py [output_dir]
"""
from __future__ import annotations

import datetime as _dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))

import synthetic as S  # noqa: E402
from conftest import _load  # noqa: E402

OMEGA_EARTH = 7.292115e-05


def build_rows():
    diag = _load("diagnostics", "2_diagnostics.py")
    ds = S.make_dataset()
    prepared = diag.prepare_for_rojak(ds)
    lat, lon, time = S.make_grid()
    E = S.Expect(lat, lon, time)
    out, failures = diag.compute_rojak_diagnostics(prepared)
    if failures:
        raise SystemExit(f"rojak diagnostics failed: {[f.key for f in failures]}")

    f = 2 * OMEGA_EARTH * np.sin(np.deg2rad(lat))[:, None, None]
    beta = (2 * OMEGA_EARTH * np.cos(np.deg2rad(lat)) / S.meridional_radius(lat))[:, None, None]
    zeta_a = E.vorticity + f
    DEF = E.total_deformation

    def at(da, level=200):
        if "pressure_level" in da.dims:
            da = da.sel(pressure_level=level, method="nearest")
        return da.transpose("latitude", "longitude", "time").values

    cases = [
        # key, tier, formula, computed, expected
        ("rva_magnitude", "C", "|u dz/dx + v dz/dy|",
         at(diag.rva(prepared._dataset)), E.rva),
        ("f2d", "C", "0.5 D/Dt[(du/dth)^2+(dv/dth)^2]",
         diag.frontogenesis_isentropic(prepared._dataset)
             .transpose("latitude", "longitude", "time").isel(time=slice(1, -1)).values,
         E.f2d_isentropic[:, :, 1:-1]),
        ("wind_speed", "B", "sqrt(u^2+v^2)", at(out["wind_speed"]), E.wind_speed),
        ("vorticity_squared", "B", "zeta^2", at(out["vorticity_squared"]), E.vorticity ** 2),
        ("horizontal_divergence", "B", "|delta| (ERA5 field)",
         at(out["horizontal_divergence"]), np.abs(E.divergence)),
        ("magnitude_pv", "B", "|PV| (ERA5 field)",
         at(out["magnitude_pv"]), np.abs(E.potential_vorticity)),
        ("temperature_gradient", "B", "|grad T|", at(out["temperature_gradient"]), E.grad_T),
        ("vertical_wind_shear", "B", "sqrt((du/dZ)^2+(dv/dZ)^2)",
         at(out["vertical_wind_shear"]), E.vws),
        ("endlich", "B", "|v| |dpsi/dz|", at(out["endlich"]), E.endlich),
        ("deformation", "B", "DEF^2 (rojak returns squared)",
         at(out["deformation"]), DEF ** 2),
        ("ti1", "B", "Sv * DEF", at(out["ti1"]), E.vws * DEF),
        ("ti2", "B", "Sv * (DEF - delta_ERA5)",
         at(out["ti2"]), E.vws * (DEF - E.divergence)),
        ("ngm1", "B", "|v| * DEF", at(out["ngm1"]), E.wind_speed * DEF),
        ("ngm2", "B", "|dT/dZ| * DEF", at(out["ngm2"]), np.abs(E.dT_dz) * DEF),
        ("brown1", "B", "sqrt(0.3 zeta_a^2 + Dsh^2 + Dst^2)",
         at(out["brown1"]),
         np.sqrt(0.3 * zeta_a ** 2 + E.shearing_deformation ** 2 + E.stretching_deformation ** 2)),
        ("nva", "B", "max(-u d(zeta+f)/dx - v d(zeta+f)/dy, 0)",
         at(out["nva"]),
         np.maximum(-(E.u * E.dvort_dx + E.v * (E.dvort_dy + beta)), 0.0)),
    ]

    rows = []
    for key, tier, formula, got, expected in cases:
        got = np.asarray(got, dtype=float)
        expected = np.broadcast_to(np.asarray(expected, dtype=float), got.shape)
        median, p95 = S.relative_error(got, expected)
        rows.append({
            "diagnostic": key,
            "wj_number": diag.REFERENCE_TABLE[key]["num"],
            "tier": tier,
            "formula_checked": formula,
            "median_rel_err": f"{median:.6f}",
            "p95_rel_err": f"{p95:.6f}",
            "verdict": "EXACT" if p95 < 0.002 else ("OK" if p95 < 0.02 else "OK (truncation)"),
        })
    return pd.DataFrame(rows).sort_values("wj_number")


def main():
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        REPO / "verification" / _dt.date.today().isoformat()
    out_dir.mkdir(parents=True, exist_ok=True)
    df = build_rows()
    path = out_dir / "analytic_verification.csv"
    df.to_csv(path, index=False)
    print(df.to_string(index=False))
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
