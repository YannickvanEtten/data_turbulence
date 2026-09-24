"""
Prosser et al. (2023) Supporting Information, Figures S4-LOG and S4-SOG,
transcribed 2026-09-16 for the final per-diagnostic scorecard.
(FINAL-CHECKS-2026-09-16)

WHY THIS FILE EXISTS
--------------------
`ada/prosser_scorecard.py` lists S4-SOG as "figure only, no stated counts".
That is true of the TEXT (Prosser states 17/21 only for MOG), but the figure
itself prints `rel=`, `abs=` and `p=` inside every one of its 21 panels,
exactly as main-text Figure 4 does. So S4-LOG and S4-SOG are two more exact
per-diagnostic targets -- 42 printed numbers that had not been used.

SOURCE
------
`Turbulence project/Articles/Supporting information si-s01 (1).docx`,
word/media/image8.emf (S4-LOG) and image9.emf (S4-SOG), in document order
rId14 and rId15, rendered to PNG with LibreOffice (`soffice --convert-to png`).
Panel letters and layout are identical to main-text Figure 4.

PRECISION -- same rules as ada/prosser_published.py
---------------------------------------------------
  rel, p        EXACT, printed in each panel.
  level_1979    APPROXIMATE, read off the red 1979 cross. Validated in
                aggregate: the mean of the 21 read-off levels is 5.338 % at LOG
                and 0.2012 % at SOG, against Table 1's exact 466.5 h / 8760 h
                = 5.325 % and 17.7 h / 8760 h = 0.2021 % -- within 0.3 % and
                0.5 %. Individual panels are good to roughly one fifth of a
                gridline; use them for ratios, never for a precise claim.

name -> (panel, rel %, p, approx 1979 level %)
"""
from __future__ import annotations

S4_LOG = {
    "negative_richardson":   ("a", 15.7, 0.3,    0.37),
    "vertical_wind_shear":   ("b", 15.8, 4e-05,  8.5),
    "colson_panofsky":       ("c", 32.6, 0.006,  0.81),
    "f2d":                   ("d", 13.5, 0.2,    0.57),
    "brown1":                ("e", 18.5, 2e-05,  6.6),
    "brown2":                ("f", 16.7, 1e-05,  9.4),
    "ti1":                   ("g", 16.6, 4e-05,  8.2),
    "ti2":                   ("h", 12.7, 0.0002, 7.55),
    "deformation":           ("i", 17.5, 0.0001, 5.75),
    "magnitude_pv":          ("j", 22.3, 0.001,  8.1),
    "vorticity_squared":     ("k", 19.3, 0.0003, 6.35),
    "temperature_gradient":  ("l", 6.6,  0.02,   10.7),
    "wind_speed":            ("m", 31.9, 0.03,   2.95),
    "endlich":               ("n", 14.1, 0.0002, 5.25),
    "ngm1":                  ("o", 26.9, 0.0005, 5.6),
    "ngm2":                  ("p", 23.9, 0.002,  1.9),
    "ubf":                   ("q", 18.4, 9e-05,  6.3),
    "horizontal_divergence": ("r", 3.3,  0.4,    4.1),
    "ncsu1":                 ("s", 20.7, 0.0009, 2.75),
    "nva":                   ("t", 16.1, 0.002,  4.95),
    "rva_magnitude":         ("u", 19.0, 0.002,  5.4),
}

S4_SOG = {
    "negative_richardson":   ("a", -13.9, 0.3,    0.0125),
    "vertical_wind_shear":   ("b", 54.9,  6e-07,  0.47),
    "colson_panofsky":       ("c", 55.1,  0.002,  0.13),
    "f2d":                   ("d", 0.1,   1.0,    0.013),
    "brown1":                ("e", 71.2,  1e-07,  0.145),
    "brown2":                ("f", 68.9,  3e-08,  0.42),
    "ti1":                   ("g", 58.5,  6e-07,  0.37),
    "ti2":                   ("h", 31.8,  0.0004, 0.37),
    "deformation":           ("i", 53.5,  5e-06,  0.145),
    "magnitude_pv":          ("j", 90.1,  2e-09,  0.45),
    "vorticity_squared":     ("k", 69.4,  5e-07,  0.08),
    "temperature_gradient":  ("l", 39.5,  6e-06,  0.40),
    "wind_speed":            ("m", 102.2, 0.03,   0.05),
    "endlich":               ("n", 36.3,  3e-05,  0.23),
    "ngm1":                  ("o", 77.6,  0.0002, 0.21),
    "ngm2":                  ("p", 42.0,  0.001,  0.051),
    "ubf":                   ("q", 43.7,  9e-05,  0.135),
    "horizontal_divergence": ("r", 8.1,   0.2,    0.117),
    "ncsu1":                 ("s", 9.7,   0.5,    0.036),
    "nva":                   ("t", 35.8,  0.002,  0.19),
    "rva_magnitude":         ("u", 34.8,  0.004,  0.20),
}

# Prosser Table 1, Annual row, 1979 hours -> % of the year. Exact.
TABLE1_ANNUAL_1979_PCT = {"light": 466.5 / 8760 * 100,
                          "moderate": 70.0 / 8760 * 100,
                          "severe": 17.7 / 8760 * 100}

assert len(S4_LOG) == 21 and len(S4_SOG) == 21
assert set(S4_LOG) == set(S4_SOG)

# The `abs=` line printed in each panel: the slope of his regression, in
# percentage points per year. Transcribed from the same two renders.
#
# SELF-CHECK, run at import: for a straight line through 42 years,
#     abs ~= rel/100 * level_1979 / 41 .
# `rel` is exact and `level_1979` is a read-off, so agreement to a few per cent
# tests BOTH transcriptions at once — a mis-read level or a mistyped abs shows
# up here rather than in a figure caption three weeks later.
S4_LOG_ABS = {
    "negative_richardson": 0.001, "vertical_wind_shear": 0.033,
    "colson_panofsky": 0.006, "f2d": 0.002, "brown1": 0.030, "brown2": 0.038,
    "ti1": 0.033, "ti2": 0.023, "deformation": 0.024, "magnitude_pv": 0.044,
    "vorticity_squared": 0.030, "temperature_gradient": 0.017,
    "wind_speed": 0.023, "endlich": 0.018, "ngm1": 0.037, "ngm2": 0.011,
    "ubf": 0.028, "horizontal_divergence": 0.003, "ncsu1": 0.014,
    "nva": 0.019, "rva_magnitude": 0.025,
}
S4_SOG_ABS = {
    "negative_richardson": -0.000, "vertical_wind_shear": 0.006,
    "colson_panofsky": 0.002, "f2d": 0.000, "brown1": 0.003, "brown2": 0.007,
    "ti1": 0.005, "ti2": 0.003, "deformation": 0.002, "magnitude_pv": 0.010,
    "vorticity_squared": 0.001, "temperature_gradient": 0.004,
    "wind_speed": 0.001, "endlich": 0.002, "ngm1": 0.004, "ngm2": 0.001,
    "ubf": 0.001, "horizontal_divergence": 0.000, "ncsu1": 0.000,
    "nva": 0.002, "rva_magnitude": 0.002,
}


def _check_abs(src, abs_map, label):
    bad = []
    for k, v in src.items():
        implied = v[1] / 100.0 * v[3] / 41.0          # rel% * level% / 41 yr
        got = abs_map[k]
        if abs(implied - got) > max(0.0006, 0.15 * abs(got)):
            bad.append(f"{label} {k}: printed abs={got:+.3f}, "
                       f"rel and level imply {implied:+.4f}")
    return bad


_PROBLEMS = _check_abs(S4_LOG, S4_LOG_ABS, "S4-LOG") + \
            _check_abs(S4_SOG, S4_SOG_ABS, "S4-SOG")
if _PROBLEMS:                       # never silent: a bad transcription is a bug
    raise AssertionError("S4 transcription is not self-consistent:\n  "
                         + "\n  ".join(_PROBLEMS))
