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
