"""
Published values from Prosser et al. (2023), transcribed for comparison.

WHERE EACH NUMBER COMES FROM, because the precision differs and pretending
otherwise is how a replication fools itself:

  rel, abs, p     EXACT. Printed inside each panel of his Figure 4 and read
                  from a 400-dpi render of page 6. `rel` is the relative change
                  in the fit from 1979 to 2020 (%), `abs` is the slope of the
                  regression line (%/yr), `p` is the Wald-test p value for the
                  slope. His caption states this.
  level_1979      APPROXIMATE, ±0.02 in the units given. Read off the position
                  of the red 1979 cross against the panel's y axis. Use it for
                  ratios and orders of magnitude, never for a precise claim.
  panel           his subfigure letter, so a disagreement can be checked
                  against the paper by eye in one step.

Table 1 (the seasonal/severity grid) is tabulated in the paper and is exact;
it lives in ada/prosser_scorecard.py, not here.

Prosser's Figure 4 is ANNUAL, MOG, North Atlantic box (36-60N, 55-10W), ERA5
at 197 hPa, n=42. Compare only against our annual MOG run.
"""

from __future__ import annotations

# name-as-we-call-it -> (panel, his title, rel %, abs %/yr, p, approx 1979 level %)
FIGURE4 = {
    "negative_richardson":   ("a", "Negative Richardson number",              -7.5, -0.000, 0.6,     0.039),
    "vertical_wind_shear":   ("b", "Magnitude of vertical shear of horizontal wind", 38.4, 0.016, 8e-07, 1.68),
    "colson_panofsky":       ("c", "Colson-Panofsky index",                   45.4,  0.003, 0.002,   0.285),
    "f2d":                   ("d", "Frontogenesis function",                   0.3,  0.000, 1.0,     0.055),
    "brown1":                ("e", "Brown index",                             45.4,  0.008, 2e-07,   0.70),
    "brown2":                ("f", "Brown energy dissipation rate",           46.2,  0.018, 5e-08,   1.60),
    "ti1":                   ("g", "Variant 1 of Ellrod's turbulence index",  42.3,  0.014, 3e-07,   1.40),
    "ti2":                   ("h", "Variant 2 of Ellrod's turbulence index",  27.0,  0.009, 3e-05,   1.37),
    "deformation":           ("i", "Flow deformation",                        37.5,  0.006, 6e-06,   0.70),
    "magnitude_pv":          ("j", "Magnitude of potential vorticity",        50.6,  0.020, 4e-07,   1.62),
    "vorticity_squared":     ("k", "Relative vorticity squared",              37.1,  0.004, 2e-05,   0.47),
    "temperature_gradient":  ("l", "Magnitude of horizontal temperature gradient", 26.0, 0.011, 2e-05, 1.75),
    "wind_speed":            ("m", "Wind speed",                              75.6,  0.005, 0.02,    0.28),
    "endlich":               ("n", "Wind speed times directional shear",      27.7,  0.006, 2e-05,   0.92),
    "ngm1":                  ("o", "Flow deformation times wind speed",       50.9,  0.010, 0.0003,  0.83),
    "ngm2":                  ("p", "Flow deformation times vertical temperature gradient", 33.8, 0.002, 0.001, 0.24),
    "ubf":                   ("q", "Magnitude of residual of nonlinear balance equation", 31.6, 0.006, 6e-05, 0.72),
    "horizontal_divergence": ("r", "Magnitude of horizontal divergence",       6.7,  0.001, 0.2,     0.54),
    "ncsu1":                 ("s", "Version 1 of North Carolina State University Index", 4.3, 0.000, 0.7, 0.070),
    "nva":                   ("t", "Negative absolute vorticity advection",   28.4,  0.005, 0.001,   0.72),
    "rva_magnitude":         ("u", "Magnitude of relative vorticity advection", 29.0, 0.006, 0.003,  0.83),
}

# Does the diagnostic use the 175/225 hPa vertical stencil, or is it computed
# entirely on the 200 hPa surface? PLAN_full_year_calibration.md §4.
#
# This split is the whole point of the per-diagnostic comparison. The project's
# standing explanation for its 12-20% level deficit is the vertical stencil
# (ours 50 hPa, Prosser's 18 hPa). If that is right, the deficit must
# concentrate in the ten stencil diagnostics and roughly vanish in the eleven
# single-level ones. If it is spread evenly, the explanation is wrong.
USES_VERTICAL_STENCIL = {
    "vertical_wind_shear": True, "negative_richardson": True,
    "colson_panofsky": True, "ti1": True, "ti2": True, "ngm2": True,
    "endlich": True, "ncsu1": True, "brown2": True, "f2d": True,
    "wind_speed": False, "deformation": False, "horizontal_divergence": False,
    "vorticity_squared": False, "magnitude_pv": False,
    "temperature_gradient": False, "nva": False, "rva_magnitude": False,
    "ubf": False, "brown1": False, "ngm1": False,
}

# The counts Prosser states in his main text, which are exact and are the
# strongest single check on Figure 4.
FIGURE4_SUMMARY = dict(n_significant=17, n_total=21, max_change=75.6,
                       n_significant_downward=0)

assert len(FIGURE4) == 21
assert len(USES_VERTICAL_STENCIL) == 21
assert sum(USES_VERTICAL_STENCIL.values()) == 10


# Panel y-axis limits, read off the same render. APPROXIMATE but useful: they
# let our panels be drawn on his scale so the two figures can be laid side by
# side and compared by eye, which is the point of matching axes at all.
FIGURE4_YLIM = {
    "negative_richardson": (0.01, 0.06), "vertical_wind_shear": (1.4, 3.0),
    "colson_panofsky": (0.2, 0.6), "f2d": (0.025, 0.085),
    "brown1": (0.68, 1.20), "brown2": (1.4, 3.0), "ti1": (1.2, 2.4),
    "ti2": (1.2, 2.0), "deformation": (0.6, 1.1), "magnitude_pv": (1.0, 3.0),
    "vorticity_squared": (0.4, 0.8), "temperature_gradient": (1.3, 2.5),
    "wind_speed": (0.0, 0.8), "endlich": (0.7, 1.5), "ngm1": (0.5, 1.7),
    "ngm2": (0.10, 0.38), "ubf": (0.6, 1.1),
    "horizontal_divergence": (0.40, 0.68), "ncsu1": (0.03, 0.11),
    "nva": (0.60, 1.20), "rva_magnitude": (0.8, 1.4),
}
assert len(FIGURE4_YLIM) == 21
