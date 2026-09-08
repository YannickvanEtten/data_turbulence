# PROMPT — full literature audit of the 21 CAT diagnostics

*Paste everything below the line into a fresh session that has the
`data_turbulence` repo and the `Turbulence project/Articles/` PDFs available.
It is written to stand alone: it assumes no memory of previous sessions.*

---

You are auditing the 21 clear-air-turbulence diagnostics implemented in this
repository against their primary literature sources. The question is narrow and
total: **for each diagnostic, is the quantity this code computes the quantity
the paper defines?**

## What you are auditing

Implementations live in `2_diagnostics.py` (hand-written) and in the vendored
`rojak` package (imported). Part of the job is establishing **which code path
actually executes** for each of the 21 — the repository contains both, and a
diagnostic may be hand-written, rojak-sourced, or hand-written *wrapping* a
rojak primitive. Audit the code that runs, not the code that looks canonical.

Primary sources, all in `Turbulence project/Articles/`:

- **Sharman et al. (2006)**, Appendix A — the definitive list, equation numbers
  A1–A40ish. This is the reference of record for most of the 21.
- **Ellrod & Knapp (1992)** — TI1 and TI2.
- **Williams (2017)** — Table 1 (the severity percentile ladder), Table 2
  (thresholds, units, and the canonical 21-item list in the canonical order).
- **Williams & Joshi (2013)** — Table 1 medians.
- **Kaplan (2005)**, **Koch & Caracena (2002)** — NCSU1 and frontogenesis
  lineage.
- **Lee et al. (2023)** — an independent ERA5 implementation of a subset.
- **Prosser et al. (2023)** — the study being replicated. His Figure 4 panel
  order is exactly Williams (2017) Table 2's row order, all 21.

**Brown (1973), *Meteorological Magazine* 102, 347–360 is NOT on disk** and is
the sole authority for the 0.3 coefficient in the Brown index. Record what
rests on it; do not invent its contents.

## Rules of evidence — read these twice

1. **Nothing in this repository counts as evidence about the literature.** Not
   docstrings, not comments, not variable names, not `FORMULA_AUDIT.md`, not
   `STATUS.md`, not any previous audit. Those are *claims to be re-tested*. Two
   sign entries in this project's own reference table have already been found
   wrong after being recorded as checked, so "it was verified before" is not a
   reason to skip a derivation.
2. **Cite equation number and page for every verdict.** "Matches Sharman" is
   not a finding; "matches Sharman (2006) A17, p. 269, including the factor of
   ½" is.
3. **Where the paper is ambiguous, say so and stop.** Record the ambiguity, the
   candidate readings, and which reading the code implements. Do not resolve an
   ambiguity by preferring whatever the code already does.
4. **A magnitude comparison cannot settle a formula question.** This project's
   thresholds are percentiles of its own data, so any constant factor — units,
   resolution, a uniform bias — cancels exactly out of the exceedance field. A
   diagnostic can be off by 8000× against a published table and still be
   correct, and one that matches on magnitude can still have the wrong shape.
   Judge formulas by derivation, not by how big the numbers are.
5. **Report the finding, not a verdict on the project.** If something is
   defensible-but-different, say defensible-but-different.

## Phase 0 — the dependency graph, before anything else

Errors here propagate, and that is the point of this audit. Build the graph
first and write it down:

Identify every **shared building block** — the quantities computed once and
consumed by several diagnostics. At minimum: vertical wind shear, the
Brunt–Väisälä frequency N², the Richardson number, the two deformation
components (shearing and stretching), relative and absolute vorticity, the
horizontal temperature gradient, potential temperature and its gradient,
divergence, the vertical wind-direction derivative, and the machinery that
converts a pressure-level derivative into an altitude derivative.

For each block, list every one of the 21 that consumes it. Then state the
**blast radius**: if this block is wrong, which diagnostics move, and in which
direction.

The Richardson number is the motivating example — it is consumed by more than
one diagnostic — but do not assume the graph; derive it from the code.

## Phase 1 — audit the shared building blocks

Do these before the 21, because a defect here explains several downstream
disagreements at once and a defect downstream may be inherited rather than
local.

For each block check, with citations:

- the defining equation and every constant in it;
- **the metric terms of every horizontal derivative.** On a latitude–longitude
  grid, ∂/∂x requires 1/(a cos φ) ∂/∂λ and ∂/∂y requires (1/a) ∂/∂φ. A missing
  or misplaced cos φ produces a **latitude-structured** error, which is
  invisible in a global median and fatal in a regional exceedance frequency.
  Check this for every gradient, curl and divergence in the code;
- whether the quantity is **computed from u, v, T** or taken from **ERA5's own
  archived field** (divergence, vorticity and potential vorticity are all
  archived). Both are defensible; they are not identical, and the choice must
  be deliberate and recorded;
- the vertical discretisation: which pressure levels, what stencil, and whether
  the hydrostatic conversion from pressure to altitude is applied consistently;
- units at every step, carried symbolically to the final result.

There is a documented history here: this project found that the vendored
library computed the Richardson number as N²/S rather than N²/S². Verify from
the current source which definition each consumer actually receives today.

## Phase 2 — audit all 21, blind

**Do not read `data_prosser/per_diagnostic_vs_prosser.csv`, `STATUS.md` §17, or
any list of "suspect" diagnostics until Phase 2 is written and saved.** The
point is that the derivation must not be steered toward finding a defect where
one is expected, nor toward clearing a diagnostic nobody suspects. Audit all 21
with equal care and equal skepticism.

For each diagnostic, in Williams (2017) Table 2 order:

| field | what to record |
|---|---|
| source | paper, equation number, page |
| formula as published | transcribed symbolically |
| formula as implemented | transcribed from the executing code path |
| variant | which variant of the diagnostic this is, where the literature has more than one, and whether the source being replicated used the same one |
| constants | every numerical coefficient, and its authority |
| sign convention | which tail indicates turbulence, and the derivation of that from the source — not from the code |
| clipping / masking | any `max(·,0)`, absolute value, or NaN handling, and whether the source specifies it |
| inputs | which fields, which levels, computed or archived |
| units | published units vs implemented units |
| verdict | MATCHES / DIFFERS / AMBIGUOUS-IN-SOURCE / UNVERIFIABLE, with the citation |

Pay particular attention to the diagnostics that are **products or residuals of
other quantities**, since they inherit every upstream defect and can also
introduce their own; and to any diagnostic whose definition involves a
**material or time derivative**, where the discrete stencil and the sampling
interval are part of the definition in practice even when the paper is silent.

## Phase 3 — reconcile with the empirical evidence

Only now, open `data_prosser/per_diagnostic_vs_prosser.csv` and `STATUS.md`
§17.5. That file compares each diagnostic's exceedance level and fitted trend
against the values printed in Prosser's own Figure 4 panels.

Then answer, explicitly:

1. Which diagnostics did Phase 2 flag that the empirical comparison did **not**?
   Those are formula defects that cancel out of a percentile-calibrated
   exceedance field — real, and lower priority.
2. Which diagnostics did the empirical comparison flag that Phase 2 found
   **clean**? These are the interesting ones: either the disagreement is
   physical (the vertical level substitution) rather than a bug, or the audit
   missed something. Say which, and why.
3. Where both agree, the case is closed and the fix is specified.

**A concept that should sharpen Phase 3.** Because thresholds are *global*
percentiles and exceedance is measured in a *regional* box, a disagreement in
exceedance *level* is a statement about that diagnostic's
region-versus-globe contrast — that is, about its **spatial structure**, not its
magnitude. Latitude-dependent terms are therefore the prime suspects: the
Coriolis parameter inside absolute vorticity, tan φ factors, and the cos φ
metric terms from Phase 1. Use this to direct Phase 3, and say for each
disagreement whether a latitude-structured error could produce it.

## Phase 4 — propagation and priority

Combine Phase 0's blast radii with Phase 2's verdicts:

- For each defect found, list every diagnostic it reaches and estimate the sign
  and rough size of the effect on each.
- Distinguish **local defects** (wrong in one diagnostic) from **inherited
  defects** (one shared block, several symptoms). An audit that fixes five
  symptoms separately when there is one cause has failed.
- Rank the fixes by (diagnostics affected) × (confidence in the finding), and
  state for each what evidence would confirm it after the fix — which number
  should move, in which direction, and by roughly how much.

## Specific questions that must be answered

These are open in the project record. Answer each with a citation or an
explicit "cannot be resolved from available sources":

1. The coefficient in the Brown index, whose only authority is a paper not on
   disk. What is it, what does it multiply, and what currently justifies it?
2. Whether the leading sign of the frontogenesis function as implemented
   matches the published form, and which frontogenesis variant the replicated
   study used.
3. Whether the one-sided clipping applied to the vorticity-advection and NCSU
   diagnostics is specified by the source or is an implementation choice.
4. Whether the Colson–Panofsky index's scale factor and units are correct,
   given it is published in non-SI units.
5. For all 21: whether the tail that indicates turbulence is the upper one,
   derived from the source rather than assumed.

## Output

Write `AUDIT_diagnostics_vs_literature.md` in the repository root:

1. **Executive summary** — defects found, ranked by blast radius, each in one
   sentence.
2. **The dependency graph** (Phase 0), as a table or diagram.
3. **The 21 audits** (Phase 2), one section each, in Williams Table 2 order.
4. **Reconciliation** (Phase 3).
5. **Propagation and priority** (Phase 4), with the confirming evidence for
   each proposed fix.
6. **Unresolved** — every ambiguity, every source not on disk, every question
   the literature does not settle.

Then propose — do not apply — the code changes, one per defect, each with its
citation. A fix applied before the audit is reviewed is a fix nobody checked.

## What not to do

Do not run jobs, download data, or recompute diagnostics; this is a reading
task. Do not fix code as you go. Do not treat a previous audit's conclusion as
a shortcut. Do not report a diagnostic as verified when the source is
ambiguous — "the paper does not specify, and we chose X" is a legitimate and
useful finding, and pretending otherwise is how the two known-wrong sign
entries survived their first audit.
