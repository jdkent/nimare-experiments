# Coverage of the combined-estimator design document

Assessment of `fd78b2ef-nimare_combined_estimator_research.md` (429 lines, 16 September 2026)
after working through its build order. Branch `claude/marginal-effect-combined-estimator` of
`neurostuff/NiMARE`; algebra in `proofs/`, measurements in `experiments/`.

**Totals.** 4 new modules, 56 tests, 21 proof files carrying 168 verified sympy claims, 4
calibration beds reproducing the document's own numerical experiments.

---

## Part 1 -- what was built, in the document's own order

| Document item | Status | Where |
|---|---|---|
| §12 step 1, correct CBES semantics | done earlier | PR #1146, 9 commits |
| §12 step 3a, scalar exact-bounds reference | **done** | `censored.py`, 18 tests |
| §12 step 3b, image-corrected predictor | **done** | `marginal.py` control variate |
| §6.2 exact interval contribution | **done** | `censored_loglik` |
| §6.2 reporting operator, retention | **done** | `retention=`, `retention_roles` |
| §12 observation record (9 states) | **done** | `ObservationState` |
| §2 practical prevalence π_δ⁺ | **done** | `practical_prevalence` |
| §7 propensity-augmented mean | **done** | `augmented_mean` |
| §6.3 step 2, small spatial blocks | **done, with a measured limit** | `blocks.py`, 13 tests |
| §6.3 step 3, exact per-study likelihood | **done** | `study_index=` |
| §10 power, assurance, ceiling | **done** | `planning.py`, 12 tests |
| §8.1 calibration (9 cells) | **reproduced** | `censored_reference_calibration.py` |
| §8.2 calibration (2 of 3 columns) | **reproduced** | `block_reference_calibration.py` |

Four of the document's own claims were sharpened rather than merely implemented:

* §10 says assurance "in general is not Power(n, E θ)" without a direction. The direction is
  derivable and **flips at the 50%-power effect**, so plug-in power is optimistic exactly where
  a study is being designed to be adequately powered. Assurance also has a ceiling, Φ(|m|/s),
  for directional rejection.
* §8.2 finds heights add little beyond a reporting indicator. That is **reproducible from the
  Fisher informations alone** (analytic 0.929 against its measured 0.943/0.909/0.917) *and* it
  is a property of a strict threshold: at a liberal cut the height carries essentially all the
  information and the interval is 27× narrower.
* §4.1's partition objection is now quantified: using `1 − p_report` for silence overstates it
  by exactly the dropped event's own probability.
* §6.1's caution that one study "does not become thousands" is a **theorem** with a number:
  per-study information about (m, τ²) is capped at diag(1/τ², 1/(2τ⁴)) whatever the voxel
  count, so K = 2/r² studies with images are needed for relative precision r on τ² — 200 for a
  tenth.

---

## Part 2 -- assumptions that could not be met, and why

### A. The paired calibration corpus (§6.4) — **not met, and it is the binding one**

The document calls this decisive and places §6.1, §6.2's calibration extension, §6.3 step 3's
real reporting rules, §12 step 5's locked test sets and §11's independent empirical tests
downstream of it.

Evidence: `experiments/calibration_corpus_feasibility.py`.

* The cached paired collection satisfies **0 of its 6** verification requirements (population,
  contrast direction, analysis type, sample, mask, statistic convention). Its coordinate table
  carries `['id','study_id','contrast_id','x','y','z','space']` and nothing else.
* Both sources it names are **reachable**, so access is not the blocker.
* Five of six requirements exist as fields across the two APIs, but the ones a censored
  likelihood needs are the sparsest: on 83 sampled unthresholded maps, degrees of freedom
  (`statistic_parameters`) is filled in **1%** of the time and contrast identity **18%**,
  against 95% for map type and 88% for subject count. `analysis_level` is 51%, and the
  one-sample/two-sample design is **not a field at all** — which §2 is explicit cannot be
  inferred from a total count.
* **Cohort identity has no field in either API.** §6.4 requires folds grouped by cohort, since
  shared cohorts across articles would otherwise leak between training and test. It cannot be
  derived; it must be audited by hand or obtained from authors.
* NeuroStore's documented `has_*` filters do not filter — every combination returns the full
  84,472 — so no matched-pair count is reported here rather than a guessed one.

**Why unmet:** this is a data-curation project, not a modelling one. Nothing in the mathematics
blocks it.

### B. The moderate-rank spatial model's between-study covariance (§6.1) — **provably not met at the targeted collection sizes**

Evidence: `proofs/information_per_study_is_bounded.py` (8 claims).

Writing a study's own level as η = m + U, every observation it produces depends on the
parameters only through η, so η is **sufficient** and the study's data are a noisy measurement
of *one draw* from N(m, τ²). Its information about the between-study parameters is therefore
capped at that draw's own, **whatever the voxel count**. Numerically confirmed over every block
count the computation can resolve: 1.2%, 10%, 29%, 58% of the ceiling at 1, 2, 5, 20 blocks.

Consequence: the document's own regimes (1 image, or 8) can give **at most 8 draws**, so τ² has
a relative standard error no better than √(2/8) = 50%. The spatial coefficients are learnable
from one study's voxels; the between-study covariance is not. It has to be borrowed, with
propagated uncertainty — which is what the document says, now with the size attached.

### C. The reporting-protocol hierarchy (§6.4) — **not learnable from tables, measured**

Evidence: `proofs/what_identifies_retention.py` (7 claims) plus its measurements.

A study's score direction in (m, τ², ρ) is fixed by its threshold and precision alone, so
studies sharing both give parallel gradients and rank-one information; three parameters need
three non-collinear (threshold, precision) pairs. That is local identification only. Measured
from the observed information and confirmed by simulation to ~15%:

| configuration | se(ρ) |
|---|---|
| 100 tables, common threshold | 1.83 |
| 100 tables, spread thresholds | 0.53 |
| 8 images + 100 tables | 0.19 |
| 20 images + 2000 tables | 0.07 |

ρ lives in [0,1]. So a free per-study detection curve is not estimable from coordinate tables at
any realistic corpus size — which is what §6.4 warns ("too flexible to learn from one table and
can absorb the effect itself"), quantified. Retention must come from images or external
calibration. Misspecifying it is not benign: §8.1's coverage collapses to **0%** at 500 tables.

### D. Whole-brain composite likelihood (§6.3 step 3) — **partially met; the correction is not built**

Evidence: `proofs/composite_block_likelihood.py` (6 claims).

Treating a study's blocks as independent understates the likelihood and **overstates the
information**, because it counts one draw of the study effect as many. Exact-to-composite
standard-error ratio: 1.000, 0.798, 0.615, 0.421, **0.204** at 1, 2, 4, 10, 50 blocks per study
— a whole-brain composite fit would report an interval five times too narrow.

The exact per-study integral is now implemented (`study_index=`), so the defect is avoidable at
small scale. What is **not** built is the sandwich or Godambe correction a genuinely whole-brain
composite fit would need; `standard_error_inflation` measures the gap but does not close it.

### E. General spatial blocks (§6.3 step 2) — **met only in the strict-threshold regime**

Evidence: `experiments/exchangeable_block_error.py`, which states a 0.01 kill condition in
advance and **fails it at 0.030**.

The block model conditions on the study effect and leaves elements independent, so every pair
correlates equally; real correlation decays. Against exact multivariate-normal orthant
probabilities at matched mean correlation, the exchangeable model **understates silence** —
overstating how often a block reports — by up to 0.030 at nine elements with slow decay and a
liberal cut. Restricted to strict cuts the worst error is 0.005; at mean within-block
correlation ≤ 0.15 it is 0.010, marginally *outside* the bound.

So it is usable for the strict thresholds published tables come from, and it is not a general
substitute for a spatial model. The conditional-multivariate-normal route the document names is
not implemented. The module's own closed form matches the exchangeable orthant probability to
1e-5, so what fails is the assumption, not its implementation.

### F. §11's comparator against the real SDM-PSI — **not met, environment**

Verified rather than assumed: no `R`, no `Rscript`, no `matlab`, no `octave`; neither
`metansue` (CRAN) nor `sdm-psi` is pip-installable. §11 says "do not use a hand-written
approximation and label it SDM-PSI", so no comparison is offered rather than a mislabelled one.

### G. §8.2's naive arm — **unreproducible, construction under-specified**

Four candidate constructions give asymptotic biases of +0.060, +0.130, +0.490 and +0.000
against the +0.256 its RMSE implies; the middle two bracket it. The two *modelled* columns
reproduce, and the correct indicator-only arm comes out unbiased, which confirms the generative
model was read correctly. Only the deliberately-wrong foil cannot be reconstructed.

### H. §12's class API, and §10's cohort validation — **scoped out, deliberately**

`ReportingCalibration.fit(...)` / `MixedEffectSize(...)` are built as functions rather than
those classes. A `ReportingCalibration` object with nothing to calibrate against would be a
shell, and the document's own step 2 puts the corpus before the model. §10's "validate
prospective sample-size recommendations on new cohorts" is the document's final criterion and
needs new cohorts.

---

## Part 3 -- the one-line summary

Everything in the document that is **mathematics or a scalar/block likelihood** is built,
derived first, and calibrated against its own numbers. Everything that needs a **matched
map-and-table corpus** is blocked on curation, and three of its assumptions turn out to have
hard quantitative limits that it states only as cautions: per-study information is capped
regardless of voxel count, retention is not estimable from tables, and the independent-block
approximation costs a factor of five in interval width.
