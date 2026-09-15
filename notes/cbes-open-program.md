# Open program: estimands and uncertainty

A working list, kept in priority order, with findings folded back in as they arrive. Written so
the state is legible without reading the whole transcript.

## Estimands

**E1. What does the point-process field estimate?** *Answered.* The conditional magnitude, with
about a quarter of attenuation at prevalence 0.25. Magnitude heterogeneity does not shift it at
all; only zero-inflation does, weakly. So the intensity formulation does not change the estimand
-- it targets the same quantity CBES's `g` targets, which no image-based meta-analysis estimates.

**E2. What does `prevalence` converge to, and is the compression invertible?** *Open, next.* It is
documented as compressed toward the middle -- a true 0.25 reads 0.49 to 0.60 -- but never
characterised. Three candidates: the fraction of studies with a non-null effect; the fraction that
*report* near the voxel, which is strictly smaller since a study can have an effect and miss its
threshold; or neither. The observed values exceed both candidates, so the compression is an
artefact rather than an estimand distinction -- but if it is an *affine* function of the truth
with stable coefficients, it is correctable, and that matters because the power argument in E5
needs prevalence on a real scale.

**E3. Does the estimand move with the distribution of sample sizes and thresholds?** *Open.* The
rate function's curvature depends on `u_k` and `sqrt(N_k)`, so a quasi-arithmetic mean generated
by it should shift when the collection's composition shifts, even with the effects held fixed.
If so there is no single estimand across a heterogeneous literature, which is worse than an
estimator bias.

**E4. What does the ALE statistic converge to?** *Open.* Used universally for region selection and
never characterised as an estimand. Candidate: a monotone function of prevalence, which would make
it the right input for the power ceiling in E5 and the wrong one for magnitude.

**E5. Is achievability a better output than magnitude?** *Open, conceptual.* Power under zero
inflation is capped at `pi + (1 - pi) * alpha`, so below prevalence 0.8 no sample size reaches 80%.
That makes prevalence, not magnitude, the decision-relevant quantity for planning -- and an
estimator that reported an achievability bound would be more useful than one reporting an
uncalibrated effect size.

## Uncertainty

**U1. Does the reported standard error cover?** *Partly done, needs redoing.* Coverage of 94.5% to
98.4% was measured against the estimator's *own* censored mixture likelihood, which tests the
arithmetic rather than the model. Against a held-out reference it has never been checked, and the
estimand mismatch in E1 guarantees it will fail there for reasons unrelated to the standard error.

**U2. What should the interval be when the scale is unidentified?** *Open, and conceptually the
most serious.* `g`'s scale is not identified from coordinates, yet `se` is finite and `z` is
reported. A finite interval on an unidentified quantity is a misrepresentation; the honest object
is either an interval on the relative map or an interval that includes the scale's own
uncertainty. `scale_interval_` exists but does not propagate into `se`.

**U3. Interval coverage in the point-process model.** *Open.* Entirely unmeasured. The observed
information is available from the same likelihood, so this is mostly bookkeeping.

**U4. Does the fixed tau-squared under-state uncertainty?** *Open.* The between-study variance is
estimated once about the naive weighted mean and held fixed, and is documented as biased low.
Biased low means intervals too narrow, which would show as under-coverage.

**U5. Uncertainty under estimand mismatch.** *Open.* If an interval for the conditional is scored
against a marginal reference, coverage fails however good the interval is. Worth separating from
U1 explicitly, because the two look identical in a coverage table.

## Reflection on success

**S1.** The criteria in `cbes-success-criteria.md` score point estimates almost exclusively. If a
joint model's benefit is precision, nothing there sees it -- which is what the exchange sweep's
seed spreads suggested.

**S2.** No criterion asks whether the *uncertainty* is honest. A method can pass every accuracy
threshold while reporting intervals that do not cover, and that is arguably the worse failure for
a literature.

**S3.** The achievability framing in E5 suggests a criterion nobody states: does the output let a
reader decide whether their planned study is possible? That is a different and more useful test
than any correlation.

---

## U2 answered: three defects in how uncertainty is reported

Read from the code and confirmed arithmetically.

### U2a. `scale_interval_` is a sample range, and its error changes sign with donor count

It is set to `(min(per_donor), max(per_donor))`. A sample range *grows* with the sample size while
uncertainty about the common scale *shrinks* like one over the square root of it, so the two
diverge in opposite directions. Drawing per-donor estimates from a normal with mean 0.60 and
standard deviation 0.08, twenty thousand times:

| donors | reported (min, max) width | honest 95% interval width | ratio |
| --- | --- | --- | --- |
| **2** | 0.0909 | 0.1781 | **0.51** |
| 3 | 0.1351 | 0.1601 | 0.84 |
| 5 | 0.1866 | 0.1323 | 1.41 |
| 10 | 0.2467 | 0.0967 | 2.55 |
| 20 | 0.2985 | 0.0692 | 4.32 |

At the **two-donor floor the estimator requires for `g_absolute`**, and which the docstring
specifically argues for, the reported interval is **half** the honest width -- it understates
exactly where it matters most. By twenty donors it overstates by 4.3-fold, which also reverses
the incentive to supply more images. The fix is a standard error of the mean across donors rather
than their range.

### U2b. `se` does not include the scale's uncertainty

`g_absolute` is the same array as `g`, so its standard error is `g`'s, which is conditional on the
scale being exactly right. The description text says the magnitudes should be read as an order of
scale, which is honest prose, but the number a reader will use does not reflect it. An interval on
a partially identified quantity that omits the identification uncertainty is a misrepresentation
however carefully the surrounding paragraph is worded.

### U2c. `g_marginal` has no standard error at all

It is emitted as `fit["g"] * fit["prevalence"]` with no variance propagated. So the one magnitude
map with a checkable reference -- the only one sharing an estimand with an image-based
meta-analysis -- is the one carrying no uncertainty. Propagating it needs the covariance of the
prevalence and the magnitude, which the observed information already contains, since the
prevalence is profiled out of it by a Schur complement.

---

## E2 answered: the compression is affine but not invertible, and worse than compressed

`prevalence_calibration` sweeps a grid of true prevalence at three configurations, 12 simulations
per cell, 24 studies, read at the truth voxel.

| N range | effect | true pi | fraction reporting | pi_hat | sd |
| --- | --- | --- | --- | --- | --- |
| 20-40 | 0.8 | 0.00 | 0.035 | **0.218** | 0.213 |
| 20-40 | 0.8 | 0.40 | 0.215 | 0.530 | 0.279 |
| 20-40 | 0.8 | 1.00 | 0.427 | 0.901 | 0.136 |
| 20-40 | 0.4 | 0.00 | 0.035 | 0.218 | 0.213 |
| 20-40 | 0.4 | 0.40 | 0.031 | 0.157 | 0.217 |
| 20-40 | 0.4 | 0.80 | 0.069 | 0.427 | 0.317 |
| 20-40 | 0.4 | 1.00 | 0.080 | **0.302** | 0.216 |
| 10-200 | 0.8 | 0.00 | 0.035 | 0.144 | 0.184 |
| 10-200 | 0.8 | 1.00 | 0.503 | 0.936 | 0.082 |

Affine fits against the two candidate truths:

| N range | effect | slope vs true pi | intercept | r | slope vs reporting | r |
| --- | --- | --- | --- | --- | --- | --- |
| 20-40 | 0.8 | 0.719 | 0.216 | 0.994 | 1.831 | 0.994 |
| 20-40 | 0.4 | **0.176** | 0.173 | **0.684** | 3.619 | 0.789 |
| 10-200 | 0.8 | 0.804 | 0.063 | 0.981 | 1.675 | 0.966 |

Three findings.

**There is a floor of 0.14 to 0.22 at a true prevalence of zero.** The map reads a fifth of
studies having an effect where none do. That is the same phenomenon as the point-process null
floor, on a different scale.

**The relationship is strikingly affine for strong effects but the coefficients are not stable.**
Correlations of 0.994 and 0.981 mean a correction exists in principle; slopes of 0.719 against
0.804 and intercepts of 0.216 against 0.063, differing only in the sample-size range, mean you
would need to know the collection's effect size and N distribution to choose it -- and knowing the
effect size removes the reason to want the map.

**At a weak effect it is non-monotone.** 0.218, 0.194, 0.157, 0.270, 0.427, 0.302 as the truth
goes 0 to 1: it reads *lower* at a prevalence of 1.0 than at 0.8. The per-simulation standard
deviation, 0.22 to 0.32, exceeds the entire range of the means. For weak effects prevalence is
close to noise.

**The caveat on all of it.** This varies prevalence *across collections* at one voxel. The
docstring claims something else -- that ordering across *voxels within one map* survives -- which
is a different quantity and could still hold. `prevalence_within_map` tests the claim as made,
with four sites of known differing prevalence inside a single fit.

This also closes the power argument. Achievability needs prevalence on a real scale to say which
side of the 0.8 ceiling a planned study sits on, and none of the above supports putting it there,
least of all for the weak effects where the question is live.
