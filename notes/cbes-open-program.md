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
