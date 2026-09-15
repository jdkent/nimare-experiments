# What counts as success

Written after a session in which several conclusions had to be withdrawn because the *metric* was
wrong rather than the estimator. Each measure below is chosen to be hard to fool in the specific
way its predecessor was fooled, and each comes with a threshold, because a criterion without one
is a description.

---

## The six measures

Ordered by how much of the estimator's claim rests on them.

### 1. Response — does the estimate move with the truth?

**Measure.** Regress the held-out truth on the estimate and report the **slope**. Not the
correlation, not the ratio of means.

**Threshold.** Slope above 0.5. A slope of 1 means a unit of estimated effect corresponds to a
unit of real effect; a slope near 0 means the map is reading something else, whatever its
correlation.

**Why this one first.** It is the single number that would have caught the central failure
immediately. A reference effect spanning elevenfold came back spanning 1.2-fold while the mean
ratio sat near 1 and the correlation looked respectable. The slope was 0.08 to 0.18 the whole
time.

### 2. Calibration — is the number right, where it matters?

**Measure.** The ratio of estimate to truth **within each decile of the truth above the median**,
reported as a profile rather than a summary.

**Threshold.** Within [0.8, 1.25] in every decile above the median.

**Why stratified.** An overall ratio is the average of a large over-estimate in the empty majority
and whatever happens where the effect is. Scoring only the top quartile is the opposite error and
is how an earlier comparison made a compressed map look calibrated.

### 3. Localisation — does it find the effects?

**Measure.** Area under the curve for recovering the truth's top decile, over covered voxels.
Rank correlation as a secondary.

**Threshold.** AUC above 0.9 where two or more images are present.

**Why not the correlation.** Plain `r` over all covered voxels is inflated by spatial smoothness
and dominated by the near-zero majority, and it inflates more for the arm with more studies, which
makes it unfair in exactly the comparisons that matter.

### 4. Error control — familywise, and named as such

**Measure.** The familywise rate, `P(any rejection)`, under global nulls. Never the mean share of
voxels rejected, and the two must never appear in the same column.

**Threshold.** Within two binomial standard errors of nominal, in *every* null tested, and the
nulls must include a **spatially clustered** one and one with sample sizes spanning two orders of
magnitude. Both are conditions under which this estimator has failed and independent noise has
not.

**Paired with power**, at a known focal effect, because controlling the error rate by never
rejecting is not success.

### 5. Convention robustness — would two literatures agree?

**Measure.** Hold the studies fixed, vary only how a paper would have tabulated them -- FDR,
voxelwise family-wise, cluster extent, crossed with reporting a cluster by its maximum or its
centre of mass -- and report the spread of measures 1 to 3 across those six.

**Threshold.** Stratum ratios agreeing within 20%.

**Why it earns a place.** This is the measure the current estimator fails worst: the same studies
give a top-stratum ratio from 0.82 to 2.29 depending only on the convention. Two collections
reporting identical effects under different conventions do not agree, and no amount of internal
validation reveals that, because it is invisible unless you vary the convention on purpose.

### 6. Value added — against the right baselines

**Measure.** Three comparisons, all on identical splits:

  * against an **image-only** meta-analysis, where images exist. Adding coordinates must not make
    measures 1 to 3 worse. This is a floor, not an aspiration, and three combination schemes have
    already failed it.
  * against the **per-coordinate information ceiling**: regress the held-out truth at a focus on
    the effect size that focus's own table reports. A pooled map should beat it; a pooled map that
    does not is extracting less than a single coordinate carries.
  * against **a convergence count** -- studies reporting within a radius -- which is the crude
    thing the estimator must justify itself against, and which currently wins on some corpora.

---

## Metrics that lie, and what they lie about

Each of these produced a wrong conclusion in this project.

**Root mean square error, on a mostly-zero truth.** Rewards shrinking toward zero. It was the only
column that improved when coordinate gating destroyed 0.36 of area under the curve. Never report
it without a slope or a stratified ratio beside it.

**The ratio of means.** Averages a large over-estimate in the empty majority against the region of
interest. Can sit at 1.00 while the map is flat.

**Correlation over all covered voxels.** Spatial smoothness generates it, and it grows with the
number of contributing studies independently of accuracy.

**Any score restricted to the top quartile.** Hides compression by construction, since the top
quartile is where a compressed map happens to be right.

**"Voxel FWE" without saying which.** The familywise rate and the mean share of voxels rejected
differ by orders of magnitude. Conflating them hid an anticonservative rate of 0.18 behind a
printed 0.0000.

**A correlation against a reference built from the studies being scored.** Manufactures
improvement. It made coordinates look beneficial alongside images when they are not.

---

## How to measure, not what

These are preconditions. A number violating one is void, not weak.

1. **The reference is independent.** Held-out subjects, or split studies with no study on both
   sides.
2. **Nothing is capped.** The number of foci is whatever survives the correction.
3. **Coordinates come from a reporting pipeline**, not a bare uncorrected height.
4. **Every corpus is reported**, including the ones where the estimator does badly. The ordering
   across corpora is itself a finding.
5. **Equal opportunity, not merely equal handicap.** Giving every arm 24 basis functions looked
   neutral and was not: the cap bound hardest on the arm with the most to express, and reversed
   the conclusion when lifted. Before comparing, check that the shared constraint is not binding
   differentially.
6. **Every comparison carries an arm with a known answer** -- an oracle, or a configuration whose
   result is already established at full scale. Four defects in this session were in the harness
   rather than the model, and every one was caught by such an arm.

---

## What measure 6 was missing

The first version of measure 6 said only that adding coordinates must not make measures 1 to 3
worse. That is a floor that licenses *permitting* coordinates; it does not justify *building* a
joint model. It would pass a joint fit scoring 0.948 against the images' 0.946 while saying
nothing about whether the complexity was earned. Three additions.

### 6a. The exchange rate

**Measure.** Sweep images against coordinate studies and read the iso-accuracy contours: how many
coordinate studies reach the accuracy of one more image.

**Why.** This is the decision a user actually faces -- chase one more shared map, or harvest
twenty more published tables. "Does the joint win" is not that question. The first estimate from
the testbed is roughly one image to fifteen coordinate studies, from 30 coordinates reaching
r 0.934 against 2 images at 0.946.

### 6b. Uncertainty, not only accuracy

**Measure.** Interval coverage and median interval width, alongside the point-estimate measures.

**Threshold.** Coverage within [0.93, 0.97] of nominal 95%, with width reported.

**Why.** Every measure above scores a point estimate. A second source may buy precision rather
than accuracy -- thirty-two studies should give narrower intervals than two or thirty even at an
identical point estimate -- and nothing here would see it. This is the most likely place for a
joint model's real benefit to be hiding, and it is currently unmeasured.

### 6c. Domain coverage

**Measure.** The fraction of the analysis volume carrying a usable estimate, by source.

**Why.** The one thing coordinates uniquely offer is territory no shared map reaches, and every
simulation here gives all arms the same support, so that contribution scores exactly zero by
construction.

### And the framing underneath

"Is the joint better" is a model-selection question and a correlation is not a model-selection
metric. **Held-out predictive likelihood** prices accuracy and uncertainty together and penalises
unearned complexity, and should be the arbiter when the point estimates are close.

---

## What success would look like, in one sentence per claim

  * **As a localiser**: AUC above 0.9 against a held-out reference, beating a convergence count.
  * **As a relative magnitude map**: slope above 0.5, stratum profile monotone, agreeing within
    20% across reporting conventions.
  * **As an absolute magnitude map**: every decile above the truth's median within [0.8, 1.25],
    which nothing measured in this project currently achieves.
  * **As an inference procedure**: familywise rate within two standard errors of nominal in every
    null including clustered and heterogeneous ones, with power reported beside it.
  * **As an addition to an image-based meta-analysis**: measures 1 to 3 no worse than the images
    alone, a stated exchange rate against images, narrower intervals at matched coverage, and a
    better held-out predictive likelihood than either source by itself.

The gap between the second and the third is the honest state of the field, not a defect peculiar
to this estimator, and the criteria should be stated so that a method may pass as a relative map
and say so, rather than quietly claiming the third.

## Metrics that lie, continued

**Coverage without width.** A first coverage run reported 100% coverage in all seven arms and
looked like a pass. It was not: in the coordinates-only arm the reported `se` was eight times
the estimator's own sampling sd while the point estimate was 23% high. The interval covered
because it covered almost everything. Coverage is only interpretable alongside the half-width as
a fraction of the effect; an interval with `1.96*se/truth` near or above 1 has certified nothing.
Report the pair, never coverage alone.

**se/sd is the wrong honesty check when the estimator is biased.** The natural diagnostic --
reported `se` against the sd of the estimate across replications -- asks whether the width
matches how much the estimator *moves*. That is the right question only if the estimator is
centred. A systematically biased estimator can be extremely stable: the coordinates-only fit had
sd 0.135 against a bias of +0.575, so `se/sd` of 1.85 looked mildly conservative while coverage
was 0.33. The quantity that predicts coverage is `se` against total error,
`RMSE = sqrt(bias^2 + sd^2)`:

```
                             bias     sd    RMSE     se   se/sd  se/RMSE  coverage
12 studies,  0 images      +0.575  0.135   0.591  0.249    1.85     0.42      0.33
12 studies, 12 images      +0.003  0.049   0.049  0.056    1.16     1.14      1.00
24 studies,  0 images      +0.496  0.088   0.504  0.119    1.35     0.24      0.00
```

`se/sd` ranks these three as 1.85 / 1.16 / 1.35 -- no signal. `se/RMSE` ranks them 0.42 / 1.14 /
0.24, in exactly the order of their coverage. Report both, but read coverage off `se/RMSE`.

**Coverage that gets worse with more studies is a bias signature.** The coordinates-only arm goes
from 0.33 coverage at 12 studies to 0.00 at 24. Nothing about the estimate deteriorated -- the
bias is the same and the sd fell. A fixed bias with a shrinking interval is the one failure mode
that *looks* like improvement on every accuracy metric (RMSE falls, correlation rises) while the
inference gets strictly worse. Any validation suite that reports only accuracy will miss it, and
this one did for months.

## Reflection: what success can still mean for `g`

The coverage measurements change what is worth aiming at, so it is worth restating the target
rather than carrying on against the old one.

Coverage of `g +/- 1.96*se` for a truth of 0.800, prevalence 1, 12 studies, 100 replications:

```
  images of 12   bias    se     coverage   1.96*se / truth
       0        +0.506  0.179     0.10          0.44
       0 *      +0.514  0.193     0.16          0.47      * peak_bias='per-study'
       2        +0.248  0.162     0.67          0.40
       6        +0.079  0.092     0.85          0.23
```

Three things follow, and none of them is a tuning problem.

**1. The documented remedy for having no images does not work on the bias.** `peak_bias`
`'per-study'` moved the bias from +0.506 to +0.514 -- nothing. This is not a contradiction of the
docstring, which says in terms that the per-study correction removes the part of the bias that
*varies between studies* and that fixing the common scale still needs images. It is a number for
that claim, and the number says the between-study part is negligible next to the common part.
Advice that stops at "use peak_bias" is advice to do nothing.

**2. The interval never reaches nominal coverage in any genuinely coordinate-based
configuration.** At six of twelve studies imaged -- a far larger share than real collections
have -- it is 0.85. The "12 of 12" arm does reach ~1.00, but it is not a CBES result at all: donor
coordinates are dropped, so with every study imaged the coordinate table is empty and the fit is
an IBMA wearing CBES's interface (see task #47). The honest reading is that the interval on `g`
has no validated operating point.

**3. Coverage degrades as studies accumulate.** 0.33 at 12 studies and 0.00 at 24, in the
coordinates-only arm of the smoke run. A fixed bias with a shrinking interval gets worse with
more data while every accuracy metric improves. So "more studies" is not a route to a working
interval; it is a route to being confidently wrong.

**So success for `g` cannot be coverage of an absolute magnitude.** Two targets remain available,
and they are not the same:

- *Relative success.* `g` is read as a pattern -- which regions are stronger than which -- and
  success is the correlation and the within-map ordering, with the absolute scale disclaimed and
  `se` either withdrawn or documented as sampling-only. This is only available if the
  coordinates-only bias is *one multiplicative constant*; if the ratio `g/truth` varies with the
  true strength, the pattern is distorted too and there is nothing to read relatively. That is
  what `is_g_a_scale_error.py` decides, and it should be decided before this target is adopted.
- *Marginal success.* `g_marginal = g * prevalence` shares an IBMA's estimand, so it is the one
  output with an external reference and can be held to coverage against one. It currently has no
  standard error at all (task #43), which is the gap to close.

What should be dropped as a target: an absolute `g` with a covering interval from coordinates.
Nothing measured over this program suggests it is reachable, and the reason is not
implementation. A reported peak height carries ~0.055 z of signal per unit of g at a 3.29 cut
(`peak_height_curve`), so the magnitude channel is almost closed before any estimator touches
it. An estimator cannot recover information the reporting practice did not emit.
