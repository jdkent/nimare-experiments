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
