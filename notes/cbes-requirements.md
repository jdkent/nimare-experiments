# CBES requirements, written as if from scratch

Scoped to a collection carrying **at least two image studies** alongside any number of
coordinate-only studies. No coordinate-only mode and no image-only mode. That scope is not a
convenience: it is what makes the central quantity identifiable at all, and it changes which
component is responsible for what.

Everything below is either a requirement or a test of one. Where a requirement is already known
to fail, it says so, with the measurement.

---

## Part A: what the estimator must do

### A1. Input contract

**A1.1** Accept a collection of studies, each contributing *either* a group-level effect-size
image (`g` and `g_var`, or a statistic plus sample size that converts to them) *or* a table of
reported coordinates with a statistic and a sample size.

**A1.2** Require at least two image studies. Refuse with a named error otherwise, not a warning
and not a degraded mode. Two rather than one because the calibration constant estimated from a
single donor has no measurable uncertainty, and a constant reported without its uncertainty will
be read as exact.

**A1.3** Accept a declared analysis mask per study, and a declared reporting threshold per
study. Both are things a paper states and the estimator otherwise has to guess.

**A1.4** Refuse, rather than silently return, any input from which the requested output cannot
be computed: no in-mask voxels, fewer than two images, a permutation null with too few distinct
arrangements.

### A2. Estimands, named before they are estimated

**A2.1** Three distinct quantities, never conflated in one map:

| symbol | meaning | which studies inform it |
| --- | --- | --- |
| `mu(v)` | effect size among the studies that have an effect at `v` | images set the scale, coordinates the pattern |
| `pi(v)` | fraction of studies with a non-null effect at `v` | all studies |
| `pi(v) * mu(v)` | effect averaged over all studies, including those with none | both of the above |

**A2.2** Every emitted map must state, in its own documentation, which estimand it is, what its
units are, and which studies contributed to it. A reader must not have to infer from a name
whether a number is on the Hedges' *g* scale.

**A2.3** The prevalence is a parameter of the model, not a by-product of spatial density. It
must be estimable separately from the magnitude and must not be a monotone function of the
number of nearby foci.

### A3. Where the magnitude comes from

**A3.1 The absolute magnitude is carried by the image studies.** Coordinate studies must not
move the magnitude scale.

*Evidence.* On the 21-study NIDM pain collection, split so the reference comes from studies the
coordinates never touched, adding coordinate studies to `k` images costs correlation with the
held-out truth at every `k` tried: +0.520 to +0.320 at one image, +0.723 to +0.589 at two,
+0.774 to +0.602 at three, +0.829 to +0.662 at five, all with the coordinates rescaled onto the
image scale. Root mean square error is worse at three counts of four. One image study beats nine
coordinate tables *plus itself*.

*Why no weighting scheme fixes it.* The per-coordinate error is a bias, not noise: the reported
`g` at a focus is 1.965 +- 0.593 against a held-out truth there of 0.519 +- 0.312, so the error
is a +1.45 offset against about 0.6 of spread, and the empirical error variance is only about
twice the nominal variance the model assigns. Noise averages out as studies are added. A bias
every coordinate study shares does not -- pooling more of them converges harder on the wrong
number.

**A3.2 Coordinate studies carry the evidence of consistency.** They inform `pi`, the coverage
that decides where anything is estimated, the between-study heterogeneity, and the permutation
null. This is what they are good for and it is not nothing: a collection of three images and
thirty tables knows far more about *where* an effect is reliable than the three images alone.

**A3.3** It follows that the null must be a null for a statistic the coordinates actually move.
A permutation null defined on a magnitude map that only the images determine is a null in which
shuffling coordinates changes nothing, and those studies contribute no inference. Whatever
statistic is tested must depend on the coordinate studies. (An earlier proposal of mine --
"coordinates feed the null but not the magnitude" -- is incoherent for exactly this reason.)

### A4. Non-reporting

**A4.1 Nothing is imputed.** Not a map, not a value, not a distribution over values. A study
that reported nothing near a voxel contributes the *probability* of having reported nothing,
and only that.

**A4.2** Silence must be distinguished from not having looked. A study that declared an analysis
mask, or that is an ROI study, contributes neither a value nor a silence outside what it
examined.

**A4.3** Silence must be distinguished from a zero effect. The mixture must allow "this study
has no effect here" and "this study has an effect too small to have been reported" to be
separate explanations, with the data choosing between them.

### A5. Invariance requirements

These are the requirements the current implementation most clearly fails, and they are the ones
worth designing against from the start.

**A5.1 Reporting-convention invariance.** Two collections reporting the same underlying effects
must produce the same magnitude, whether their authors thresholded by FDR, by voxelwise
family-wise error, or by cluster extent, and whether they tabulated each cluster by its peak or
by its centre of mass.

*Currently fails.* Holding one collection fixed and changing only the convention, the ratio of
`g` to the held-out truth where that truth is largest runs from 0.82 to 2.29.

**A5.2 Compression.** The estimate must respond to the truth. Across strata of a held-out
reference spanning elevenfold, the estimate must span materially more than it currently does.

*Currently fails.* It spans about 1.2-fold. Seven independent interventions left this untouched
(a truncated-normal selection correction, the per-study peak-bias rescaling, subtracting the
censoring floor, realistic reporting, centre-of-mass foci, oracle cluster-extent coverage, and a
widened assumed cluster), which is why A3.1 puts the magnitude on the images instead.

**A5.3 Scale invariance of the relative map.** Whatever map is designated the default read must
be invariant to the unidentified constant.

**A5.4 Exchangeability.** The null must hold fixed every nuisance structure that is not the
hypothesis: study membership, focus positions, sample sizes, per-study precision.

**A5.5 Determinism.** A given seed and input yield a bit-identical result.

### A6. Uncertainty

**A6.1** Report the calibration constant with an interval derived from the image donors, not as
a point.

**A6.2** Report a standard error that accounts for the mixture -- both which component an
observation came from and the cost of not knowing the prevalence -- and state the reference
distribution to use it with.

**A6.3** Where a voxel's magnitude is not identified, say so rather than returning the point an
optimiser stopped at. A voxel informed by one reported value and several silences has a flat
likelihood, and the returned number is an artifact of the stopping rule.

---

## Part B: what a rigorous validation must do

Part B exists because most of what went wrong in developing this estimator was measurement, not
modelling. Three protocol errors each produced a confident, wrong conclusion that survived for
days. The preconditions below are therefore not advice.

### B0. Preconditions. A result violating any of these is void, not weak.

**B0.1 Never cap the number of foci per study.** No "strongest N", no `.head(n)`. A cap fixes
the count and lets the effective threshold float to that study's Nth peak, which is a function
of its signal, so every threshold-dependent quantity becomes signal-dependent. This invalidated
a conditional-reference measurement, a threshold-correction fit, and a whole cross-dataset run.
The count of coordinates is not a free parameter; it is whatever survives the correction.

**B0.2 Extract coordinates the way papers produce them.** Multiplicity-corrected thresholds
(FDR, voxelwise family-wise, or a cluster-forming cut with a family-wise extent test computed
rather than assumed); only surviving clusters, so a study with nothing left reports nothing and
drops out; one row per cluster, at its peak or its centre of mass. Report results under more
than one convention, because papers do not agree on one and the estimator's output depends on
which was used.

**B0.3 The reference must be independent of the coordinates.** With one map per study there is
nothing to split within a study, so a reference built from the same maps that produced the peaks
conditions the comparison on the noise being measured. Use held-out subjects where per-subject
data exist, or split the studies so no study contributes to both sides.

**B0.4 Distinguish a familywise rate from a per-voxel rate.** `P(any rejection)` and the mean
share of voxels rejected differ by orders of magnitude and are easy to swap. State which is
meant at every reported number. A reported "voxel FWE of 0.0000" turned out to be a per-voxel
share where the familywise rate was about 0.10.

**B0.5 Report the configuration.** A measurement of a mixed collection made with the
calibration switched off is a measurement of a configuration the documentation tells users not
to use.

### B1. Correctness of the machinery, against independent implementations

| test | criterion |
| --- | --- |
| **B1.1** Kernel weights all 1 reduces to DerSimonian-Laird | matches PyMARE to floating-point tolerance |
| **B1.2** Censored mixture likelihood | matches an independently written brute-force maximisation over a grid |
| **B1.3** Analytic score and curvature | match finite differences of an independently written log-likelihood at every cell of a `(pi, mu)` grid |
| **B1.4** Observed-information standard error | interval coverage within [0.93, 0.97] of nominal 95% across prevalences, cutoffs and study counts |
| **B1.5** Determinism | bit-identical output for a repeated seed |

### B2. Recovery of the estimands, in simulation where the truth is set

**B2.1** Simulate the whole reporting process: a true effect field, a study-level draw, a
sampling draw, a within-study threshold deciding whether a peak is reported. A point-based
simulator cannot test anything about peak height and must not be used for it.

**B2.2** The simulated statistic must be on the scale it claims. A *t*-scale noncentrality
written into a field and read as a *z* produced an estimate of 8.8 from a true 1.6 and inverted
the apparent direction of the selection bias.

**B2.3** Vary `pi` and `mu` independently and recover both. Criterion: `mu` within 15% at each
level of `pi`; `pi` monotone in the truth with a rank correlation above 0.9. If only ordinal
recovery is achieved, the map must be documented as ordinal.

**B2.4** Vary sample size across studies by two orders of magnitude. Nothing may depend on the
roster's homogeneity.

### B3. Magnitude, against an independent truth

This is the test the estimator exists to pass, and it needs per-subject data.

**B3.1 Held-out-subject design.** Split subjects: some become synthetic studies, each thresholded
and reduced to a coordinate table or kept as an image; the rest, never used to make a coordinate,
give the truth. Selection is then statistically independent of what it is compared against.

**B3.2 Acceptance criteria**, stratified by the truth, because an average ratio hides everything:

- ratio of estimate to truth within [0.8, 1.25] in every stratum above the truth's median;
- the regression of truth on estimate has a slope above 0.5, so the estimate responds to the
  effect rather than to the threshold;
- correlation with the held-out truth above +0.6 over covered voxels.

**B3.3 Convention invariance (A5.1).** Run B3.1 under at least three reporting conventions.
Criterion: the stratum ratios agree within 20% across conventions. The current implementation
spans 0.82 to 2.29 and fails this.

**B3.4 The coordinate-dilution regression test.** With `k` images fixed, adding coordinate
studies must not reduce correlation with the held-out truth by more than 0.02. This is the test
that motivated A3.1, and it must stay in the suite so a future change cannot quietly reintroduce
the problem.

**B3.5 Ceiling check.** Regress the held-out truth at a focus on the effect size the table
reports there, with no estimator in between. This bounds what any method can extract, and the
estimator should be scored against that bound rather than against perfection. Measured on pain:
slope 0.08 to 0.18, so one coordinate explains 5% to 9% of the variance in the effect at its own
location -- and the pooled map beats it, which is the honest way to report the magnitude result.

### B4. Inference

**B4.1** Familywise error rate, `P(any rejection)`, under a global null, at voxel level and at
each cluster level. Criterion: within two binomial standard errors of nominal. With 100
simulations that is 0.05 +- 0.044.

**B4.2** The global null must include a **spatially clustered** variant, not only independent
noise. Clustered nulls are where the current voxel-level rate reads 0.100 +- 0.039 against
nominal 0.05, while cluster-size and cluster-mass read 0.000 and 0.050. An estimator validated
only on unclustered nulls will not have found this.

**B4.3** Heterogeneous sample sizes, spanning 10 to 1000. This is the condition under which an
across-study shuffle rejected at 96.7% against a nominal 5%.

**B4.4** Power at a known focal effect, reported at voxel and cluster level, alongside the error
rates. An estimator that controls error by never rejecting is not acceptable.

**B4.5** Refusal behaviour: a collection with too few distinct arrangements returns `p = 1` with
a warning, and this is tested, not assumed.

### B5. Real-data behaviour

**B5.1** At least three independent corpora, differing in how much the studies have in common.
Results must be reported for all of them, including the ones where the estimator does badly --
the ordering across corpora is itself the finding.

**B5.2** ROI and partial-coverage studies present in at least one corpus, with the declared
analysis mask honoured.

**B5.3** Mixed conventions within one collection, which is the realistic case and the hardest
one for A5.1.

**B5.4** A collection where the images and the coordinates disagree, to confirm the disagreement
surfaces rather than averaging away.

---

## What this scope buys, and what it does not

Requiring two images makes the absolute scale identifiable and lets A3.1 put the magnitude where
the unbiased data is. It does not fix A5.2: a voxel covered only by coordinate studies is still
estimated from values that carry a shared, threshold-driven bias, and the honest response is to
say which voxels those are rather than to report a number that looks like the others.

The estimator that results is best described as an image-based meta-analysis whose coverage,
consistency and inference are extended by the coordinate literature -- not as a coordinate-based
effect-size estimator that can also read images.
