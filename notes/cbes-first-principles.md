# What the problem actually is, and what the literature already solved

## 1. The generative process, stated honestly

A study `k` measures a latent effect field `delta(v)` with `N_k` subjects, producing a smooth
statistic field roughly `delta(v) * sqrt(N_k)` plus correlated noise. It thresholds at `u_k` --
by height, by FDR, or by a cluster-extent test -- and tabulates the *surviving topological
features*: one location and one value per peak or per cluster. What reaches a meta-analysis is
that table.

So the question is not "how do I average the reported effect sizes". It is: **what in a table of
thresholded topological features is a sufficient statistic for `delta(v)`?**

## 2. The measurements answer it, and the answer is not the heights

Two independent measurements, both on the current estimator's own inputs.

*Peak height barely moves with the signal.* In `peak_height_curve`, at a cut of 3.29 the mean
`|peak|` goes 3.84 -> 3.90 -> 4.86 as the true field mean goes 0 -> 1 -> 3. A whole unit of
effect buys 0.055 of reported height.

*Peak count moves a great deal.* Over the same sweep the number of peaks per field goes
45.7 -> 123.1 -> 439.3. Roughly tenfold, against 1.26-fold in the height.

`information_ceiling` says the same thing from the real data: regressing held-out truth on the
effect size a table reports at that focus gives a slope of 0.08 to 0.18, so one coordinate
explains 5% to 9% of the variance in the effect at its own location.

**The information about `delta` is in the intensity of the point process, not in its marks.**
That is forced by the data, not chosen. Every failure recorded in `cbes-measurements.md` follows
from having built the likelihood on the marks: the compression, the convention-dependence, the
seven rejected corrections, and the fact that `pi * mu` validates better than `mu` -- prevalence is
a crude nonparametric estimate of the intensity, which is the quantity that carries the signal.

## 3. Random field theory supplies the missing link

The mapping from effect size to peak rate is known. Cheng and Schwartzman derived explicit
formulas for `E[M_u]`, the expected number of local maxima above `u`, for smooth isotropic
Gaussian fields with zero mean; Zhao, Cheng and Schwartzman extend them to *non-zero rotationally
symmetric mean functions*, which is precisely the signal-bearing case
(<https://arxiv.org/abs/2301.04830>).

That gives `E[M_u] = g(delta(v) * sqrt(N_k), u_k, smoothness)`, monotone in the signal and
therefore invertible. Three things the current model fights become model components:

  * **`sqrt(N_k)`** enters the link, rather than being a nuisance that makes the reported
    magnitude scale like `(threshold + overshoot) / sqrt(N)`;
  * **`u_k`** enters the link, so FDR versus family-wise versus cluster-extent thresholding is a
    known offset rather than the uncontrolled 0.82-to-2.29 swing measured on one collection;
  * **smoothness** enters the link, and is estimable from the foci themselves.

The current censoring term instead models the reported *height* as a truncated Gaussian. That is
the wrong observation model for a local maximum, which is why fitting it returned 0.26 for a true
0.5 and 3.6 for a true 2.0.

## 4. This is necessarily a spatial, not a mass-univariate, model

An intensity cannot be estimated at a voxel in isolation: a study contributes zero or one focus
there. The rate only exists by borrowing across space, which means a spatial basis -- splines, a
Gaussian process, or a log-Gaussian Cox process -- and a joint fit. A mass-univariate estimator,
which is what CBES is, cannot represent the quantity that carries the signal.

## 5. What the literature has already built

| work | what it does | what it does not do |
| --- | --- | --- |
| **CBMR** (Yu, Lobo, Riedel, Bottenhorn, Laird, Nichols, *Biostatistics* 2024; multi-group version in *Imaging Neuroscience* 2025) | Generative coordinate-based meta-regression. Spline-parameterised log-intensity, Poisson / negative binomial / clustered NB / quasi-Poisson to handle over-dispersion, study-level covariates. **Already shipped inside NiMARE** as `nimare.meta.cbmr`. | Estimates intensity, not effect size in interpretable units. No selection model linking intensity to `delta`. |
| **Bayesian LGCP regression** (Samartsidis et al., *JRSS-C* 2019) | Foci as a log-Gaussian Cox process, log intensity as a Gaussian process, with regression on study covariates. | Same: intensity, not effect size. |
| **BSPC / hierarchical point process** (Kang, Johnson, Nichols, Wager) | Bayesian spatial point classification; multi-type hierarchical point process for meta-analysis. | Same. |
| **SBLFR** (Montagna, Wager, Barrett, Johnson, Nichols, *Biometrics* 2018) | **Joint CBMA and IBMA.** Foci as a doubly-stochastic Poisson process; images as noisy observations of the same underlying smooth spatial function; sparse latent factor structure for parsimony. | Framed around consistent activation and task prediction, not calibrated effect size. |
| **Selective peak inference** (Davenport and Nichols, *NeuroImage* 2019) | Resampling correction for the selection bias at local maxima, giving unbiased raw units, Cohen's *d* and partial *R*. Explicitly separates "survived a threshold" from "is a local maximum". | Needs the full statistical image, so it cannot be applied to a coordinate-only study. |
| **ES-SDM** (Radua et al., *Eur Psychiatry* 2011) | Combines statistical parametric maps with peak coordinates on an effect-size scale. | Imputes unreported voxels, which is the thing this project rules out, and the scale is asserted rather than earned. |

The architecture in the requirements document -- images carry the magnitude, coordinates carry
the consistency -- is essentially Montagna et al.'s joint model, published in 2018. That is
reassuring about the direction and it means the build is an adaptation rather than an invention.

## 6. The gap worth filling

Nobody has put the **selection mechanism** inside the intensity of a **joint** image-plus-
coordinate model whose latent field is in **effect-size units**.

  * the point-process literature models where foci fall, deliberately not how large the effect
    is, so it never needs `u_k` or `N_k` in a link function;
  * the effect-size literature (ES-SDM) gets units by imputing;
  * the selective-inference literature (Davenport and Nichols) debiases correctly but needs the
    image.

A model that closes this would be:

1. a latent field `delta(v)` on a spline or Gaussian-process basis, in Hedges' *g*;
2. **image studies** as Gaussian observations of `delta(v)` with known variance -- these anchor
   the scale, which is why the two-image floor is the right scope;
3. **coordinate studies** as an inhomogeneous point process with
   `log lambda_k(v) = log rho(v) + log E[M_{u_k}](delta(v) * sqrt(N_k), smoothness_k)`,
   over-dispersed per CBMR's finding that plain Poisson does not fit;
4. inference via CBMR's existing machinery -- penalised splines, parametric bootstrap when foci
   are sparse, which its authors found necessary below about 200 foci per group.

## 7. What this changes about the requirements document

It resolves the incoherence flagged in A3.3. The requirements said coordinates should stop
feeding the magnitude, and then noted that this leaves them contributing nothing to a null
defined on the magnitude map. Under an intensity model that dilemma disappears: coordinates
inform `delta` through *where and how often* they appear, which is real information about the
magnitude, while contributing nothing through their heights, which is where the bias lives.
Coordinates and images become two observation models on one latent field rather than two
competing measurements of the same number, and the scale mismatch that `peak_bias` exists to
patch stops being a free constant -- `u_k` and `N_k` are in the link.

## 8. Honest costs

This is a much larger build than the current estimator, and three things could sink it.

**The link needs `u_k` and per-study smoothness.** Papers report thresholds unevenly and
smoothness almost never. Smoothness is estimable from the foci, but that is an assumption to
validate, not a given.

**Cluster-extent reporting breaks the peak-count link.** One focus per surviving cluster means
the observed count is a count of *clusters*, not of maxima, so `E[M_u]` is the wrong functional
and a cluster-size distribution is needed instead. Measured on pain, this is the common case.

**Over-dispersion is real.** CBMR needed negative binomial and clustered negative binomial
rather than Poisson, which means the intensity is not the whole story and study-level clustering
of foci has to be modelled too.

None of these is a reason not to try. They are the things a validation plan has to cover, and
Part B of the requirements document already covers most of them.

---

## 9. Testing section 2, and finding it half wrong

`intensity_vs_height.py` builds three signals from the *same* coordinate tables and scores each
against held-out HCP subjects: the kernel-weighted mean of the reported effect sizes, which is
what `g` is built from; the number of distinct studies with a focus within a radius, which is the
intensity and what `prevalence` approximates; and the raw foci count.

| design | scheme | signal | r, covered | r, top decile |
| --- | --- | --- | --- | --- |
| 20 x 16 | cluster | height | +0.226 | +0.540 |
| 20 x 16 | cluster | study count, 10 mm | **+0.604** | +0.552 |
| 30 x 12 | cluster | height | +0.305 | +0.504 |
| 30 x 12 | cluster | study count, 10 mm | **+0.544** | +0.520 |
| 20 x 16 | FDR | height | **+0.625** | +0.749 |
| 20 x 16 | FDR | study count, 10 mm | +0.422 | +0.024 |
| 30 x 12 | FDR | height | **+0.647** | +0.707 |
| 30 x 12 | FDR | study count, 10 mm | −0.011 | +0.114 |

**Section 2's claim is regime-dependent, not general.** Under cluster-extent reporting, about six
foci per study, the count beats the height by roughly 2.5-fold, as claimed. Under FDR, two to
three hundred foci per study at a cut near 2.7, the height wins decisively and the count collapses
to zero or below.

The mechanism is saturation. Sparse, selective reporting places a focus only where the effect is
strong, so location carries the signal while the heights sit near the peak of whatever cluster
they came from and carry almost nothing. Permissive reporting places foci nearly everywhere, so
the count stops discriminating effect from region size, while the heights now span a real range.

The corrected statement is not "the information is in the intensity" but **which channel carries
the information is set by the detection function**. That is the quantity Ogata and Katsura
estimate rather than assume, and it is the argument for fitting a detection curve instead of
choosing a channel: the curve says how to weight the two.

A second result needs following up. Under cluster reporting a plain count of studies reporting
within 10 mm reaches +0.604, against +0.427 for CBES's `g` and +0.542 for `pi * g` on the same
design. The splits differ, so this is approximate, but the gap is large enough that a matched
comparison is owed -- and if it holds, the magnitude map is being beaten by a convergence map on
the very quantity it exists to estimate.

---

## 10. A one-dimensional testbed, and what it caught in three minutes

`micro.py` puts the whole question on a 4096-point line: a latent effect field, twenty studies
that each observe it with noise, threshold it and report maxima, and a truth that is known
exactly rather than held out. Three tests run in under three seconds together. A whole-brain
equivalent takes tens of minutes and needs HCP downloads.

The loop immediately earned its keep by finding two defects **in the test code**, one of which
was misdiagnosed first.

*Defect one, real but not the culprit.* The noise was normalised by each realisation's own
empirical standard deviation, which gives the field a Student-t marginal on the effective degrees
of freedom. Harmless near the centre, wrong in the tail, and the tail is where every question
here lives. Fixed by normalising with the filter's own coefficients. It did not change the
calibration, so the first diagnosis was wrong.

*Defect two, the actual culprit.* The truncated-normal check fed `|z|` into a likelihood written
for signed values. At a true mean of 0.5 and a cut of 3.29, 2.8% of exceedances are negative and
folding them onto the positive side shifts the observed mean by about 0.19; at a true 2.0 there
are essentially none. That asymmetry was the entire discrepancy -- the estimator returned 0.80
for a true 0.5 and 1.98 for a true 2.0. The equivalent 3-D script used signed values throughout,
so its results stand.

### Calibration: partial, and the failure is informative

| true mean | cut | MLE on peaks | MLE on draws |
| --- | --- | --- | --- |
| 0.5 | 3.29 | 0.500 | 0.519 |
| 2.0 | 3.29 | 2.466 | 2.009 |

The draws recover the truth, so the likelihood and the field are right. The peaks over-estimate at
high signal, as in 3-D. But 3-D gave **0.26** for a true 0.5 where this gives 0.50, so the
under-estimation at low signal does not reproduce. A 1-D local maximum need only beat two
neighbours; a 3-D one competes with twenty-six, so selection is far weaker here. Anything this
testbed says about the low-signal regime is therefore unvalidated and should be checked in 3-D.

### T1: it is cluster collapsing, not density

| rule | cut | foci/study | r, height | r, count |
| --- | --- | --- | --- | --- |
| every maximum | 4.50 | 7.3 | **0.868** | 0.700 |
| every maximum | 5.00 | 5.3 | 0.571 | **0.750** |
| one per cluster | 3.50 | 9.4 | 0.794 | **0.859** |
| one per cluster | 4.50 | 5.5 | 0.538 | **0.782** |
| one per cluster | 5.00 | 4.1 | 0.590 | **0.859** |

Reporting every maximum keeps the height informative as the cut rises. Collapsing each cluster to
one representative destroys it -- 0.93 down to 0.59 -- while the count holds at 0.68 up to 0.86.
This is the mechanism behind the 3-D flip, and it corrects section 9: density is a symptom, the
reporting rule is the cause. A paper that tabulates one row per cluster has thrown away the
magnitude channel before the meta-analysis ever sees it.

### T3: the ecology claim, repeated too strongly, now measured

Negative log-likelihood relative to the best scale, true scale 0.8:

| data | a = 0.4 | a = 0.6 | a = 0.8 | a = 1.0 | a = 1.5 |
| --- | --- | --- | --- | --- | --- |
| coordinates only | 14.0 | **0.0** | 9.8 | 32.8 | 98.6 |
| plus one image | 449.9 | 114.6 | **0.0** | 95.6 | 1204.8 |

The coordinate-only profile is *not* flat, so the scale is not strictly unidentifiable -- the
claim borrowed from the presence-only literature was too strong, as already suspected. But its
minimum sits at 0.6 against a true 0.8, so it is identifiable and **biased**, and the curvature is
shallow. One image moves the minimum onto the truth and sharpens it by roughly an order of
magnitude.

That is a better argument for the two-image scope than either version before it: images are not
needed to make an unidentifiable problem identifiable, they are needed because the
coordinate-only likelihood is biased and nearly flat, which is exactly what produces an estimate
that tracks the reporting convention instead of the effect.

### A modelling trap found by the micro-test: a free per-study intensity deletes the counts

T4's first version could not recover the effect scale even when handed each study's true
reporting threshold -- 0.616 against a true 0.8 -- which meant the harness was being measured
rather than the idea. The cause was not the crude rate function but the likelihood around it.

An inhomogeneous Poisson log-likelihood is `sum_i log lambda(x_i) - integral lambda`. Write
`lambda_k = C_k * r_k(x)` with a free intensity constant per study, and `C_k` profiles out
analytically to leave a **multinomial over locations**: every count disappears and only the
*shape* of the intensity survives. That is exactly the channel argued in section 2 to carry the
signal, deleted by a modelling choice that looks innocuous.

Sharing one `C` across studies restores it:

| estimator | intensity constant | recovered scale | error |
| --- | --- | --- | --- |
| oracle, true per-study cut | per-study | 0.616 | −0.184 |
| oracle, true per-study cut | shared | 0.688 | −0.112 |
| cut inferred from smallest reported | shared | 0.729 | −0.071 |

The scale is then identified by *how the count changes with each study's threshold and sample
size*, rather than each study's foci count being a free parameter that explains itself.

This is a constraint on any point-process formulation of the problem, including the CBMR-style
frameworks the literature note recommends building on: if study-level covariates are allowed to
absorb each study's overall rate, the model becomes a shape-only estimator and cannot see effect
size at all. It will still fit, and it will look right.

The oracle arm remains biased at 0.688, so the high-threshold rate approximation
`exp(-(u - m)^2 / 2)` is still mis-specified and the soft-versus-hard comparison in that table is
not yet readable. Replacing it with a real expected-maxima formula is the next step before T4 can
answer the question it was built for.

---

## 11. Option C works, and the rate function was never the blocker

T4 stalled because I treated "write the expected-maxima rate" as a hard derivation. It is not:
the testbed contains the exact field generator, so the peak-height law can be **measured**. Draw
zero-mean fields, take every interior local maximum, keep the sorted heights, and the density of
reported maxima where the field's mean is `m` and the cut is `u` is

    rho_max * [ Sbar(u - m) + Sbar(u + m) ]

with `Sbar` the survival of that height law, the upper term for maxima clearing `+u` and the
lower for minima clearing `-u`. Exact for this field rather than a high-threshold approximation,
so a failure downstream belongs to the architecture and not to the rate. In three dimensions the
same law comes from Cheng and Schwartzman, with the non-zero-mean extension of Zhao, Cheng and
Schwartzman, or from the data's own smoothness -- measuring it here is a shortcut for the
testbed, not one that hides a hard step.

With that in place, T5 fits one latent scale from images as Gaussian observations and coordinates
as a point process, with the intensity constant shared across studies:

| images | coordinates | estimator | recovered scale | error |
| --- | --- | --- | --- | --- |
| 2 | 12 | images only | 0.790 | −0.010 |
| 2 | 12 | **coordinates only** | **0.833** | **+0.033** |
| 2 | 12 | joint | 0.791 | −0.009 |
| 2 | 30 | coordinates only | 0.836 | +0.036 |
| 5 | 12 | coordinates only | 0.806 | +0.006 |

**Coordinates alone recover the scale to within 1% to 5%.** Every scheme tried before this had
them two- to fourfold out. The difference is that this reads the intensity -- where and how often
foci appear, given each study's threshold and sample size -- and never touches the reported
value. It is the first result in this project where a coordinate-only magnitude is not badly
biased, and it contradicts the pessimism the rest of these notes built up from CBES's behaviour.
That pessimism was about modelling the marks, and was right about the marks.

Three limits before this is over-read.

*It estimates one scalar with a known spatial shape.* That is far easier than recovering the
field, which is the problem that matters. T6 tests the field.

*The joint does not beat images alone*, 0.791 against 0.790, because two images already pin a
single scalar. Any gain from joining has to appear in the field, where images are absent.

*One dimension understates selection*, per the calibration check in section 10, so the numbers
will not transfer even if the direction does.

### Smoothness is not a blocker: the shared intensity constant cancels it

The intensity needs a peak-height law, the law depends on smoothness, and no method exists for
estimating smoothness from coordinates alone -- a real gap in the literature. T7 asks whether the
gap matters, by generating at one smoothness and fitting with a law measured at another.

| assumed sigma | versus truth | maxima per point | recovered scale | error |
| --- | --- | --- | --- | --- |
| 2.0 | x0.5 | 0.0954 | 0.867 | +0.067 |
| 3.0 | x0.75 | 0.0643 | 0.867 | +0.067 |
| 4.0 | x1 | 0.0485 | 0.874 | +0.074 |
| 6.0 | x1.5 | 0.0324 | 0.870 | +0.070 |
| 8.0 | x2 | 0.0243 | 0.873 | +0.073 |

A fourfold error in smoothness moves the scale by 0.007. This is structural rather than lucky.
`rho_max` enters the intensity multiplicatively and the shared constant is profiled as
`C = n_total / integral`, so `C * rho_max` is invariant and the maxima density cancels
identically. What is left is the *shape* of the peak-height survival, which for a standardised
Gaussian field depends on the dimensionless ratio of spectral moments rather than on the full
width at half maximum.

So the model does not need a smoothness estimate. That is the second payoff from sharing the
intensity constant across studies, the first being that a free per-study constant profiles the
counts away entirely.

The caveat is narrow and worth keeping. Sigma was varied within one kernel family, where the
autocorrelation *shape* is identical by construction, so what is shown is invariance to the scale
of smoothness and not to the shape of the spectrum. A genuinely different spectral shape --
non-Gaussian smoothing, unsmoothed data -- would still matter. That is a far weaker requirement
than knowing each study's FWHM, and far more stable across a literature.

## 12. How this differs from CBMR, stated correctly after being stated wrongly

CBMR's predictor, read from `nimare/meta/cbmr/predictor.py`, is separable:
`log_intensity_by_pattern` returns a spatial log-intensity and `moderator_effect` returns "the
scalar linear predictor contributed per experiment", combined as
`lambda_k(x) = exp(spline(x)) * w_k`. One spatial shape, scaled per study.

My first statement of the difference was that non-separability -- each study having a *different*
intensity shape, because `u_k` and `sqrt(N_k)` enter inside the nonlinearity -- is what identifies
the effect scale, the way varying detection across observers identifies the Gutenberg-Richter law.
T8 tests that by removing the variation.

| collection | recovered scale | error | curvature of the profile |
| --- | --- | --- | --- |
| homogeneous, one threshold and one sample size | 0.762 | −0.038 | 368 |
| varying sample size only | 0.775 | −0.025 | 491 |
| varying sample size and threshold | 0.817 | +0.017 | 565 |

**The homogeneous collection still identifies the scale.** Variation is not necessary. It helps --
the profile is 54% sharper and the bias flips from −0.038 to +0.017 -- but the mechanism is
something simpler: with a *known* threshold and sample size, the expected count above that
threshold is already a monotone function of the effect, so counting identifies it even when every
study is identical. Cross-study variation buys precision, not identifiability.

The difference from CBMR survives, restated. It is not separability as such; it is that CBMR's
link contains no threshold and no sample size, so the effect field is not a parameter of the model
at all and there is nothing to invert back to. That is deliberate: CBMR answers where foci occur
and whether that depends on a covariate, and answers it well. The proposal keeps its likelihood
family, its spline bases, its negative binomial and clustered negative binomial handling of
over-dispersion and its inference, and replaces a free log-intensity spline with an effect field
pushed through a known reporting mechanism.

One CBMR feature would actively break it: study-level covariates give each experiment a free
multiplicative rate, and sample size is among the covariates its paper names. That absorbs the
counts, which is the channel the effect scale rides on.

## 13. The field test, and why it corrects section 11 and the smoothness result

T5 recovered the effect *scale* from coordinates alone to within 1% to 5%, with the spatial shape
known. T6 frees the shape -- 24 basis coefficients -- which is the problem that matters, and the
result does not carry over: coordinates alone reach a shape correlation of +0.44 with the scale
inflated 5.7-fold, and the joint fit is indistinguishable from the images by themselves, +0.689
against +0.686.

Three probes to find out why, two of which refuted their own hypothesis.

*Not the optimiser.* Started at the truth the fit walks away, to a ratio of 2.79, and reaches a
**better** negative log-likelihood than the truth has -- 2834 against 2968. The likelihood prefers
an inflated field, so the optimiser was doing its job.

*Not cluster-collapsed reporting.* With every study reporting every maximum the ratio is still
2.43. Collapsing makes it worse, 3.43, but is not the cause.

*It is the profiled intensity constant.*

| intensity constant | start | r | ratio |
| --- | --- | --- | --- |
| profiled out | truth | +0.681 | 2.43 |
| profiled out | flat | +0.328 | 4.25 |
| fixed from the maxima density | truth | +0.701 | **1.76** |
| fixed from the maxima density | flat | +0.751 | **1.77** |

Setting `C = n_total / integral` makes the likelihood depend only on the *shape* of the intensity
and throws away the absolute expected count. With a known spatial shape that costs nothing, since
a single amplitude cannot hide. With a free field it is fatal: the coefficients reproduce any
intensity shape at any amplitude. Fixing the constant cuts the inflation to 1.76 and, unexpectedly,
makes the fit start-independent -- the profiled version reached entirely different optima from a
flat and a truth start, so it was badly conditioned as well as unidentified.

**This withdraws the smoothness result in section 10's follow-up.** I concluded the model needs no
smoothness estimate because the maxima density cancels against the profiled constant. That
cancellation is real and it *is* the loss of level identification. The two are one knob:

  * profile the constant -- smoothness does not matter, the effect level is unidentified;
  * fix it from an estimated density -- the level is identified, a smoothness estimate is required.

T7 looked benign only because it fitted one amplitude with a known shape, where the level could
not hide. On a free field the trade-off bites, and the requirement to estimate smoothness comes
back. In the two-or-more-images scope that is answerable from the images, which is one more thing
the scope buys.

A residual inflation of 1.76 remains with the constant fixed. Candidates not yet separated: the
basis represents the truth only to r = 0.679, the ridge penalty is arbitrary, and the rate
function is still an approximation at the cluster-reporting studies.

### T6 complete: the joint never beats the images

| images | coordinates | estimator | r | rmse | ratio |
| --- | --- | --- | --- | --- | --- |
| 2 | 30 | images only | +0.686 | 0.157 | 1.25 |
| 2 | 30 | coordinates only | +0.442 | 0.983 | 5.74 |
| 2 | 30 | joint | +0.689 | 0.156 | 1.29 |
| 5 | 30 | images only | +0.694 | 0.155 | 1.22 |
| 5 | 30 | coordinates only | +0.501 | 0.570 | 4.16 |
| 5 | 30 | joint | +0.695 | 0.155 | 1.23 |

Adding thirty coordinate studies to two images moves the shape correlation from +0.686 to +0.689.
The point-process formulation does not repeat the *harm* the earlier schemes did -- pooling,
rescaled pooling and gating all made things measurably worse, and this does not -- but it does not
help either. On this corpus, with these study counts, coordinates modelled correctly are neutral
rather than beneficial once images are present.

That is a meaningfully better place to be than "actively harmful", because a neutral channel can
become useful where images are absent, which is the case none of these designs can test. But it
is not the result section 11 was heading toward.

### The basis was the binding constraint, and T6's verdict is suspended

A sweep of the penalty and the basis, with the intensity constant fixed:

| basis | penalty | strength | ceiling r | r | ratio |
| --- | --- | --- | --- | --- | --- |
| 24 | ridge | 0.01 | 0.679 | +0.751 | 1.77 |
| 24 | roughness | 1 | 0.679 | +0.675 | 1.65 |
| 24 | roughness | 30 | 0.679 | +0.532 | 1.82 |
| 24 | roughness | 300 | 0.679 | +0.472 | 1.96 |
| 48 | ridge | 0.01 | **0.954** | **+0.936** | **1.47** |
| 48 | ridge | 1 | 0.954 | +0.934 | 1.46 |
| 48 | roughness | 1 | 0.954 | +0.922 | 1.47 |
| 48 | roughness | 30 | 0.954 | +0.837 | 1.52 |
| 48 | roughness | 300 | 0.954 | +0.677 | 1.80 |

The roughness hypothesis is refuted: it degrades the correlation and does not fix the inflation,
so CBMR's penalty is not the answer to this particular problem. The basis is: doubling it lifts
the achievable correlation from 0.679 to 0.954 and the achieved one from 0.751 to 0.936, while the
inflation falls to 1.47.

**This probe contains no images at all, so +0.936 is coordinates alone** -- against the +0.442 and
+0.501 T6 reported for the same arm. T6 gave every arm 24 basis functions and profiled the
constant, and both choices bit hardest on the arm carrying the most spatial information to
express. Its conclusion that the joint never beats the images was measured under a handicap that
was not neutral between arms, and is suspended pending the rerun.

The lesson is the same one the harness bugs kept teaching: a shared handicap is not a fair
comparison when the arms differ in how much they need the thing being limited.

Completing the sweep at the larger basis settles the penalty question: a hundredfold change in
ridge strength moves the correlation from 0.936 to 0.934 and the inflation from 1.47 to 1.46, and
roughness at matched strength is the same. Heavy penalties only hurt. The residual inflation is
**insensitive to regularisation**, so it is not overfitting and CBMR's roughness penalty is the
right tool for a different problem.

Basis size explains the correlation; something else sets the level. The named suspect is the
threshold: these probes infer each study's cut as the smallest reported value, the signal-dependent
floating cut that has caused trouble throughout this project. T5 suggested that is benign when the
shape is known, and the field setting is exactly where it would stop being benign.

## 14. The rerun overturns T6, and with it an argument made against SDM-PSI

With 48 basis functions and the intensity constant fixed:

| images | coordinates | estimator | r | rmse | ratio |
| --- | --- | --- | --- | --- | --- |
| 2 | 30 | images only | +0.946 | 0.089 | 1.38 |
| 2 | 30 | **coordinates only** | **+0.934** | 0.107 | **1.37** |
| 2 | 30 | joint | +0.948 | 0.088 | 1.37 |

Thirty coordinate studies recover the field almost as well as two images, against the +0.442 and
5.74-fold inflation the handicapped version reported for the same arm. The joint is marginally
ahead of the images, +0.948 against +0.946, which is within noise but is at least no longer a
loss.

Two things follow.

**The residual inflation was never the coordinate model's.** It is 1.38 in the images-only arm and
1.37 in the coordinate arm, so it belongs to the basis, the ridge and taking absolute values in the
scoring -- my harness, common to every arm. Three probes were spent chasing it as a property of the
point process.

**The information-ceiling argument was over-extended, including against SDM-PSI.** The measured
ceiling -- one tabulated coordinate explaining 5% to 9% of the variance in the truth at its own
location -- bounds what a single coordinate's *value* carries. It does not bound a method that
reads the whole table, as this project's own point-process model reaching +0.934 demonstrates.
Any argument of the form "no method can recover the magnitude, because the per-coordinate ceiling
is low" is wrong, and it was used here against both CBES and SDM-PSI.

What survives against a method that takes reported peak values at face value is narrower and
still stands: those values carry a shared bias of about +1.45 against 0.6 of noise, which does not
average down. The remedy is to model the selection rather than the value, which is what the
intensity formulation does.

## 15. The exchange rate, and two checks that should have come first

### One image is worth about fifteen coordinate studies

Correlation with the truth, spread across three seeds in brackets, basis ceiling 0.954:

| | 0 coords | 5 coords | 15 coords | 40 coords |
| --- | --- | --- | --- | --- |
| **0 images** | -- | 0.890 [0.012] | 0.927 [0.009] | 0.933 [0.002] |
| **1 image** | 0.926 [0.017] | 0.930 [0.016] | 0.932 [0.020] | **0.945 [0.006]** |

One image alone reaches 0.926; fifteen coordinate studies alone reach 0.927. That is the number a
user actually needs, and it says chasing one shared map is worth roughly fifteen table-harvests --
a different research strategy from the one "coordinates are contaminating" would suggest.

**At one image, coordinates clearly help.** Adding forty takes the correlation from 0.926 to 0.945
and tightens the seed spread from 0.017 to 0.006, three times more stable. That is the precision
benefit predicted to be hiding outside the point-estimate measures, and it is the first time the
joint has beaten a component on anything.

It also explains why T6 saw nothing: it used two and five images, which is past the point where
the benefit has saturated. **Coordinates help in proportion to how few images there are**, which
is the practically relevant regime, since most collections have none to three.

### The null-field check, which should have been first

Every test of this model had been handed real signal. On a field of pure noise it returns a mean
of 0.052 and a maximum of 0.173, against a signal field whose mean is near 0.35. With a
deliberately misspecified rate law, 0.060 and 0.185.

So it does not hallucinate, with either rate law. There is a small positive floor near 0.05,
which is expected -- six or seven noise peaks per study are still reported and the model must
explain them with something -- but it sits an order of magnitude below real signal. This was the
cheapest possible check and it came after every celebratory result rather than before them.

## 16. What the point-process estimand actually is

The model fits one field common to every study, so the question of what that field means when
studies differ was never asked. `estimand.py` asks it: studies either have the effect or have
none, with prevalence `pi`, and those that have it carry a study-level perturbation of size `tau`.

| prevalence | tau | ratio to mu | ratio to pi*mu |
| --- | --- | --- | --- |
| 1.00 | 0.0 | 1.68 | 1.68 |
| 1.00 | 0.5 | 1.70 | 1.70 |
| 0.50 | 0.0 | 1.37 | 2.73 |
| 0.50 | 0.5 | 1.31 | 2.63 |
| 0.25 | 0.0 | **1.28** | **5.13** |
| 0.25 | 0.5 | 1.30 | 5.21 |

Prevalence falls fourfold and the estimate moves only 1.68 to 1.28, where the marginal would have
gone to 0.42. **The point-process field is essentially the conditional magnitude**, with about a
quarter of attenuation at low prevalence -- not a compromise between the two, which is what a
Jensen argument had predicted and what an earlier note claimed outright as "the marginal". Both
were wrong. And `tau` does nothing at all: 1.68 against 1.70, 1.28 against 1.30. Magnitude
heterogeneity does not shift the estimand; only zero-inflation does, and weakly.

So the intensity formulation **does not change the estimand**. It targets the same conditional
quantity CBES already targets, which is the quantity no image-based meta-analysis estimates.

### Which estimand each method targets

| estimand | meaning | estimated by |
| --- | --- | --- |
| `mu` | effect among studies that have one here | CBES `g`, the point-process field |
| `pi` | fraction of studies with an effect here | CBES `prevalence` |
| `pi * mu` | effect averaged over all studies, zeros included | CBES `g_marginal`, **and every IBMA** |
| foci density | rate of reported peaks | ALE, MKDA, CBMR |

Every IBMA estimator in NiMARE -- DerSimonian-Laird, Hedges, weighted least squares, the
likelihood estimators -- pools per-study effect maps as draws around one pooled mean, so a study
with no effect at a voxel enters the average as a zero and the expectation over studies is
`pi * mu`. On a real Hedges' g scale, since each per-study `g` comes from an image with a known
sample size.

Two consequences.

**`g_marginal` shares its estimand with an IBMA**, differing only in an unidentified scale. That
is why it validated best against held-out references all session: not a better estimator, the only
CBES map pointed at the same quantity as the reference.

**`g` has no image-based reference at all.** No IBMA estimates the conditional. Worse, the
conditional is arguably not identifiable from images either, because computing it needs each study
classified as having-an-effect-or-not at each voxel, which is a thresholding decision and
reintroduces the selection the estimator exists to correct.

Part of the apparent inflation of `g` is therefore mismatch rather than error, worth `1/pi`, which
on pain is 1.4 to 1.8. It does **not** explain the compression: measured prevalence rises 0.556 to
0.702 across strata, so `1/pi` falls only 1.80 to 1.42, a factor of 1.27 against an observed ratio
fall of 9.4.

### What can be claimed

*Supported.* Relative statements within one map, shape correlation 0.89 to 0.94 across every
condition including misspecified rate laws and low prevalence. Localisation. "When a study finds
this region, how large is the effect it finds", which is what the conditional means. And --
the one no other estimator offers -- **predicted reporting**: given a threshold and a sample size,
how many foci a study should report near a voxel, which is literally what the model is fit to and
is directly usable for planning, replication judgements and simulating coordinate data.

*Unsupported.* An absolute Hedges' g comparable to an IBMA, which is a category error before any
bias enters. Cross-collection magnitude comparison when thresholds or sample sizes differ, since
the estimand's definition moves with the rate function's curvature. Separating prevalence from
magnitude. Power calculations for a new study, which need the marginal on a real scale.
