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

## 17. The one intuition that has survived everything: counts carry the signal, heights do not

Stated as a single claim, because it now organises almost every result in this program:

> In a coordinate table, the information about an effect is in *whether and where* foci appear,
> not in *how large* the reported statistics are. Every method that used counts and locations has
> worked; every method that used reported magnitudes has failed, and failed in the direction the
> selection predicts.

The evidence, gathered for other reasons and pointing the same way each time:

| measurement | what it says |
| --- | --- |
| peak height moves 0.055 z per unit g at a 3.29 cut, while the count moves ~10x | the height channel carries a few percent of what the count channel does |
| one coordinate explains 5-9% of the variance in the truth at its own location | a single mark is nearly uninformative |
| coordinates-only `g` is +63% biased; coverage 0.10 at 12 studies, 0.00 at 24 | the height-based estimator does not converge on the truth, it converges on the threshold |
| `peak_bias='per-study'` moves the bias +0.506 -> +0.514 | the between-study part of the height bias is negligible; the common part is everything |
| adding coordinates to images degrades magnitude (r 0.774 -> 0.614) | height information is not merely weak, it is contaminated, and the contamination is shared so it does not average away |
| the point-process fit recovers the scale to 1-5% and the field to r 0.934 | the count/location channel is strong enough to recover a magnitude *indirectly* |
| `prevalence` is the only output with no external competitor and it is count-driven | what holds up is what counts |

**Why this is not obvious a priori.** A reported peak height looks like a measurement of effect
size -- it is on the right scale, it has a sample size attached, it converts to a Hedges' g. The
trouble is that it is the maximum of a field over a region selected *because* it was large, so
its distribution is pinned near the threshold almost independently of the truth. The floor at
mu = 0 is 3.84 at a 3.29 cut and 3.90 at mu = 0.5; the truth moved by half a standard deviation
and the reported height moved by 0.06. Whereas the *number* of surviving clusters went from 45.7
to 64.6 -- a 41% change. The signal was never in the value.

**What follows for the design.** CBES is built the other way round: a censored likelihood on
reported magnitudes, with locations serving only as support for a smoothing kernel. That places
the model's whole estimating equation in the channel that carries least, and uses the strong
channel only to decide which voxels to touch. The point-process formulation (task #40) inverts
that -- intensity estimated from where and how often foci appear, with magnitudes ignored or used
only as weights -- and it is the only approach here that has recovered an absolute scale from
coordinates alone.

**What follows for the deliverable.** If this is right, `prevalence` is the feature and `g` is a
vestige of the wrong premise. `prevalence` is the count channel expressed directly; its known
faults (a floor near 0.2, compression, non-monotonicity at weak effects, ordering reliable only
on average) are calibration problems in a quantity that is at least identified. `g`'s fault is
that the data barely constrain it. Those are not the same kind of problem and they should not be
presented as equally fixable.

**The check that would falsify it.** A count-only estimator should then beat the height-based
magnitude map at recovering the *pattern* of the truth, and a naive count of studies with a
nearby focus should be competitive with the fitted `prevalence`. Both are queued
(`prevalence_vs_naive.py`, task #38). If the naive count loses badly to the fitted prevalence,
the likelihood is extracting something real from the heights after all and this claim is too
strong.

## 18. What identifies prevalence separately from magnitude: heterogeneous power

Write down what the data at a voxel actually are, once the intuition of section 17 is taken
seriously and the reported heights are set aside as nearly uninformative. For study *k* and voxel
*v* the informative observation is a binary one:

    R_kv = 1 if study k reported a focus within r of v, else 0

and its probability factors into having the effect and detecting it:

    E[R_kv] = pi_v * D(mu_v, n_k, u_k)

with `D` the probability that a study of sample size `n_k` reporting at threshold `u_k` produces
a surviving cluster there given a true effect `mu_v`. This is exactly an **occupancy model with
imperfect detection**: `pi` is occupancy, `D` is the detection function, and the observation is a
detection/non-detection record.

That framing imports a known and decisive result. In an occupancy model, occupancy and detection
are **separately identified only through repeat visits or through covariates that shift detection
while leaving occupancy alone** (MacKenzie et al. 2002 is the canonical statement; it is why
occupancy designs mandate repeat surveys). With a single visit and no detection covariates, only
the product `pi * D` is identified -- any `pi` can be traded against any `D` that preserves it.

Map that back:

| occupancy model | coordinate meta-analysis |
| --- | --- |
| site | voxel |
| occupancy `pi` | fraction of studies with the effect |
| detection probability | probability a study's map survives thresholding there |
| repeat visits to a site | multiple studies covering the voxel |
| detection covariates | sample size `n_k`, reporting threshold `u_k`, smoothness |

**So `pi` and `mu` in CBES are separately identified only by the spread of sample sizes and
reporting thresholds across the studies that cover a voxel.** With a homogeneous roster --
every study 20 subjects at p < 0.001 corrected -- the likelihood has one identified combination
and the split into prevalence and magnitude is determined by the model's functional form rather
than by the data, however many studies there are.

This explains, in one mechanism, four findings that were recorded separately as puzzles:

- **`g_marginal` is the best-behaved magnitude map.** It is the closest thing the estimator emits
  to the identified combination. It was written up as "two biases happen to cancel"; the reason
  the cancellation keeps recurring is that the product is what the data pin down.
- **`prevalence` has a floor near 0.2 and is non-monotone at weak effects.** A weakly identified
  parameter is pushed around by the prior implicit in the likelihood's shape, and non-monotonicity
  is what a ridge in the likelihood looks like when you profile it.
- **The `pi`/`mu` split is unstable across collections while the product is stable.** Textbook
  weak identifiability.
- **Peak heights cannot rescue it.** The height channel carries ~5% of the count channel's signal,
  so it cannot supply the second equation the separation needs.

**The prediction, which is E3 and is now sharp.** Separation quality should improve with the
*spread* of `n` and `u` across studies covering a voxel, and should be near-absent when the
roster is homogeneous. Specifically: fit the same true `(pi, mu)` under a roster with `n` fixed at
30, then with `n` drawn over 15-120, then with the reporting threshold varying too. If the theory
holds, the fitted `pi` should be badly biased and nearly flat in the truth under the fixed roster
and should track the truth increasingly well as the spread grows -- while `pi*mu` stays about as
good in all three. That is a clean separation of "the product is identified" from "the factors
are", and it is decisive either way.

**What it would give the user if it holds.** A diagnostic the estimator can compute before
fitting: the spread of sample sizes among studies contributing at each voxel. Where that spread
is small, report the product and refuse the split. That is a far better guard than a warning
about "regimes", because it is computed from the collection in hand.

**And a design direction.** If identification comes from heterogeneous power, then the estimator
should be *told* each study's power rather than inferring it from a reported height -- which
means the sample size and the threshold are the load-bearing inputs, and the statistic column is
close to decoration. That is testable too: refit with the reported magnitudes replaced by their
per-study mean, destroying all within-study height information, and see how much is lost.

### 18a. Precision on the last claim: weakly identified, not unidentified

Section 18 says a homogeneous roster leaves only the product identified. That overstates it, and
the overstatement matters because it predicts the wrong result for the test.

There is a second identification route that survives a homogeneous roster: the *shape of the
reported height distribution above the common threshold*. Reported values are a truncated sample,
and a truncated normal's shape depends on how far its mean sits below the cut, so the observed
heights do carry information about `mu` independently of how many studies reported. That is
precisely the channel the censored likelihood is built on.

It is a weak channel, not an absent one. Two measurements bound it: reported height moves 0.055 z
per unit of true g at a 3.29 cut, and a truncated-normal MLE fitted to genuine local maxima
returns 0.257 for a true 0.5 and 3.6 for a true 2.0 — biased in both directions, which is why
that model was tried and rejected here. So the residual route is both weak and misspecified for
local maxima.

The corrected claim, and the one the test should be read against:

> With a homogeneous roster the prevalence/magnitude split is **weakly** identified, through the
> shape of the truncated height distribution alone — a channel carrying a few percent of the
> count channel's signal and misspecified for maxima. Heterogeneous power adds a second,
> stronger channel. The product is well identified throughout.

Predicted shapes, so the result can falsify rather than accommodate: under a fixed roster the
fitted prevalence should have a small but non-zero slope on the truth — poor, not flat. The slope
should rise as `n` and then `u` are allowed to vary. The marginal's slope should be respectable in
all three and change least. A *flat* fixed-roster slope would mean the height channel contributes
nothing at all, which is stronger than anything measured so far and would be worth knowing; a
fixed-roster slope already near 1 would kill the theory outright.

## 19. An estimand that is identified, interpretable, and in the strong channel

The audit has been negative for a while, so it is worth asking what a coordinate table *can*
support, designed forward rather than patched.

From section 18, the quantity the data directly measure is

    P_k(v) = pi_v * D(mu_v, n_k, u_k)
           = the probability that study k reports a focus within r of voxel v.

This is the *reporting probability*, and it has properties none of the current outputs have.

**It is identified.** It is the expectation of an observed Bernoulli. No scale constant, no image
donor, no separation of prevalence from magnitude. Everything measured in this program that went
wrong went wrong in factoring `P` into `pi` and `D`; `P` itself is what the counting gives you.

**It is in the strong channel.** It is built from whether and where foci appear -- the channel
that carries roughly ten times what the heights do (section 17).

**It is interpretable without a convention.** "In a study of 30 subjects reporting at p < 0.001
with cluster-extent correction, the probability of a focus within 10 mm of this voxel is 0.42."
That is a sentence a reader can check against their own experience and a reviewer can argue with.
Compare `g = 1.31` on a scale the coordinates cannot identify, or `prevalence = 0.47` that must be
read ordinally and is right in half of maps.

**It answers the question that was actually asked.** "How can I use this to power my next study?"
has no good answer in terms of `g` or `pi`: power under zero inflation needs both factors
separately, and both are badly estimated. It has a direct answer in terms of `P`: standardise to
the planned design and read off the chance of reporting a focus there. And unlike a power
calculation from `g`, it does not require the magnitude scale to be right -- it requires only that
the detection function be interpolated over the range of `n` and `u` the collection contains.

**Standardisation is the whole trick, and its limits are honest ones.** Fit `D` with `n_k` and
`u_k` as covariates, then evaluate at a stated reference `(n0, u0)`. The map is then "reporting
probability for a reference study", and the extrapolation is legitimate only across the range of
power the collection actually spans -- which is the same heterogeneity that identifies the fit.
A collection of twenty 20-subject studies at one threshold can report `P` at that design and
should refuse to extrapolate to N = 100. That is a real restriction, and it is *stateable*, unlike
"read it ordinally".

**How it differs from ALE and MKDA.** Those produce a convergence statistic whose units are the
kernel's, referred to a null. They answer "is there more agreement here than chance", which is a
hypothesis test. `P` answers "how often would a study like mine find this", which is an estimate
with a scale. The two are complementary, and the second is the one a reader wants when the answer
to the first is yes. That is also the honest version of the claim the PR currently makes for `g`:
an effect-size-like quantity on coordinates -- except this one is identified.

**What it gives up.** It is not an effect size. It will not combine with an IBMA, it does not
answer "how big is the effect", and it depends on the reporting conventions of the literature it
was fitted to -- if the field's thresholds shift, `P` at a fixed reference shifts with them. Those
are genuine losses. They are smaller than the loss of reporting a magnitude that is 40% high with
an interval that covers 10% of the time.

**Cheap test of the proposal, before believing any of it.** Simulate collections with known `pi`
and `mu` and a spread of `(n, u)`; fit a logistic or complementary-log-log detection model for
`P_k(v)` with `n` and `u` as covariates; evaluate at a reference design; compare against the
simulated truth `P` at that design, computed exactly. Score calibration (does a predicted 0.4
happen 40% of the time) and the standardisation error as the reference moves away from the
collection's centre of mass. If that calibrates while `g` does not, the proposal is worth putting
to the maintainer as an additional output -- not a replacement, since the null and the map
machinery are already built and would carry it.

## 20. There is a window of detectability, and both factors are identified inside it

Section 18 said the prevalence/magnitude split is identified by the spread of study power, and
18a weakened that to "weakly identified". Both were too coarse. The measurement
(`reporting_probability.py`: 200 collections, 24 studies each, `n ~ U(20,60)`, four thresholds,
and an occupancy likelihood fitted on *exact* detection records with the *true* detection
function -- so a best case, a ceiling on what CBES could do) says something sharper.

```
        true mu     0.00    0.20    0.40    0.60    0.80
      median mu   -0.325   0.208   0.425   0.611   0.769
        mean mu   -0.103   0.247   0.544   0.758   3.306
   fits saturated    0.01    0.04    0.06    0.09    0.26
        true pi     1.00    1.00    0.75    0.75    0.50
      median pi    0.000   0.133   0.687   0.767   0.542
   dD/dmu at truth   0.000   0.341   1.876   2.195   0.546
```

**Both factors are recovered well where the detection gradient is large.** At the three middle
sites the medians are 0.425/0.611/0.769 against true 0.40/0.60/0.80 and 0.687/0.767/0.542 against
true 0.75/0.75/0.50. That is not weak identification; that is a working estimator. And the
failures line up exactly with `dD/dmu`, the slope of the detection probability in the magnitude:
1.876 and 2.195 where it works, 0.000 / 0.341 / 0.546 where it does not.

**The failures are at the two ends, for two different reasons.**

- *Below the window*: at `mu = 0` there is nothing to detect and `pi` is genuinely unidentifiable
  -- median `pi` comes back 0.000 for a true 1.00, which is the correct answer to an
  ill-posed question, not an error. At `mu = 0.2` the gradient is 0.341 and `pi` reads 0.133.
- *Above the window*: at `mu = 0.8` with `n ~ 40` the effect clears the cut almost always, so
  detection saturates, `mu` is unidentified from above and runs away -- 26% of fits have it
  pinned -- while `pi` absorbs the level. That is why the *mean* `mu` is 3.306 while the median
  is 0.769.

So: **prevalence and magnitude are separately identified in a window of detectability, roughly
where `mu*sqrt(n)` is within about 1.5 of the reporting threshold.** A spread of sample sizes and
thresholds widens the window, because different studies place the same `mu` at different points
on their own detection curves -- which is the grain of truth in section 18, now with a
mechanism and a location.

**This explains CBES's prevalence behaviour, quantitatively and in one piece.** Every pathology
recorded separately is a position relative to the window:

| observed | position |
| --- | --- |
| floor near 0.2 at a true zero, and non-monotone at weak effects | below the window: no detection gradient, so the fit is driven by the likelihood's shape |
| `prevalence` 0.994 at a site with peak z = 4.4 in the coverage bed | above the window: detection saturated, prevalence absorbs the level |
| `prevalence` 0.68 for a true 1.0 on held-out HCP | inside or near the window, where it is merely compressed |
| ordering right on average but exactly right in only half of maps | voxels within one map sit at different points, so they are not comparably identified |

That last row is the important one for the docstring. A map spans magnitudes, so it spans the
window -- some voxels are below it, some above, some inside. Comparing prevalence across voxels
compares quantities that are identified to different degrees, which is why the ordering is
reliable on average and unreliable in any one map. The fix is not calibration; it is a
per-voxel statement of whether the voxel is in the window at all.

**Correcting myself twice over.** I first read the mean `mu` of 3.306 against a true 0.8 as a
likelihood ridge and checked whether `(0.80, 0.50)` and `(3.31, 0.585)` make the same predictions.
They do not -- they differ by up to 0.28 over the design range, so the ridge story was wrong. The
mean was simply a bad summary of a distribution with a 26% saturated tail. Two lessons, both
mine: a heavy-tailed parameter should have been summarised by its median from the start, and
"weakly identified" was a vague label that stopped me looking for the structure, which turned out
to be a clean one.

**And the proposal in section 19 survives its first test.** Reporting probability standardised to
a reference design, estimated/true:

```
  ref N  in range   mu=0.2      mu=0.4      mu=0.6      mu=0.8
     15        NO   0.01/0.01   0.05/0.03   0.17/0.13   0.25/0.21
     30       yes   0.02/0.01   0.12/0.10   0.39/0.37   0.41/0.43
     45       yes   0.03/0.03   0.20/0.20   0.56/0.58   0.50/0.49
     60       yes   0.05/0.04   0.28/0.32   0.66/0.68   0.54/0.50
    100        NO   0.10/0.10   0.46/0.57   0.75/0.75   0.58/0.50
    200        NO   0.22/0.32   0.63/0.74   0.77/0.75   0.58/0.50
```

In range the largest error is 0.04. Out of range it degrades to 0.11 and does so smoothly and in
a consistent direction (under-estimating at larger N, because the saturated fits cap at their own
`pi`). That is the behaviour the proposal needed: right where the collection supports it, visibly
and disclosably wrong outside. Note that the *product* is well estimated at `mu = 0.8` (0.41/0.43,
0.50/0.49, 0.54/0.50) at exactly the site where the factors are worst -- which is the whole point.

## 21. The z-to-g conversion rests on a degrees-of-freedom assumption nobody reports

This came out of my own convention bug, but it is a property of the estimator and it is worth
having on its own. When a study reports a z-map peak, CBES maps z back to a t on `n - 1` degrees
of freedom -- treating the reported z as a p-value-preserving image of a t, which is what most
neuroimaging software produces -- and then converts the t to an effect size. Right in principle;
the difficulty is that reported peaks sit far into the tail, where that map is steep.

Effect size recovered from a reported z at `n = 30`, varying only the assumed residual df:

```
reported z   df=29   df=60   df=120   df=1000   spread
      3.30   0.653   0.626    0.614     0.604    1.08x
      4.00   0.830   0.776    0.752     0.733    1.13x
      5.00   1.133   1.009    0.959     0.918    1.23x
      6.00   1.522   1.273    1.178     1.105    1.38x
```

**The assumption is load-bearing and the sensitivity grows with the reported height** -- so it is
worst exactly where the winner's curse is worst, and in the same direction. A z of 6 read with
df = 29 when the map's effective df was really in the hundreds over-states the effect size by
38%.

The effective df of a published z-map is often not `n - 1`. Variance smoothing raises it, and
FSL's FLAME does that deliberately; a mixed-effects analysis has its own; software differs in
what it writes into a z-map; and papers seldom state it. So a collection is not merely uncertain
about this -- it is systematically likely to have a higher effective df than the sample size
implies, which biases the magnitude upward.

Two modest responses, neither requiring new modelling: say so in the docstring next to the
existing magnitude caveats, and allow an explicit per-study df in the metadata so a caller who
knows it is not forced to let the sample size stand in. A study reporting a *t* is unaffected --
the conversion from t is direct -- which is a small argument for preferring `stat_column="t_stat"`
where a collection offers both.
