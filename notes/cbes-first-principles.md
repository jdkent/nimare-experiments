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
