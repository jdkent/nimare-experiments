# The same problem in other fields

Stated without neuroimaging vocabulary, the problem is: **recover a latent continuous field from
thresholded, selectively reported point patterns whose reporting threshold varies by observer and
is not reliably known, given that a few observers report everything.** Several fields have that
problem and have been working on it for decades. Three of them change the design.

---

## 1. Ecology: integrated species distribution models

The closest structural match, and it is close enough to be embarrassing.

Presence-only records -- museum specimens, citizen-science sightings -- are abundant, spatially
biased, and carry no information about where nobody looked. Planned presence-absence surveys are
scarce, unbiased, and know their own effort. An **integrated SDM** fits both to one latent
intensity surface: the presence-only records as a *thinned* Poisson process, where the thinning
function is the observation bias, and the survey data with a binomial likelihood, in a single
joint likelihood.

Map it across: presence-only records are coordinate studies, planned surveys are image studies,
the thinning function is the reporting threshold, the latent intensity is the effect field.

**What transfers, and it is the important one.** Presence-only data alone *cannot* separate the
intensity from the sampling effort -- this is a known structural non-identifiability, not a
small-sample problem. The planned surveys are what break it. That is a far better justification
for the two-image floor than the one in the requirements document, which argued from wanting an
interval on a calibration constant. The real reason is that without unbiased data the effect
field and the reporting propensity are not separately identifiable at all, and every attempt to
recover the magnitude from coordinates alone was fighting that rather than a modelling error.

**Second thing that transfers.** Fithian, Elith, Hastie and Keith (*Methods in Ecology and
Evolution* 2015, <https://doi.org/10.1111/2041-210X.12242>) assume the sampling bias is *shared
across species* and borrow strength across them to estimate it. The analogue is a reporting
function shared across studies or across collections, estimated inside the joint likelihood.
This project already tried that idea and abandoned it as task 32, "does one f(threshold) correct
g across collections" -- but as a post-hoc rescaling of a fitted map, which is why it failed.
Inside a joint likelihood it is a different and much better-posed proposition.

Mature software exists: `intSDM` and `PointedSDMs` in R, over INLA.

## 2. Seismology: catalogue completeness

An earthquake catalogue is a list of events above a detection threshold that varies with the
network, the time of day, and the ongoing aftershock sequence. Estimating the Gutenberg-Richter
law from it requires handling that.

Ogata and Katsura (1993) do not estimate the threshold and then discard everything below it.
They write the observed magnitude-frequency distribution as *the product of the true exponential
law and a smooth detection probability*, and estimate both jointly -- detection as a cumulative
normal or logistic in magnitude, with parameters `mu` (50% detection) and `sigma` (the width of
the partial-detection band). Completeness is then a derived quantity, `m_c = mu + n*sigma`, not
an input.

**What transfers.** This estimator currently infers each study's threshold from its smallest
reported value, which is a hard cut *and* a function of the study's signal -- the floating-cut
problem that has invalidated several measurements here. Replacing it with a per-study **soft
detection function estimated jointly with the effect field** fixes both. It also absorbs the
thing that produced the 0.82-to-2.29 convention swing: a smooth detection curve does not need to
know whether the authors used FDR, voxelwise family-wise error, or a cluster-extent test, because
all three produce a monotone increase in reporting probability with the statistic and the curve
just fits it. A hard threshold cannot represent cluster-extent reporting at all.

## 3. Econometrics: inference on winners

Andrews, Kitagawa and McCloskey (*Quarterly Journal of Economics* 2024,
<https://doi.org/10.1093/qje/qjad043>) give median-unbiased estimators and confidence intervals
that are valid *conditional on which target was selected*, for exactly the structure here: you
report a quantity because its estimate was the largest, so the estimate is biased and the usual
interval does not cover.

**What transfers, mostly as a warning.** Their central practical finding is that purely
conditional procedures, while valid, can be extremely imprecise, which is why they propose
*hybrid* procedures combining conditional and projection intervals. That retrospectively explains
the truncated-normal experiment recorded in `cbes-measurements.md`: fitting a conditional
likelihood to selected peaks returned 0.26 for a true 0.5 and 3.6 for a true 2.0. That is not a
coding failure, it is the known behaviour of conditional inference when the conditioning event
carries most of the information. The hybrid construction is the principled fix if conditional
inference on peak values is wanted at all.

## 4. Statistical genetics: the winner's curse toolbox

Genome-wide association studies have the identical problem -- effect sizes for the SNPs that
passed a genome-wide threshold are exaggerated -- and, like meta-analysis, usually have only
summary statistics rather than individual data.

**What transfers, and it is a redirection.** The systematic review by Forde and colleagues
(*PLOS Genetics* 2023, <https://doi.org/10.1371/journal.pgen.1010546>) evaluates the whole
toolbox and reports that **conditional-likelihood methods perform poorly**, while bootstrap and
empirical-Bayes shrinkage estimators are competitive. That is a pre-existing, independent verdict
against the approach tried here, and a specific pointer to what to try instead. There is an R
package, `winnerscurse`, implementing the alternatives on summary statistics alone.

## 5. Astronomy: Eddington bias and sub-threshold statistics

When source counts rise steeply toward faint fluxes, noise scatters more sources up across the
detection limit than down, so the observed bright counts are inflated -- Eddington bias, named in
1913. Submillimetre astronomy corrects it by forward-modelling the count distribution with a
population prior and deconvolving, iteratively.

**What transfers.** Two things. First, the correction is done by forward-modelling the *counts*
with an assumed population distribution, not by debiasing individual measurements -- the same
move as putting the effect into the intensity rather than the marks. Second, "astronomy below the
survey threshold" uses the statistics of the *undetected* population, fitting the pixel intensity
histogram to constrain counts below the detection limit rather than discarding that region. The
analogue is using the absence of reported foci quantitatively, which the censoring term attempts
but at the level of individual values rather than of the count distribution.

Extreme value theory's peaks-over-threshold framework is the formal home for all of this, and it
handles non-stationary and covariate-dependent thresholds directly.

---

## What this collectively says to do

1. **Justify the two-image floor by identifiability, not by wanting an interval.** Ecology says
   the coordinate-only problem is structurally non-identifiable. That is a stronger statement
   than anything in the current documentation and it explains the whole sequence of failures.
2. **Replace the hard, inferred threshold with a jointly estimated soft detection function per
   study.** Ogata and Katsura, 1993. This is the single most actionable import and it targets the
   convention-dependence directly.
3. **Estimate the reporting function by borrowing across studies inside the joint likelihood**,
   rather than fitting a correction to an already-fitted map. Fithian et al.
4. **Do not return to conditional-likelihood debiasing of peak values.** Two independent
   literatures say it underperforms, and a measurement here already showed it failing.
5. **If peak values are debiased at all, use bootstrap or empirical-Bayes shrinkage**, and for
   the image studies specifically, Davenport and Nichols' resampling correction, which needs the
   full image and is therefore exactly suited to the two-plus images this scope guarantees.
