# CBES: what was measured, and where

Findings that used to live in `nimare/meta/cbma/effectsize.py` docstrings and comments. The
code keeps the decision; this keeps the evidence. Script names are in `experiments/`.

## Magnitude is not recoverable from coordinates

`peak_information` on the 21 NIDM pain studies: reported peaks average z = 3.639 against a
null-peak expectation of 3.63 at the same threshold, an excess of **+0.009**, where the effect
actually present at those locations would have produced +1.10. The peak heights carry
essentially no effect-size information. This is a fact about underpowered fMRI reporting, not
about the estimator, and it means no method built on peak heights — CBES, SDM-PSI, or a
selection model exploiting the threshold pile-up — can recover absolute magnitude.
(`measure_bias.py`, `signal_vs_noise.py`)

Against the pain images CBES reaches rho = 0.84 (ALE 0.18) and overestimates about twofold. On
an 11-study NeuroVault collection of a weak contrast — 35 peaks reaching 11% of the brain — it
overestimates by **10.7x** and the rank correlation falls to 0.13. Twofold is the friendly case.
(`second_collection.py`, `nv_thresh.py`)

Both validations are partly circular: the coordinates are manufactured by thresholding the same
images used as truth, with the reporting model CBES assumes.

## SDM-PSI comparison (`compare_sdm.py`)

Same peaks, same 21 images as truth, image mean = 0.493.

| method | r all | r top | mag ratio | mean abs est |
|---|---|---|---|---|
| SDM-PSI | 0.800 | 0.730 | 0.47 | 0.242 |
| CBES, known threshold | 0.775 | 0.814 | 1.99 | 0.967 |
| CBES + rho | 0.787 | 0.809 | 2.06 | 0.975 |

Both miss by ~2x in opposite directions — SDM deflates (its imputation bounds unreported voxels
below significance), CBES inflates. SDM emits absolute Hedges' g without saying the scale is
unearned. CBES wins on the top quartile, which is the identified claim.

## Image calibration of the scale (`images_x_peaks.py`, `which_image.py`)

Donors held out of the truth so no image is scored against itself. Pooling the images into one
fit made the scale drift with donor count (0.610, 0.590, 0.516 for 1/2/5). One ratio per donor,
pooled, holds it: 0.610, 0.615, 0.622. Thinned to 20 peaks/study it still slides 0.57, 0.53,
0.51 — hence `scale_interval_`.

`peak_bias_scale="reference"` borrows from a NeuroVault corpus matched on sample size, and made
recovery of known truth 1.5-2x worse when the collection was not typical of that corpus.
Sample size predicts magnitude (corr(log N, log magnitude) = -0.395, 28% error reduction);
*spatial* similarity does not (0.26% of variance). (`rate_vs_image.py`, `similarity_test.py`)

## Single-study voxels (`one_study_voxel.py`, `profile_cost.py`, `gate_check.py`, `gate_sweep.py`)

The likelihood there is a **plateau, not a ridge**: with one reported peak and five or more
silent studies, every mu in [0, 2] is within 2 log-likelihood units of the profile maximum,
because one report cannot separate a moderate effect from a false positive at mu = 0. Two
reporting studies tighten it to [0.16, 0.95]. The specific bound is configuration-dependent —
varying sigma and cutoff across studies gave [0.05, 1.43] — but the conclusion is not.

So the EM never converges there: g moved up to 0.24 between max_iter 25 and 200, and 0.06
between 400 and 800. Stopping on the likelihood makes the value reproducible and nothing more.

How common: a function of the peak count, not the collection.

| peaks/study | brain covered | k=1 of covered | median k |
|---|---|---|---|
| 3 | 19.5% | 70.1% | 1 |
| 10 | 51.7% | 57.5% | 1 |
| 20 | 77.2% | 39.9% | 2 |
| all (~130) | 99.9% | 0.9% | 7 |

A profile MLE would cost 9-13 mu evaluations against the EM's 25, but the interval endpoints
add two root-finds for 31-35 total — so it buys reproducibility, *not* speed. And se does not
protect z at those voxels: median se 0.81 at k=1 vs 0.58 at k>=5, and 13-19% still clear
|z| = 2. Validity there rests on the permutation null, which is valid for any statistic.

## Inference (`permutation_null_validate.py`, `few_images_null.py`, `selection_none_fpr.py`)

Uncorrected rate 0.022-0.057 against nominal .05 over 20 simulations at 200 permutations, across
global nulls, foci confined to a quarter of the mask, and 1-5 image studies. Power at focal
g = 0.8 across 30 studies: 18/20 at voxel FWE. Conservative throughout, lowest where an image
study's sign flip gives the null only a few states.

The relocation null was removed: it is not conditional on multiplicity, so a convergent site hit
the p floor under both nulls. Permuting whole rows also broke the multiplicity invariant —
relabelling put two foci of one study on a voxel and dropped the null's study count from 30 to
~20. Only the value columns are permuted. (`real_geometry_null.py`, `relocation_domain.py`)

## tau2 and the standard error (`tau2_hksj.py`)

Kernel-weighted Paule-Mandel and REML are indistinguishable from DL at these study counts
(bias at true tau2 = 0: DL +0.030/+0.023/+0.016 at k = 3/5/10; PM +0.031/+0.024/+0.017; REML
+0.028/+0.022/+0.015). The textbook downward-bias concern does not reproduce — zero-truncation
induces an offsetting upward bias. No change warranted.

HKSJ is a clear win on coverage of nominal-95% intervals:

| k | tau2 | model SE | HKSJ |
|---|---|---|---|
| 3 | 0.20 | 86.4% | 95.9% |
| 5 | 0.20 | 89.2% | 95.4% |
| 10 | 0.20 | 92.0% | 95.3% |

Its df must be Kish's n_eff, not sum(w): only Kish's form is invariant to rescaling the weights,
and a voxel reached by distant foci has weights summing well under one, which sends the df to
zero (read as 99.3% coverage before the fix). Same trap broke my PM target `sum(w) - 1` and made
weighted REML unbounded below.

HKSJ only reaches the `selection_model="none"` path: the zero-inflated model discards the pooled
inverse-variance SE and reports the censored likelihood's curvature instead.

## coverage_radius (`radius_threshold.py`)

Swept 8-34 mm. The effect-size map barely notices at realistic peak counts (at 10 peaks/study,
r_top 0.108-0.109, magnitude 3.74-3.78). On dense tables a wider radius trades magnitude for
pattern: r_top 0.39 to 0.51, magnitude 1.70 to 2.51.

`prevalence` is what depends on it: 0.165 to 0.438 at 20 peaks/study, and on dense tables it
**saturates** (0.997 at the 20 mm default, exactly 1.0 by 26 mm). That is a dense-table
artefact — at realistic counts prevalence sits at 0.08-0.21 — but prevalence figures measured on
peaks thresholded out of whole images should not be quoted.

## Performance (`profile_cbes.py`, `kernel_share_big.py`, `f32_censor.py`, `fused_proto.py`)

Profile of a 63-study whole-brain fit: censoring kernel 21%, E step 19%, coverage sets 11%,
derivatives 9%. An earlier docstring claimed the kernel was 63%; it is not. Everything is
memory-bound, so Amdahl caps any single-kernel win around 1.25x.

Rejected with measurements:
- **float32 censoring kernel**: 1.43x on the kernel at 8M pairs, but a fit holds only ~0.5M
  pairs per call whatever the study count, giving 1.03-1.06x overall. Not worth single precision
  in the likelihood. (One pair of fits read 1.40x — that was warm-up.)
- **Fused numba EM sweep**: agrees with numpy to 1e-13, but only 1.36x serially; the 3.75x came
  from eight threads, and the permutation null already parallelises whole fits.
- **SQUAREM**: made it worse, 26 to 76 censoring evaluations.
- **Empirical-Bayes Beta shrinkage of pi**: a prior centred on the collection mean pulled
  correct estimates (0.490, true 0.5) down to 0.315; shrinking toward zero compounded across
  sweeps (median 0.397 to 0.077 at k=1). Discarded on real-data evidence.

Kept: padded-grid dilation in `_coverage_entries` (2.4-3.2x, entry-for-entry identical),
geometry caching, likelihood-based EM stopping (3% at max_iter 25, 15% at 200).

`selection_model="none"` is **9.9x** faster (23.1s to 2.3s), not the "3x" once documented —
without silent pairs there are no coverage sets, no E step, and no plateau, so the EM converges.
The whole 10x is the selection correction.

## Threshold inference

`infer_threshold_from_minimum` assumes height thresholding with one local max per cluster.
Cluster-extent thresholding lifts the inferred threshold by 0.19-0.57 z, and most modern fMRI
papers use cluster-extent or TFCE. `threshold="pooled-min"` applies the most liberal study's cut
to everyone. (`thresh_est.py`, `thresh_check.py`, `radius_threshold.py`)

## The within-analysis null (2026-09-14/15)

Shuffling reported values across the whole focus table let a large-N study's peak land on a
small-N study's voxel. Restricting the shuffle to within an analysis fixed it: under a global
null with sample sizes from 10 to 1000, the uncorrected rate went from 96.7% to 0.0125, and it
no longer depends on how widely precision varies (0.0113 at N 20-40). Images take the same
within-study rearrangement rather than a sign flip, so both sides randomize one hypothesis.
(`within_analysis_null_rates.py`)

The cost is power. At a focal g = 0.8 across 30 studies: 0.710 uncorrected, **0.030 voxel FWE**.
The observed max |z| sits only ~19% above the null median, because permuting within a study
keeps every large value in the map and merely scatters it, while a maximum over thousands of
voxels catches whichever scattered voxel drew well. Cluster size and mass are statistics of
concentration and should survive better; under the global null they read 0.000 and 0.050 against
a voxel rate of 0.100. (`cluster_power.py`)

**A measurement error worth remembering.** `within_analysis_null_rates.py` first reported the
mean *share* of voxels rejected, not the familywise rate. That reads 0.0000 where the rate is
about 0.10, and it made the test look far more conservative than it is.

## Finding 7 and the magnitude (2026-09-14/15)

Four instruments in a row were wrong in ways that each produced a confident conclusion:

1. The likelihood already integrates to one over reported-and-silent, so it is a censored
   likelihood and adding a truncation normaliser double-charges the selection. Prevalence
   collapsed to its floor and z fell to 0.026.
2. The default point simulator has no peak-height inflation at all, so it cannot test a
   peak-height correction.
3. The field simulator wrote `signal * scale + noise` as a Z when that is the *t* statistic's
   noncentrality. Converting its own peaks returned g = 8.8 against a true 1.6, and it inverted
   the apparent direction of the bias. Fixed with `t_to_z`; the natural reporting density also
   fell from 28-29 peaks/study to 13.8, inside the realistic range.
4. Every real-data reference estimated the wrong quantity. The inverse-variance pooled image
   map estimates pi*mu; CBES estimates mu. `g * prevalence` matches it (ratio 1.17, r 0.575)
   where `g` alone reads 2.3-9.9x high.

Corrected simulator, prevalence 1: ratio to truth 0.23 / 0.68 / 0.97 / 0.93 / 0.95 / 0.95 at
true g = 0.2 / 0.4 / 0.6 / 0.8 / 1.2 / 1.6. Above ~0.6 the required scale is 1.05 +- 0.02 across
a 2.7-fold range, which is not what an unidentified constant looks like.

`peak_information` does not discriminate: excess z is +0.03 at true g of 0.2, 0.4 *and* 0.6,
where the ratio runs 0.23 to 0.97. Real collections sit in that same band, so excess z near zero
says nothing about whether the magnitude is right.

The conditional reference built by conditioning a study on clearing its own threshold is
inflated and nearly flat in true signal: it reads 1.09-1.42 while the pooled truth varies
tenfold. CBES moves 0.42-1.13 over the same range, so it is the more responsive instrument.
(`matched_by_signal.py`)

Across five collections and four thresholds the CBES/reference ratio rises monotonically with
the inferred cut in every one, `ratio = 0.328 + 0.203*cut` (r 0.651, residual sd 0.126), and
`log(peaks)` does not predict it (r -0.238, p 0.34). But the roster shrinks as the threshold
rises, so the within-collection slope is confounded with study count.
(`threshold_correction_fit.py`)

**What none of this could settle**, and why the HCP work exists: with one map per study, "this
study has an effect here" is not observable independently of the effect's size. Collection 4337
has per-subject maps for 787 subjects across 23 contrasts, so synthetic studies can be built
from one set of subjects and the truth measured on a disjoint set.
(`hcp_index.py`, `hcp_heldout_validation.py`)

## The magnitude map reads the reporting threshold, not the effect (four corpora)

`peak_selection_event.py` settles what the abandoned truncation patch got wrong. Conditioning a
reported value on `|g| > c` is a *single-draw* selection event, and the values in a coordinate
table are not single draws -- they are local maxima of a smooth field, which clear a cut because
they are the largest thing in their neighbourhood. Applying the truncated-normal correction to
peaks drawn from a field with a true mean of 0.5 returns 0.257; at a true mean of 2.0 it returns
3.6. The 0.259 the simulator fixture produced under that patch was the model speaking, not a bug,
so the patch was reverted rather than tuned.

`peak_height_curve.py` measures what is left to work with. Mean `|peak|` against the true field
mean, smoothed white noise at unit variance:

| cut | floor at true mean 0 | slope 0 to 1 | slope 1 to 3 |
| --- | --- | --- | --- |
| 2.50 | 3.09 | 0.208 | 0.752 |
| 3.29 | 3.84 | 0.055 | 0.481 |
| 4.00 | 4.53 | 0.028 | 0.249 |
| 5.00 | 5.52 | 0.044 | 0.052 |

At a cut of 3.29 a whole unit of true effect moves the reported height by 0.055. The height is
the cut plus an overshoot that has forgotten the signal, which is what the estimator's own
`peak_bias` warning has been saying.

`cross_dataset_floor.py` asks whether this shows up outside the held-out HCP design, which is
the only question that decides anything. Studies are split in half: one half becomes thresholded
peak coordinates for CBES, the other half is pooled by inverse variance for the truth, so the
reference is never conditioned on the noise that produced the peaks.

| corpus | what the studies share | r(g, truth) | r(pi*g, truth) | ratio in the top stratum |
| --- | --- | --- | --- | --- |
| HCP held-out subjects | one task, one population (pi = 1) | +0.50 to +0.69 | -- | 1.37 to 1.57 |
| NIDM pain, 21 studies | one construct, different labs | +0.14 to +0.22 | +0.34 to +0.39 | 1.9 for g, 1.2 for pi*g |
| passive viewing, 12 collections | a cogatlas label | -0.078 | +0.057 | 3.7 |
| go/no-go, 12 collections | a cogatlas label | +0.097 | +0.028 | 3.1 |
| emotional regulation, 10 collections | a cogatlas label | -0.098 | +0.014 | 3.0 |

On pain the truth spans 11x across its strata (0.084 to 0.961) while `g` spans 1.2x (1.50 to
1.82); on two of the three NeuroVault groups `g` *falls* as the truth rises. Subtracting the
floor does not repair it: `g - u/sqrt(N)` is +0.78 against a truth of 0.41 on pain and +0.75
against 0.07 on passive viewing.

Read the NeuroVault rows with their caveat. Those groups share a paradigm *label*, not a
contrast, so the held-out reference is itself nearly flat -- its top stratum is only 0.24 to
0.35 -- and a near-zero correlation there is partly an absent common effect rather than purely
the estimator's failure. What they do establish is that no magnitude claim can be supported on
that kind of corpus.

The ordering across the three is the informative part: the magnitude map degrades exactly as
spatial agreement between the studies degrades. That is expected of any coordinate-based
estimator. What does not survive is the absolute scale, in every corpus including the one where
prevalence is 1 by construction.

`pi * mu` -- the `g_marginal` map that was removed -- is uniformly the better estimate of the
held-out truth on pain: correlation +0.34 against +0.22 for `g`, and a top-stratum ratio of 1.21
against 1.89. It does not rescue the NeuroVault groups, where nothing does.

## Redone with realistic reporting: the level moves with the convention, the shape never moves

The measurements above extracted peaks at a fixed uncorrected cut and kept the strongest ten per
study. Both are wrong, and the second is wrong in a way that had already invalidated two earlier
results: a cap fixes the count and lets the effective threshold float to the study's tenth
strongest peak, which is a function of its signal. `reporting.py` replaces it -- corrected
thresholds, only surviving clusters, one row per cluster, nothing capped -- and every corpus was
re-run through it.

NIDM pain, 21 studies, split half so the reference comes from studies the coordinates did not:

| scheme | focus | foci/study | ratio of g | top stratum | r(g) | r(pi*g) |
| --- | --- | --- | --- | --- | --- | --- |
| FDR q=0.05 | max | 108 | 3.23 | 1.49 | +0.331 | +0.194 |
| FDR q=0.05 | centre of mass | 112 | 2.52 | 0.82 | +0.236 | +0.146 |
| voxelwise FWE | max | 13 | 3.53 | 1.81 | +0.238 | +0.095 |
| voxelwise FWE | centre of mass | 14 | 3.19 | 1.61 | +0.164 | +0.227 |
| cluster extent | max | 6 | 4.58 | 2.29 | +0.277 | +0.483 |
| cluster extent | centre of mass | 6 | 2.49 | 1.15 | +0.135 | +0.320 |

Reporting the coordinates properly is worth a lot: r for `g` under FDR rises from +0.217 under
the capped extraction to +0.331, so the earlier protocol was penalising the estimator rather
than exposing it.

Two things follow, and they point in opposite directions.

**The level is a property of the reporting convention, not of the effect.** Across the six rows
the truth in the top stratum is essentially constant -- 0.92 to 1.02 -- while `g` there reads
0.82 to 2.29. Cluster extent with a centre-of-mass focus, which is both the most common practice
and the one that produces a plausible six-row table, happens to land at 1.15. That is the
closest to unbiased anything has reached against an independent reference, and it is a
coincidence of convention: change the convention and the same data give 2.29.

**The shape never improves.** Under cluster extent with centre of mass the held-out truth spans
0.083 to 0.935 across its strata, an eleven-fold range, and `g` spans 0.894 to 1.076, a
1.2-fold one. Every other row is the same. No reporting scheme, no focus convention, no
selection model and no `peak_bias` setting has moved that compression.

A centre-of-mass focus is not selected for being extreme, so it should carry no winner's curse,
and it does measurably reduce the level -- 4.58 to 2.49 under cluster extent. It does not touch
the compression. Whatever flattens the map is therefore not peak selection.

### What the realistic pipeline does to the model's assumptions

It repairs one and breaks two.

*Repaired.* One focus per cluster makes a coordinate a sparse region marker, which is what the
20 mm coverage design assumes it is, instead of one of a hundred near-duplicate maxima.

*Broken.* The censoring term's selection event is `|g| >= c` at a height threshold. Under cluster
extent the event is "this cluster was large enough", and the cluster-forming cut of 3.09 is
nowhere near the smallest value that actually appears. Told 3.09, the model believes small values
were reportable and under-corrects; left to infer the cut from the smallest reported value, that
cut floats with the signal -- the capping problem again, produced by nature this time.

*Broken, and the likely mechanism for the compression.* A large cluster yields one focus, and
everything in it beyond the coverage radius is entered as the study having been *silent* there.
A bigger true effect makes a bigger cluster and so a larger misread fraction, which pushes the
estimate down hardest exactly where the effect is largest. `silence_misread.py` measures it.

### The silence misreading is real, large, and not the cause

`silence_misread.py` confirms the mechanism exists. Under cluster-extent reporting with one
focus per surviving cluster, 46.3% of the voxels a pain study found significant lie outside the
coverage radius of that cluster's single focus, and the estimator enters every one as the study
having been silent there. It scales with cluster size as predicted -- 54% and 61% for the two
studies with the largest clusters, 0.7% for the study whose largest cluster is 63 voxels -- and
under voxelwise family-wise thresholding, where surviving blobs are tiny, it is only 2.2%.

`cluster_extent_coverage.py` then removes the error and nothing else, by giving each study an
analysis mask covering everything except its own significant-but-uncovered territory. This is
the *oracle* version of importing cluster extent: the true significant set from the full image,
which no reported number could beat.

| truth stratum | truth g | g as now | ratio | g, silence fixed | ratio |
| --- | --- | --- | --- | --- | --- |
| 0-50% | 0.084 | 1.749 | 20.74 | 1.744 | 20.68 |
| 50-75% | 0.238 | 1.869 | 7.86 | 1.864 | 7.84 |
| 75-90% | 0.393 | 1.894 | 4.82 | 1.887 | 4.80 |
| 90-99% | 0.650 | 2.030 | 3.13 | 2.011 | 3.10 |
| 99-100% | 0.960 | 2.205 | 2.30 | 2.177 | 2.27 |

Correlation with the truth goes +0.209 to +0.194. Reclassifying 46% of significant territory
moves the estimate by about 1%, in the wrong direction for the correlation.

Why so little: a voxel deep inside a large cluster has no focus within kernel reach either, so
its estimate is set entirely by *other* studies' values. Cancelling one study's silence vote
barely moves a voxel several other studies are already speaking about.

`extent_sphere_fidelity.py` closes the line. Across 119 surviving clusters the median is 131
voxels, whose volume-matched sphere has a 13 mm radius -- smaller than the 20 mm coverage
already in use. Placed at the reported maximum it recovers 40% of the cluster against the
fixed sphere's 66%, because clusters are elongated and the maximum is not their centroid. An
extent-driven radius would cover *less* of the typical cluster than the current fixed one while
buying at most the 1% above.

The mask device works and is verified: masks load, are keyed per analysis, and the fitted values
move. The negative result is the model's, not the harness's.

### Assuming a cluster size, in the version that could have worked

Cancelling the misread silence was the weak form of the idea. The strong form lets the assumed
cluster carry the study's reported *value*, so the focus speaks for the whole cluster instead of
a 13 mm neighbourhood of its peak. `fwhm` already does exactly that, with the coverage radius
following at twice its size, and it needs no reported extent -- which matters, because papers do
not give one reliably.

| fwhm | coverage | r | 0-50% | 50-75% | 75-90% | 90-99% | 99-100% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10 mm | 20 mm | +0.209 | 20.74 | 7.86 | 4.82 | 3.13 | 2.30 |
| 15 mm | 30 mm | +0.177 | 21.23 | 7.64 | 4.79 | 3.18 | 2.28 |
| 20 mm | 40 mm | +0.150 | 21.17 | 7.51 | 4.71 | 3.15 | 2.24 |
| 25 mm | 50 mm | +0.138 | 20.50 | 7.25 | 4.56 | 3.07 | 2.19 |

The held-out truth spans elevenfold across those strata. Widening the assumed cluster from 20 mm
to 50 mm moves the top-stratum ratio from 2.30 to 2.19 and costs a third of the correlation.
Both forms of the assumption are ruled out.

### The running tally of what does not fix the compression

1. The truncated-normal selection correction. Wrong event for a local maximum: it returns 0.257
   for a true mean of 0.5 and 3.6 for a true mean of 2.0.
2. `peak_bias='per-study'`, the estimator's own documented remedy. Ratio 3.98 to 3.83.
3. Subtracting the censoring floor `u/sqrt(N)`. Fixes part of the level, none of the shape.
4. Reporting properly instead of capping. Worth a lot to the correlation, +0.217 to +0.453, and
   nothing to the shape.
5. Centre-of-mass foci, which carry no winner's curse. Level 4.58 to 2.49, shape unchanged.
6. Oracle cluster-extent coverage, using the true significant set. About 1%.
7. An assumed cluster size driving the kernel. Nothing, and the correlation falls.

Seven independent changes, no movement. That points at the input rather than the model, which
`information_ceiling.py` tests directly by comparing the reported magnitudes against the
held-out truth at the same locations, with no estimator in between.

### The ceiling is the coordinates, not the estimator

`information_ceiling.py` bypasses the fit entirely. For every focus a paper would print it takes
the effect size the table reports -- the only magnitude the estimator ever receives about that
location -- and the held-out truth at the same voxel, from studies that played no part in
producing the focus.

| scheme / focus | reported g at the focus | held-out truth there | r | regression of truth on reported |
| --- | --- | --- | --- | --- |
| cluster / max | 1.965 (sd 0.593) | 0.519 (sd 0.312) | +0.222 | 0.117 x reported + 0.289 |
| cluster / centre of mass | 1.257 (sd 0.366) | 0.556 (sd 0.302) | +0.102 | 0.084 x reported + 0.450 |
| FDR / max | 0.954 (sd 0.343) | 0.235 (sd 0.198) | +0.307 | 0.178 x reported + 0.066 |

A whole unit of reported Hedges' g buys 0.08 to 0.18 units of real effect, and the intercept
carries most of the value. One tabulated coordinate explains 5% to 9% of the variance in the
true effect at its own location.

CBES *beats* that ceiling: +0.277 against the per-focus +0.222 under cluster/max, +0.331 against
+0.307 under FDR/max. Pooling across studies already extracts more than any individual
coordinate carries. The estimator is not underperforming its input -- the magnitude is not in a
coordinate table, and no change to the likelihood can put it there.

(The +1.000 correlation with the reporting study's own map at the focus is the sanity check that
the extraction and the conversion agree, since the reported value *is* that map's value there.)

### Landed

`nimare/meta/cbma/effectsize.py` now states the compression, the dependence on reporting
convention, the per-coordinate ceiling, and the list of remedies measured and rejected, so the
warning is no longer only about an unidentified constant.

### Held-out HCP under realistic reporting, and why pi*mu works

Re-running the held-out HCP design through `reporting.py` rather than a fixed cut with a cap.
Cluster-extent thresholding, one focus per surviving cluster at its centre of mass:

| design | 0-50% | 50-75% | 75-90% | 90-99% | 99-100% | overall | r |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20 subjects x 16, `g` | 10.67 | 3.24 | 1.90 | 1.25 | 1.27 | 1.66 | +0.519 |
| 20 subjects x 16, `pi*g` | 4.58 | 1.56 | 1.07 | 0.95 | 1.27 | 1.13 | +0.663 |
| 30 subjects x 12, `pi*g` | 4.23 | 1.36 | 0.87 | 0.81 | 1.16 | 1.00 | +0.649 |

Above the truth's 75th percentile `pi*g` is within 10% to 20% of a truth measured on subjects
the estimator never saw, and the pain top-stratum figure of 1.23 replicates here at 1.16 to 1.27.

It does not work for the reason the map's name gives. Every synthetic HCP study is drawn from
one population, so the true prevalence is exactly 1 and `pi*g` should equal `g`. Instead
`prevalence` comes back near 0.68. Multiplying by it is shrinking an inflated magnitude by a
data-driven factor, not averaging over studies that have no effect -- two documented biases
cancelling, `g` inflated upward by peak selection and `prevalence` compressed downward toward
the middle of its range. Worth using, not worth extrapolating.

Below the median it is still four to five times the truth, so the compression is intact.

### Landed

`g_marginal` reinstated in `nimare/meta/cbma/effectsize.py` with all three caveats stated, and
the PR description rewritten: the within-analysis null, the refusal of degenerate collections,
the corrected familywise error rates, the no-cap validation protocol, and the magnitude limits
above replacing the numbers that came from the capped extraction.

## One image beats ten coordinate tables, and mixing them costs

`do_coordinates_help.py` scores three estimates of the same held-out half of NIDM pain: inverse-
variance pooling of `k` image studies alone, CBES with those `k` as images and the rest of the
working half as coordinate tables, and CBES on the whole working half as tables. Eight splits,
cluster-extent reporting, one focus per surviving cluster.

| images | estimate | r | ratio | rmse |
| --- | --- | --- | --- | --- |
| 1 | images only | +0.520 | 1.62 | 0.446 |
| 1 | images + 9 coordinate studies | +0.083 | 1.92 | 0.666 |
| 2 | images only | +0.723 | 1.50 | 0.346 |
| 2 | images + 8 coordinate studies | +0.545 | 1.35 | 0.345 |
| 3 | images only | +0.774 | 1.34 | 0.276 |
| 3 | images + 7 coordinate studies | +0.614 | 1.25 | 0.285 |
| 5 | images only | +0.829 | 1.22 | 0.210 |
| 5 | images + 5 coordinate studies | +0.697 | 1.24 | 0.244 |
| 0 | 10 coordinate studies | +0.218 | 4.51 | 1.587 |

Adding coordinates costs correlation at every `k` and root mean square error at three of four.
A single image study beats nine coordinate tables *plus itself*, +0.520 against +0.083, and beats
ten coordinate studies alone by +0.52 to +0.22 at a quarter the error. The only thing the tables
improve is the level, and only slightly.

This is what a shared bias predicts and a shared noise would not. Per coordinate on this
collection the reported g at a focus is 1.965 +- 0.593 against a held-out truth of 0.519 +-
0.312: the error is a +1.45 bias against about 0.6 of noise, and the empirical error variance is
only about twice the nominal variance the estimator assigns. Noise averages out with more
studies. A bias every coordinate study shares does not -- pooling more of them converges harder
on the wrong number, and mixing them with unbiased images drags the pool toward it in proportion
to their weight.

### Corrected: the calibration explains the anomaly, not the verdict

The table above ran the mixed fit with `peak_bias=None`, so coordinates entered on their
inflated scale against images on the true one with nothing reconciling them -- the configuration
the class docstring warns against for mixed collections, which ship `peak_bias='per-study'` with
`peak_bias_scale='images'` for exactly this case. Rerun with both settings:

| images | estimate | r | ratio | rmse |
| --- | --- | --- | --- | --- |
| 1 | images only | +0.520 | 1.62 | 0.446 |
| 1 | mixed, uncalibrated | +0.083 | 1.92 | 0.666 |
| 1 | mixed, scaled to images | +0.320 | 1.66 | 0.506 |
| 2 | images only | +0.723 | 1.50 | 0.346 |
| 2 | mixed, uncalibrated | +0.545 | 1.35 | 0.345 |
| 2 | mixed, scaled to images | +0.589 | 1.40 | 0.373 |
| 3 | images only | +0.774 | 1.34 | 0.276 |
| 3 | mixed, scaled to images | +0.602 | 1.29 | 0.331 |
| 5 | images only | +0.829 | 1.22 | 0.210 |
| 5 | mixed, scaled to images | +0.662 | 1.24 | 0.274 |

The calibration accounts for the one genuinely strange number -- at a single image the mix scored
below both of its parts, and rescaling lifts it from +0.083 to +0.320. It accounts for nothing
else: at two images and above the two mixed rows are within 0.02 of each other, and images alone
still win at every count on both correlation and error. Mixing costs about 0.13 to 0.17 of
correlation and that is not an artefact of the scale.

What this changes is the remedy, not the finding. An earlier note here proposed that coordinates
stop feeding `g` while continuing to feed the null -- which is incoherent, because the
permutation null *is* a null for the `g` map: if shuffling coordinate magnitudes does not move
`g`, those studies contribute nothing to inference either. The coherent version leaves pooling,
coverage and the null alone and redefines the one map that is already conditional on images:
`g_absolute`, currently the mixed fit with its constant pinned, measured worse than an
inverse-variance pool of the image studies by themselves.

Three limits on the result: eight splits of one collection; every pain study is whole-brain, so
the tables added no spatial coverage that the images did not already have, which would not hold
against ROI images; and pain's median N of 16 makes a single image noisy, though still unbiased.

What this does *not* score is what the coordinate pathway uniquely provides -- localisation
inference under the within-analysis null, whose error rates are valid, and `prevalence`. Neither
depends on the magnitude scale, so the finding bears on where `g` should come from, not on
whether coordinates belong in the estimator.

### Why the earlier mixed-collection result said the opposite

`compare_mixed.py` and `mixed_weight.py` concluded that coordinates *improve* the magnitude when
images are present:

```
                     procedure    r all    r top      mag
                5 images alone    0.756    0.756     0.76
        CBES, coordinates only    0.788    0.804     2.04
       CBES, 5 images + coords    0.896    0.908     0.98
```

Four assumptions were wrong. One is fatal.

**The reference contained the studies being scored.** Both scripts state it in their own output:
`truth from all 21 images`. The five image studies are part of that truth, and the sixteen
coordinate studies' peaks were extracted from maps that are also in it. "Adding sixteen
coordinate studies improves prediction" therefore means, for the most part, that telling the
estimator about sixteen studies helps it predict a reference built from those sixteen studies.
It is close to tautological, and it is precondition B0.3 of the requirements.

**Scoring was restricted to the top quartile of the truth.** `mag = 0.98` is a median ratio over
the strongest 25% of voxels. That is the one stratum where the ratio really is near 1 -- measured
later at 1.15 to 1.23 -- while the bottom half runs four to ten times the truth. Restricting to
the top quartile makes the compression invisible by construction.

**One realisation.** A single fixed assignment of pain_01 to pain_05 as the images, no splits and
no error bars, against eight random splits in the current design.

**`r all` was taken over every voxel with a nonzero truth**, where spatial smoothness alone
generates correlation, and the arm carrying twenty-one studies' spatial information gains more
from that than the arm carrying five.

The old experiment already contained its own refutation. Its weight sweep:

```
              5 images alone    0.756    0.756     0.76
         + coords, weight x1    0.896    0.908     0.98
         + coords, weight x5    0.849    0.858     1.45
        + coords, weight x29    0.790    0.802     2.05
```

As coordinate influence grows the correlation falls monotonically and the scale inflates toward
the coordinate-only value of 2.04. That is the dilution signature. It was read as "coordinates
are not being swamped, so the measured mis-weighting has no practical consequence" -- the
question the script was built to answer -- rather than as evidence that coordinate influence
degrades the estimate.

One difference between the two designs is not yet accounted for and should not be claimed as
settled: the old extraction took every local maximum above a fixed threshold, the current one
takes one focus per surviving cluster, about six per study. The current coordinate arm carries
less information, and some of the gap may be that rather than the reference. Testable, untested.

## Error rates under the within-analysis null, with the familywise metric actually correct

The numbers previously recorded here came from a log written *before* the familywise fix, whose
"voxel FWE" column read 0.0000 and 0.0001 -- the second of which is arithmetically impossible as
a mean of zeros and ones over 100 draws, which is what gave it away. The script was corrected and
never rerun. 100 simulations per cell, 200 permutations, nominal 0.05, binomial standard error
0.022.

| cell | studies | uncorrected (per-voxel) | voxel FWE (familywise) |
| --- | --- | --- | --- |
| global null, N 20-40 | 20 | 0.0113 | 0.070 |
| global null, N 10-1000 | 20 | 0.0125 | 0.060 |
| global null, N 10-1000 | 12 | 0.0120 | 0.060 |
| global null, 2 foci/study | 20 | 0.0007 | **0.180** |
| power at focal g = 0.8, N 20-40 | 30 | 0.710 | 0.030 |
| power at focal g = 0.8, 2 foci/study | 30 | 0.640 | 0.020 |

The first three cells hold their level. The heterogeneous-N cells matter most: that is the
condition under which the old across-study shuffle rejected at 96.7%.

**The fourth cell is a defect.** Twenty studies reporting two foci each give a familywise rate of
0.18, about six standard errors above nominal, and the guard does not fire -- 2^20 arrangements
clears its 10^4 threshold easily. Swapping two values inside a study is a tiny perturbation, so
the max-statistic distribution those arrangements generate is far narrower than the count
implies. The pairing is diagnostic: 0.0007 uncorrected, ultra-conservative, against 0.18
familywise. A too-narrow permutation distribution of the maximum puts the familywise cutoff too
low. The guard counts arrangements when it should measure the spread of the attained maxima.

Few foci per study is the realistic regime under cluster-extent reporting, so this is not a
corner case, and the separately measured clustered global null -- voxel 0.100 +- 0.039 against
0.000 for cluster size and 0.050 +- 0.028 for cluster mass -- is probably the same cause.

### The familywise defect is a degenerate null, not the tail fit

`why_fwe_fails.py` separates the two mechanisms that could produce a familywise rate of 0.18 at
two foci per study, by rerunning the same simulations with the generalized Pareto tail
approximation switched off.

| foci/study | tail fit | voxel FWE | null max CV | distinct maxima |
| --- | --- | --- | --- | --- |
| 2 | on | 0.150 | 0.032 | 6 |
| 2 | off | 0.150 | 0.032 | 6 |

Identical. The tail fit is not implicated, and Winkler et al. (2016) validate that method and
recommend it for familywise-corrected p-values in any case. The cause is the null itself: across
200 permutations the maximum statistic takes **six distinct values**, with a coefficient of
variation of 0.032. Swapping two values inside a study does not move the maximum. No method of
reading a p-value off a distribution that flat can control the error rate.

That fully specifies the fix. It has to be a refusal, not a different p-value calculation, and
the statistic to refuse on is already computed during correction: the number of distinct attained
maxima and their spread. The current guard counts 2^20 arrangements and waves through a null that
attains six values.

### Coordinates lose on localisation too, and the design is the realistic one

The mixing experiment only ever scored magnitude, which left open whether coordinates lose the
number but win the map. Rescored with a rank correlation and the area under the curve for
recovering the truth's top decile:

| images | estimate | r | rank r | AUC |
| --- | --- | --- | --- | --- |
| 1 | images only | +0.520 | **+0.514** | **0.820** |
| 1 | plus 9 coordinate studies | +0.320 | +0.287 | 0.750 |
| 2 | images only | +0.723 | **+0.706** | **0.905** |
| 2 | plus 8 coordinate studies | +0.545 | +0.630 | 0.868 |
| 3 | images only | +0.774 | **+0.761** | **0.920** |
| 3 | plus 7 coordinate studies | +0.614 | +0.687 | 0.889 |
| 5 | images only | +0.829 | **+0.817** | **0.936** |
| 5 | plus 5 coordinate studies | +0.697 | +0.742 | 0.914 |

Images alone win every column of *this* table at every count.

**CORRECTED (later session), because the sentence that used to stand here over-claimed.** It read
"they lose both, less badly for the map" -- concluding about magnitude from a table that contains
no magnitude metric. `r`, `rank r` and `AUC` are all invariant to the scale and offset of the
estimate: Pearson to any affine map, Spearman and AUC to any monotone one. Demonstrated by
multiplying an estimate by 0.6 and adding 0.4, which destroys its magnitude and leaves all three
identical to four decimals while `ratio` moves 1.378 to 1.606. So this table measures **spatial
pattern and ordering only**, and is structurally incapable of rewarding or punishing the level.

The magnitude columns are in the table above, and they split. Decomposing
`rmse^2 = (mean error)^2 + var(error)` on a truth mean of 0.519:

```
 img                 estimate  ratio   rmse  mean err  sd of err
   2              images only   1.50  0.346    +0.260      0.229
   2  mixed, scaled to images   1.40  0.373    +0.208      0.310
   3              images only   1.34  0.276    +0.176      0.212
   3  mixed, scaled to images   1.29  0.331    +0.151      0.295
   5              images only   1.22  0.210    +0.114      0.176
   5  mixed, scaled to images   1.24  0.274    +0.125      0.244
```

At two and three images the mixture has the **smaller mean error** and the **larger spread**, so
images win `rmse` because the variance penalty exceeds the bias gain -- not because the level is
better. At five images the mixture loses both.

So the accurate scoreboard: **images win pattern decisively at every count and total error at
every count; the mixture wins the level at one to three images.** Which is what the sentence two
paragraphs up already said -- "the only thing the tables improve is the level, and only slightly"
-- and the summary sentence contradicted it. The lesson is narrow and worth keeping: a table of
scale-invariant metrics cannot support a conclusion about scale, however many columns it has.

**The design is the realistic one, which I nearly mis-stated.** The image arm takes `work[:k]` and
the coordinate arm takes `work[k:]` -- *different* studies, with the truth from a third set never
used by either. So this is not "degrade an available image to coordinates". It is three images
against three images plus seven additional coordinate-only studies, scored on eleven held-out
studies. More than tripling the study count, in coordinate form, makes every metric worse.

Two caveats survive and they are where the remaining case for coordinates lives.

*Same territory.* Every pain study is whole-brain and shares a construct, so the coordinate
studies cover ground the images already cover. Coordinates reaching regions, tasks or populations
the images do not would add coverage this design cannot see. That is the strongest untested
argument for including them.

*Seven, not seventy.* A real coordinate meta-analysis has thirty to a hundred studies. The
bias-versus-noise decomposition says more will not help, since a shared bias does not average
down, but that is an extrapolation from seven.

### SDM-PSI: what can and cannot be said

The only head-to-head in this repository is void. `compare_sdm.py` scored against a truth built
from all 21 pain images -- the same images whose maps produced the coordinates both methods were
given -- extracted peaks at an uncorrected z of 3.29, and capped. All three preconditions in
`PROTOCOL.md` are broken, so the figures it produced (SDM 0.47 of the truth, CBES 1.99) should not
be quoted, and they have been removed from the pull request.

What survives is structural, and it is not one verdict but three.

*Against ALE, for localisation.* Nothing measured here contradicts SDM-PSI being a sound
coordinate-based method, and it has author validation against pooled individual data behind it.

*Against an image-based meta-analysis when two or more images exist.* The evidence against
coordinate-derived magnitude is a property of the input rather than of any estimator: reported
peak values carry a shared bias of about +1.45 against 0.6 of noise, and a shared bias does not
average down with more studies. SDM-PSI anchors its imputation on those same values, so it
inherits it.

*On absolute Hedges' g from coordinates.* One tabulated coordinate explains 5% to 9% of the
variance in the truth at its own location. Imputation interpolates under a model; it cannot
manufacture information the table does not contain. The ceiling is method-independent.

Two things it does better than this estimator. It uses the bound that unreported voxels sit below
threshold -- evidence of absence, which the measurements here identify as real information and
which CBES's censoring term currently contributes almost nothing toward. And multiple imputation
propagates uncertainty more honestly than a point estimate.

One thing it does worse. It emits a number that looks like an image-based result without
signalling that the scale is unearned. Deflating is the safer direction to be wrong in than
inflating, but the output does not tell a reader either way.

**A corrected comparison is blocked on tooling, not on design.** SDM-PSI is installed at
`/tmp/claude-0/sdm` with its pain inputs, but the driver that sequenced `sdm_core pp`, `mean` and
`mi` lived in a previous session's scratchpad and is lost, and the argument syntax is
undocumented. The design needed is otherwise settled: split the 21 studies, give ten to both
methods as cluster-extent coordinates with no cap, hold eleven back for the truth, and reduce the
imputations from 50 and raise the threads from 1 so the run finishes.

### Coordinates as a spatial prior: the third combination scheme to fail

If coordinates cannot supply magnitude but do carry location, the natural division of labour is to
let the image estimate carry the scale and let coordinate density say where to trust it. Tested
against studies neither arm saw, ten splits:

| images | estimator | r | rank r | AUC | rmse |
| --- | --- | --- | --- | --- | --- |
| 1 | images only | +0.482 | +0.375 | **0.826** | 0.359 |
| 1 | shrunk by coordinate density | +0.418 | +0.247 | 0.662 | 0.272 |
| 1 | hard coordinate gate | +0.413 | +0.246 | 0.662 | 0.290 |
| 1 | ORACLE gate on the truth | +0.750 | +0.856 | 0.902 | 0.265 |
| 5 | images only | +0.734 | +0.588 | **0.943** | 0.168 |
| 5 | shrunk by coordinate density | +0.338 | +0.166 | 0.585 | 0.288 |
| 5 | hard coordinate gate | +0.340 | +0.166 | 0.585 | 0.282 |
| 5 | ORACLE gate on the truth | +0.860 | +0.886 | 0.955 | 0.148 |

The gate costs 0.16 to 0.36 of area under the curve, and costs more the better the image estimate
is -- the more there is to damage. The one column that improves, root mean square error at a
single image, is shrinkage toward zero flattering itself against a mostly-zero truth.

**The oracle row matters as much as the failure.** A gate that knew where the truth was would lift
the area under the curve to 0.955 and the rank correlation to 0.886, so spatial gating is a good
idea and coordinate density is simply a bad gate. With five to nine coordinate studies the count
is mostly zero, one or two, so `c / (c + 2)` multiplies by 0, 0.33 or 0.5: it zeroes large regions
where the truth is moderate and nobody happened to place a focus, and keeps regions where somebody
reported a noise peak.

That is three schemes tested and failed -- naive pooling, pooling with the coordinates rescaled
onto the image scale, and gating. Coordinates degrade an image-based estimate however they are
attached to it. What survives is partition: coordinates answer where images are absent, and are
reported as their own quantity rather than modulating the image estimate.

### Where coordinates help: through their silence, not their values

The whole-brain scoreboard says images win pattern and total error while the mixture wins the
level. `where_do_coordinates_help.py` asks *where*, scoring the same split-half NIDM pain design
per voxel with two images held fixed, eight splits, and reporting signed error because the
absolute one can be gamed by shrinkage.

By decile of the held-out truth:

```
truth decile    voxels  |err| images  |err| mixed  signed images  signed mixed
           0     23520         0.212        0.177         +0.211        +0.176
           3     23512         0.175        0.134         +0.139        +0.095
           6     23512         0.209        0.151         +0.099        +0.025
           7     23520         0.228        0.165         +0.091        -0.003
           8     23520         0.251        0.205         +0.104        -0.023
           9     23520         0.262        0.303         +0.137        +0.006
```

At the top decile, where the `|x|` upward bias of an absolute value is negligible for both arms
and the comparison is clean, **the mixture is essentially unbiased (+0.006) against +0.137 for
two images alone** -- and its absolute error is *worse* (0.303 against 0.262). Bias down twentyfold,
variance up. That is the same bias-for-variance trade the `ratio` and `rmse` columns encode,
localised to the voxels where it can be read without an artefact.

By how many coordinate studies actually reported near the voxel -- and this is the surprise:

```
coord studies    voxels  |err| images  |err| mixed  signed images  signed mixed
            0    217916         0.204        0.163         +0.133        +0.061
            1     15545         0.260        0.254         +0.169        +0.127
            2      1525         0.318        0.318         +0.226        +0.212
            3       198         0.333        0.311         +0.273        +0.228
```

**The bias reduction is largest where no coordinate study reports at all** -- 54% at count 0,
against 6% at count 2. It shrinks monotonically as more coordinate studies report nearby. That is
the opposite of what the magnitude channel would predict, and it identifies the mechanism.

A coordinate study that reported nothing near a voxel still enters the **censoring term**: its
silence is evidence the effect there is small, and that corrects the upward bias of a two-image
estimate over the 218,000 voxels (87% of the mask) where no focus lands. Where a study *did*
report, its peak value carries the shared +1.45 inflation measured on this collection, which
pushes the estimate back up and cancels most of the benefit.

So the coordinate corpus **helps through its silence and hurts through its values**, and the two
act in opposite directions at every voxel. The net is positive on the level only because silence
covers most of the brain.

The design consequence is sharper than "coordinates dilute". The value of a coordinate table in a
mixed collection is almost entirely in the *censoring* term -- the fact that a study looked and
reported nothing -- and almost none of it is in the reported magnitudes. An estimator that used
coordinate tables only as presence/absence evidence, discarding the peak heights entirely, would
capture the benefit measured here and avoid the cost. That is a testable design and it has not
been tried.

### Silence-only coordinates dominate the shipping mixed configuration

If the corpus helps through silence and hurts through values, the design follows: keep the tables
in the censoring term, drop their magnitudes from the pooled mean. NIDM pain, 2 images held fixed,
8 splits, truth is the held-out half.

```
          estimate       r  rank r    AUC  mean err  err at top    rmse
       images only  +0.723  +0.706  0.905    +0.198      +0.221   0.346
     mixed (ships)  +0.589  +0.600  0.874    +0.160      +0.185   0.373
      silence only  +0.643  +0.676  0.884    +0.072      -0.020   0.277
  coordinates only  +0.218  +0.205  0.650    +1.451      +1.230   1.587
```

**Silence-only is better than the shipping configuration on all six point estimates, and
paired-significant on five of them.** Paired across the eight splits:

```
  vs mixed (ships)              difference  paired p        vs images only   difference  paired p
              r                     +0.054    0.0355                     r      -0.080    0.0009
         rank r                     +0.075    0.0035                rank r      -0.030    0.0527
            AUC                     +0.010    0.3453                   AUC      -0.021    0.0378
       mean err (closer to 0)       -0.088    0.0005              mean err      -0.126    0.0000
     err at top (closer to 0)       -0.205    0.0001            err at top      -0.241    0.0000
           rmse (lower)             -0.096    0.0004                  rmse      -0.069    0.0021
```

The AUC difference against the shipping configuration is **not** significant (p 0.35), so "wins on
all six" overstates it -- the honest claim is five of six, with the top-decile AUC indistinguishable.
Everything else against the shipping mix is significant and in the same direction.

Against images-only the trade is significant on both sides: better on all three magnitude
columns (p 0.0021 or below) and worse on `r` (p 0.0009) and `AUC` (p 0.0378), with `rank r`
borderline (p 0.0527).

**One caveat on all these p-values.** Eight splits drawn from 21 studies overlap heavily -- each
split uses 10 of them -- so the splits are not independent and the effective sample size is below
eight. The paired tests are therefore optimistic about significance, though the effect sizes on
the error columns (-0.205 at the top decile) are large enough that the direction is not in doubt.

**Against images-only it is a real trade.** Magnitude decisively better -- `rmse` 0.277 against
0.346, and top-decile bias **-0.020 against +0.221**, essentially unbiased -- and pattern worse,
`r` 0.643 against 0.723. So the answer depends on whether `g` is read as a number or as a map,
which is the same fork the metric audit exposed.

**This is not the flattening idea that was withdrawn.** That replaced every peak height with a
constant in a *coordinates-only* fit and took `r` from 0.230 to 0.04 with the top-decile AUC to
chance. Checked by reading the script rather than trusting the note: it builds studies with
`"points"` only and no `"images"` key, so the magnitudes were the sole source of localisation
there. Here the images supply the pattern. I proposed this variant before verifying that
distinction, which was the wrong order.

**The implementation cost is small, and not by accident.** The channels are already separate:
`_accumulate` builds the pooled mean from the focus table, while `_apply_selection_model` takes
the roster from `sample_sizes` and the silence geometry from the table. Filtering inside
`_accumulate` alone removes the magnitudes and leaves every study's silence intact, so the probe
is a subclass in `experiments/silence_only_coordinates.py` and nothing in the estimator changed.

Two consequences worth weighing rather than asserting. `peak_bias` and `peak_bias_scale` become
meaningless in this mode -- there are no coordinate magnitudes for a scale to act on -- which also
retires the calibration path, the all-donor crash that lived in it, and the two-donor requirement
for `g_absolute`. And `prevalence` is untouched, having never been a magnitude question.

The first implementation attempt crashed inside `_resolve_peak_bias_scale`, which fits the
coordinates alone to calibrate. That was informative rather than incidental: with the magnitudes
discarded there is nothing to calibrate, so the arm belongs at `peak_bias=None` and should never
reach that path.

### Prediction: silence-only should survive at one image where the shipping mix collapses

Asked whether the design works with a single image. Writing the reasoning down before running it.

The two-image floor exists for one reason: `_MIN_SCALE_DONORS = 2`, because
`peak_bias_scale="images"` reads the peak-height scale off the donors and one donor barely pins
it. **In silence-only mode there is no scale to pin** -- the coordinate magnitudes are discarded,
so `peak_bias` is off and the calibration path is never entered. The rationale for the floor
evaporates with it.

And the existing one-image numbers are where the shipping mix looks worst:

```
 images  estimate                       r  ratio   rmse
      1  images only               +0.520   1.62  0.446
      1  mixed, uncalibrated       +0.083   1.92  0.666
      1  mixed, scaled to images   +0.320   1.66  0.506
```

A single image plus nine coordinate tables scores +0.083 -- below *both* of its parts, which the
notes flagged as the one genuinely strange number in that table, and which the calibration only
partly rescues (+0.320 against +0.520 for the image alone).

So the prediction, in two parts:

1. **Silence-only should not collapse at one image.** The collapse is attributed to a scale fitted
   from a single donor, and there is no scale here.
2. **The magnitude gain should be *larger* at one image than at two.** One image has a ratio of
   1.62 and a mean error near +0.32, against 1.50 and +0.26 at two -- a noisier, more
   upward-biased baseline for silence shrinkage to correct. At two images the top-decile error
   went +0.221 to -0.020; at one it starts further from zero, so there is more to recover.

The failure mode that would refute both: if the pattern correlation at one image falls toward the
coordinates-only floor of +0.218, the silence geometry is carrying the map rather than the single
image, and one image is not enough to anchor it. That is the thing to look at first.
