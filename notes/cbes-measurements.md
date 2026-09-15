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
