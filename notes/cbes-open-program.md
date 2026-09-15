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

---

## U2 answered: three defects in how uncertainty is reported

Read from the code and confirmed arithmetically.

### U2a. `scale_interval_` is a sample range, and its error changes sign with donor count

It is set to `(min(per_donor), max(per_donor))`. A sample range *grows* with the sample size while
uncertainty about the common scale *shrinks* like one over the square root of it, so the two
diverge in opposite directions. Drawing per-donor estimates from a normal with mean 0.60 and
standard deviation 0.08, twenty thousand times:

| donors | reported (min, max) width | honest 95% interval width | ratio |
| --- | --- | --- | --- |
| **2** | 0.0909 | 0.1781 | **0.51** |
| 3 | 0.1351 | 0.1601 | 0.84 |
| 5 | 0.1866 | 0.1323 | 1.41 |
| 10 | 0.2467 | 0.0967 | 2.55 |
| 20 | 0.2985 | 0.0692 | 4.32 |

At the **two-donor floor the estimator requires for `g_absolute`**, and which the docstring
specifically argues for, the reported interval is **half** the honest width -- it understates
exactly where it matters most. By twenty donors it overstates by 4.3-fold, which also reverses
the incentive to supply more images. The fix is a standard error of the mean across donors rather
than their range.

### U2b. `se` does not include the scale's uncertainty

`g_absolute` is the same array as `g`, so its standard error is `g`'s, which is conditional on the
scale being exactly right. The description text says the magnitudes should be read as an order of
scale, which is honest prose, but the number a reader will use does not reflect it. An interval on
a partially identified quantity that omits the identification uncertainty is a misrepresentation
however carefully the surrounding paragraph is worded.

### U2c. `g_marginal` has no standard error at all

It is emitted as `fit["g"] * fit["prevalence"]` with no variance propagated. So the one magnitude
map with a checkable reference -- the only one sharing an estimand with an image-based
meta-analysis -- is the one carrying no uncertainty. Propagating it needs the covariance of the
prevalence and the magnitude, which the observed information already contains, since the
prevalence is profiled out of it by a Schur complement.

---

## E2 answered: the compression is affine but not invertible, and worse than compressed

`prevalence_calibration` sweeps a grid of true prevalence at three configurations, 12 simulations
per cell, 24 studies, read at the truth voxel.

| N range | effect | true pi | fraction reporting | pi_hat | sd |
| --- | --- | --- | --- | --- | --- |
| 20-40 | 0.8 | 0.00 | 0.035 | **0.218** | 0.213 |
| 20-40 | 0.8 | 0.40 | 0.215 | 0.530 | 0.279 |
| 20-40 | 0.8 | 1.00 | 0.427 | 0.901 | 0.136 |
| 20-40 | 0.4 | 0.00 | 0.035 | 0.218 | 0.213 |
| 20-40 | 0.4 | 0.40 | 0.031 | 0.157 | 0.217 |
| 20-40 | 0.4 | 0.80 | 0.069 | 0.427 | 0.317 |
| 20-40 | 0.4 | 1.00 | 0.080 | **0.302** | 0.216 |
| 10-200 | 0.8 | 0.00 | 0.035 | 0.144 | 0.184 |
| 10-200 | 0.8 | 1.00 | 0.503 | 0.936 | 0.082 |

Affine fits against the two candidate truths:

| N range | effect | slope vs true pi | intercept | r | slope vs reporting | r |
| --- | --- | --- | --- | --- | --- | --- |
| 20-40 | 0.8 | 0.719 | 0.216 | 0.994 | 1.831 | 0.994 |
| 20-40 | 0.4 | **0.176** | 0.173 | **0.684** | 3.619 | 0.789 |
| 10-200 | 0.8 | 0.804 | 0.063 | 0.981 | 1.675 | 0.966 |

Three findings.

**There is a floor of 0.14 to 0.22 at a true prevalence of zero.** The map reads a fifth of
studies having an effect where none do. That is the same phenomenon as the point-process null
floor, on a different scale.

**The relationship is strikingly affine for strong effects but the coefficients are not stable.**
Correlations of 0.994 and 0.981 mean a correction exists in principle; slopes of 0.719 against
0.804 and intercepts of 0.216 against 0.063, differing only in the sample-size range, mean you
would need to know the collection's effect size and N distribution to choose it -- and knowing the
effect size removes the reason to want the map.

**At a weak effect it is non-monotone.** 0.218, 0.194, 0.157, 0.270, 0.427, 0.302 as the truth
goes 0 to 1: it reads *lower* at a prevalence of 1.0 than at 0.8. The per-simulation standard
deviation, 0.22 to 0.32, exceeds the entire range of the means. For weak effects prevalence is
close to noise.

**The caveat on all of it.** This varies prevalence *across collections* at one voxel. The
docstring claims something else -- that ordering across *voxels within one map* survives -- which
is a different quantity and could still hold. `prevalence_within_map` tests the claim as made,
with four sites of known differing prevalence inside a single fit.

This also closes the power argument. Achievability needs prevalence on a real scale to say which
side of the 0.8 ceiling a planned study sits on, and none of the above supports putting it there,
least of all for the weak effects where the question is live.

## E2b — does the ordinal reading of `prevalence` hold *within* one map? (answered)

E2 varied prevalence across collections and read one voxel. The docstring makes a
different claim: that ordering across voxels *within a single map* survives the bias.
Tested directly (`experiments/prevalence_within_map.py`): four sites in one fit at
(-32,0,0), (-10,0,0), (12,0,0), (34,0,0) with true prevalences 0.25/0.50/0.75/1.00,
24 studies, U=3.2905, 4 mm localisation jitter, 16 simulations.

```
effect 0.8
     site prevalence     0.25     0.50     0.75     1.00
  fraction reporting    0.198    0.375    0.656    0.828
   fitted prevalence    0.360    0.473    0.874    0.930
           fitted sd    0.194    0.153    0.126    0.073
  Spearman within a map: +0.862 (sd 0.169); exact ordering in 50% of maps

effect 0.4
     site prevalence     0.25     0.50     0.75     1.00
  fraction reporting    0.044    0.086    0.112    0.141
   fitted prevalence    0.134    0.365    0.440    0.540
           fitted sd    0.186    0.277    0.247    0.319
  Spearman within a map: +0.586 (sd 0.341); exact ordering in 12% of maps
```

Findings.

1. **The ordinal claim is real but weak, and only at a strong effect.** +0.862 mean
   Spearman says the ranking is mostly right on average; exact ordering in half the
   maps says that in any *one* map you publish, a four-site ranking is as likely to be
   scrambled somewhere as not. "Ordering survives" is too strong a phrase for a
   coin-flip.
2. **At a weak effect the ordinal reading is not usable.** +0.586 with sd 0.341 and 12%
   exact ordering. The sd is the thing: individual maps land anywhere from anti-ordered
   to ordered. A user comparing two voxels in one weak-effect map is reading noise.
3. **The map is compressed at both ends, differently at each.** Strong effect: 0.25
   reads 0.36 (inflated, the E2 floor showing up within the map) and 1.00 reads 0.93
   (deflated, the HCP result). So the *range* is compressed from 0.75 wide to 0.57, and
   the compression is not a single monotone transform of the truth — its direction
   flips across the range. Any user who reads a prevalence difference as a magnitude of
   difference is over-reading by ~25% at the ends and under-reading in the middle.
4. **Fitted sd shrinks with prevalence** (0.194 → 0.073 strong, roughly flat weak).
   Rare sites are both biased up and noisier, which is the worst combination for the
   use the docstring invites: picking out which of two regions is *less* consistent.

Action: the docstring's prevalence caveat needs to say the ordinal reading is reliable
only on average over many maps and only for strong effects, not that "ordering
survives". Logged as task #44.

## E6 — the image/coordinate exchange rate depends on effect size (answered)

Rerun of the exchange sweep at a weak effect (contrast: the earlier sweep was strong).
Correlation with the truth, 3 seeds, basis ceiling 0.954:

```
                  0 coords        5 coords       15 coords       40 coords
     0 img        --          0.807 [0.025]    0.844 [0.018]    0.894 [0.012]
     1 img    0.804 [0.052]    0.832 [0.057]    0.829 [0.085]    0.893 [0.013]
     2 img    0.897 [0.004]    0.894 [0.015]    0.894 [0.009]    0.900 [0.027]
     4 img    0.933 [0.010]    0.922 [0.009]    0.925 [0.005]    0.925 [0.002]
```

1. **The exchange rate is not a constant of the method.** Strong effect: 1 image ≈ 15
   coordinate studies, and the 2nd image beat 40. Weak effect: 1 image ≈ 5 coordinate
   studies (0.804 vs 0.807), and the 2nd image ≈ 40 coordinate studies (0.897 vs
   0.894). Images are worth *less* relative to coordinates when the effect is weak,
   because a single weak image is itself noisy while coordinate counts still
   accumulate.
2. **There is one cell where coordinates genuinely help alongside images**: 1 image,
   0.804 → 0.893 with 40 coordinate studies (+0.089, far outside the [0.052] seed
   spread). This is the first positive joint result that survives an honest reference.
   It sits *outside* the scope the user set (≥2 images), which is the right scope: at 2
   images the gain is 0.897 → 0.900, inside noise, and at 4 images coordinates
   *cost* 0.008.
3. **So the earlier flat "coordinates never help alongside images" was measured at one
   effect size and one image count.** Corrected statement: coordinates help only when
   the image arm is information-starved — one image, weak effect. Once two or more
   images are present, at either effect size, coordinates add nothing and eventually
   subtract. The scoping decision (≥2 images) therefore lands exactly on the regime
   where the joint model has no advantage to offer, which is worth saying out loud
   rather than hoping the exchange table reads better at some other effect size.

## The testbed calibration I skipped, and what it cost

My own standing rule is that every fast testbed needs a calibration check before its numbers
mean anything. I set the first coverage bed up by guess -- a single blob of peak g 0.6 at 4 mm
voxels smoothed with sigma 2 voxels -- and launched 100 replications of nine arms before
measuring what regime that was. It was this one:

```
measured field FWHM         19 mm   (intended ~8)
RFT extent demand           70 voxels
studies reporting anything  18%
reporting studies per fit   2.1 of 12
replications unfittable     30%
```

So a nominally 12-study meta-analysis was a 2-study one, and the 70% that could be fitted were
selected for having happened to produce signal -- which inflates the very bias the run was
measuring. Two compounding errors, both invisible in the output, which read as a clean table.

The smoothing was the root cause and the arithmetic is elementary: smoothing white noise with a
Gaussian of sd `sigma` gives a field whose autocorrelation sd is `sigma*sqrt(2)`, so
`FWHM ~ 3.3*sigma` voxels, not `2.355*sigma`. At 4 mm voxels, sigma 2 is a 27 mm field. Real
studies report 6-10 mm. And the extent threshold scales steeply with smoothness, so the error
did not degrade the bed gently -- it moved it into a different experiment.

Calibrated bed (`calibrate_coverage_bed.py`), swept rather than guessed:

```
sigma  FWHM mm  min ext  peak g  % reporting  foci if any
  0.8      7.6     11.8    0.80         1.00         3.67
  1.0      9.5     19.8    0.80         0.97         2.55
  1.5     14.1     47.6    0.80         0.70         1.48
  2.0     19.4     92.6    0.80         0.30         1.11
```

sigma 0.8 with peak g 0.8 is the regime: every study reports, a table lists 3-4 clusters, which
is what real coordinate tables look like. Everything below that row is a different problem.

Rule tightened: a calibration check is not a smoke test that the script runs. It is a printed
measurement of the *reporting regime* -- fraction of studies reporting, foci per reporting study,
measured smoothness, extent demand -- refused unless it matches published practice, before any
estimate is read out of the bed.

## A prediction, recorded before the answer arrives

The coverage run's residual bias falls steeply with the number of image donors: +0.506 at 0,
+0.248 at 2, and (in the smoke run) +0.063 at 6 and +0.003 at 12, out of 12 studies. Two
mechanisms fit, and they are not the same thing.

**Donor count.** The scale constant is a median over per-donor ratios, so with 2 donors it is a
median of two numbers and noisy. Then the bias is an estimation error in the scale and should
fall roughly as the donor count grows, independent of how many coordinate studies sit alongside.

**Image fraction.** The pooled estimate weights image and coordinate contributions, so with 2 of
12 studies imaged, ten peak-height studies drag the estimate up; with 12 of 12 there is
essentially no coordinate contribution left to bias it. Then the bias tracks the *share* of
studies with images and the scale constant is incidental.

The arms already queued discriminate them. "24 studies, 6 images" holds the donor count at 6
while halving the image fraction from 50% to 25%.

- If donor count drives it: bias at 6-of-24 should match 6-of-12, about +0.06.
- If image fraction drives it: bias at 6-of-24 should be clearly worse than 6-of-12 and land
  somewhere around +0.12 to +0.20, between the 17% and 50% rows.

Which it is matters for advice. Under the first, "get two or three images" is sound guidance.
Under the second, what matters is that images are a large *share* of the collection, and two
images among twenty studies buys very little -- which would make the scoping decision
(at least 2 images) necessary but nowhere near sufficient.

## U1/U3/U4 — interval coverage against a known truth (answered)

Bed: 5 sites, calibrated regime (7.6 mm FWHM, every study reports, 3-4 foci each), prevalence 1
at every site so the conditional and marginal estimands coincide, truth 0.800 at the read-out
voxel, 100 replications, cluster-extent FWE reporting with one max-statistic focus per cluster.

```
arm                              mean g    bias  mean se  sd of g  se/sd  cover  half/truth
12 studies,  0 images, tau 0.0    1.306  +0.506    0.179    0.104   1.72   0.10       0.44
12 studies,  0 images, per-study  1.314  +0.514    0.193    0.133   1.45   0.16       0.47
12 studies,  2 images, tau 0.0    1.048  +0.248    0.162    0.101   1.61   0.67       0.40
12 studies,  6 images, tau 0.0    0.879  +0.079    0.092    0.082   1.12   0.85       0.23
12 studies, 12 images, tau 0.0    0.810  +0.010    0.056    0.056   1.00   0.96       0.14
12 studies,  0 images, tau 0.3    1.486  +0.686    0.413    0.355   1.16   0.58       1.01
```

**The pooling and the standard error are correct.** The 12-of-12 row has `se/sd` of exactly 1.00
and coverage 0.96 against a nominal 0.95. That is as good as a calibration check gets, and it
means the local random-effects pooling, the observed-information SE and the read-out path all
work. Whatever is wrong is not in that code.

**But that row is not a CBES result.** Donor coordinates are dropped, so with every study imaged
the coordinate table is empty and the fit is an IBMA behind CBES's interface -- a configuration
the estimator refuses outright when reached directly (task #47). It belongs in the table as the
*reference* row: this is what the machinery does when the inputs are honest.

**The failure is entirely in the coordinate channel, and it is bias, not width.** Every arm's
`se/sd` is at or above 1, so no interval is too narrow for the estimator's own variability.
Coverage nonetheless falls to 0.10 because the estimate is 63% high. `se/RMSE` -- reported width
against total error -- is 0.30 at 0 images, 0.60 at 2, 0.81 at 6, 1.00 at 12, and that is the
column that tracks coverage.

**U4 answered, and not the way I expected.** The fixed, low-biased tau2 does *not* produce
under-coverage: at tau 0.3 the coordinates-only `se` widens from 0.179 to 0.413 and `se/sd` stays
at 1.16. Coverage rises from 0.10 to 0.58 -- but by width, not by accuracy: `1.96*se/truth`
reaches 1.01, so the interval is as wide as the effect it is estimating. The heterogeneity arm is
the cleanest example in the whole program of coverage improving while the output gets less
useful.

**U5 answered.** Mismatch and width are separable here and the separation is complete: width is
right everywhere (`se/sd` ~1-1.7), and every coverage failure is attributable to bias. So there
is no version of "tighten the SE" or "loosen the SE" that fixes `g`; only removing the bias does,
and the only thing measured to remove it is images -- lots of them.

## E7 — which axis controls the `g_marginal` cancellation? (in progress)

The docstring is already honest that the cancellation is regime-dependent: "That the two happen
to cancel ... is why this is worth using and also why it should not be trusted beyond the regimes
it has been measured in." What it does not say is *which* axis moves it, and a user cannot act on
"beyond the regimes measured" without knowing what to look at in their own collection.

The coverage bed supplies the missing observation. There, as in the HCP design, the true
prevalence is exactly 1 -- but `prevalence` comes back at 0.994, not 0.68, so there is no
downward push to offset `g`'s +0.5 and `g_marginal` inherits the whole bias. The two beds differ
in *reporting density*: the coverage bed's read-out site is strong and nearly every study reports
a focus within the kernel, while in the HCP design most studies were silent at any given voxel.

Mechanism, stated so it can be wrong: `prevalence` is pushed down only by silence, and `g` is
pushed up by peak selection regardless. Where reporting is dense, `prevalence` saturates near 1,
the downward push disappears, and `g_marginal` equals `g` -- including its bias. Where reporting
is sparse the two offset. If that holds, `g_marginal` is least trustworthy at the strongest and
most consistently reported voxels, which are the ones readers actually interpret, and the
warning a user needs is about reporting density rather than about "regimes".

`marginal_cancellation.py` crosses true magnitude with true prevalence in one fit -- strong and
universal, strong and uncommon, weak and universal, weak and uncommon -- so reporting density
varies for reasons the truth tells apart, and scores `g_marginal` against the true marginal
effect `pi * mu`. A cancellation that is structural shows small marginal bias at all four sites;
one driven by density shows the bias tracking the reporting-density column.

## E8 — the prediction resolved: bias is the coordinate channel's weight share (answered)

Recorded before the discriminating arm ran: if the residual bias tracks the *image fraction* it
should be clearly worse at 6-of-24 than at 6-of-12, "somewhere around +0.12 to +0.20"; if it
tracks the *donor count* it should match 6-of-12 at about +0.06.

Measured at 6 of 24: **+0.176**. Same six donors, 2.2x the bias of 6-of-12. The donor-count
mechanism is dead; the image fraction is the variable.

It is stronger than that. The whole bias curve is one formula. If the coordinate channel carries
a fixed bias `b0` and is pooled by inverse variance against an unbiased image channel worth `r`
coordinate studies per study, the surviving bias is `b0` times the coordinate channel's weight
share:

```
    bias(f) = b0 * (1 - f) / ((1 - f) + r*f)
```

Fitted on the five measured points with `b0` pinned by the f=0 row, so one free parameter:

```
  image fraction  measured  weight-share model   resid
            0.00     0.506               0.506  -0.000
            0.17     0.248               0.243  +0.005
            0.25     0.176               0.181  -0.005
            0.50     0.079               0.079  -0.000
            1.00     0.010               0.000  +0.010
```

`b0 = 0.506`, `r = 5.40`. Max residual 0.010 on biases spanning 0.496.

What this buys.

1. **The coordinate channel is not being corrected by anything; it is being diluted.** The fit has
   no term for a correction, and adding images changes only the weight share. That is exactly why
   `peak_bias='per-study'` moved nothing: it operates inside a channel whose bias is common, and
   the model says only the channel's *weight* matters.
2. **The bias is predictable before any fit.** A user who knows their image fraction knows roughly
   how high the magnitude will run. Getting bias under 10% of the effect needs `f >= 0.5` -- half
   the collection imaged -- and at `f = 0.17` it is still +31%. The "at least 2 images" scoping is
   necessary and, on this evidence, nowhere near sufficient: two images among twenty studies is
   `f = 0.10` and leaves about +0.37.
3. **`r = 5.4` is a magnitude exchange rate**, and it is a different quantity from the pattern
   exchange rate measured earlier (1 image ~ 15 coordinate studies at a strong effect, ~5 at a
   weak one). Same order, different question -- worth keeping distinct rather than quoting one
   number for both.
4. **`b0` belongs to the reporting regime, not the estimator.** It is set by the threshold, the
   smoothness and the sample sizes, which is why no code change has moved it.

Design consequence worth putting to the maintainer: the estimator already knows `f` and could
warn with a number rather than a paragraph -- "6 of 24 studies supply images; on simulated
collections in this reporting regime the magnitude runs about 20% high at that share" -- instead
of the current qualitative warning that the scale will be "far too large".

Caveat on scope: `r` is specific to this bed's kernel, sample sizes and foci-per-study. The
functional form should carry; the constant should be re-measured before it is quoted anywhere.

### A second prediction from the weight-share model, recorded before the test reports

The model's `f` is a *weight* share, not a study count, and the two channels do not cover the
same voxels. An image contributes at every voxel; a coordinate study contributes only within its
kernel support of a reported focus. So the effective image share is **local**: high where few
studies reported a nearby focus, low where many did.

That predicts two things the test now running (`is_g_a_scale_error.py`, five sites at true weights
1.00 / 0.875 / 0.75 / 0.625 / 0.50) should show.

1. **With images present, the bias is worse at densely-reported sites.** Strong sites attract
   foci from most studies, so the coordinate channel holds a larger share of the local weight
   there and keeps more of its bias. The ratio `g/truth` should therefore *rise* with the true
   strength, not stay flat.
2. **A handful of images can make the relative map worse while making the absolute scale
   better.** At `f = 0` the coordinate bias applies everywhere; adding two images removes it
   unevenly -- most where reporting is sparse, least where it is dense -- so the *pattern* is
   distorted by an amount that did not exist before. If so, the 2-image arm should have a larger
   ratio spread across sites than either the 0-image or the 12-image arm: a U shape, worst in
   the middle.

Prediction 2 would also retire a puzzle. An earlier NIDM pain measurement had images+coordinates
correlating *worse* with a held-out truth than images alone (r 0.774 -> 0.614 at 3 images), which
I recorded as "coordinates degrade both magnitude and localisation" without a mechanism. Uneven
local dilution is a mechanism, and it predicts the effect is largest exactly where the image
count is small -- which is where that measurement was taken.

If instead the ratio is flat across sites in every arm, the coordinates-only bias is one
multiplicative constant, `g` is a legitimate relative map, and the "read it relatively" defence
in the PR stands as written. Either answer is worth having; they imply different documentation.

### The complete coverage table, and the weight-share fit on all of it

All twelve arms, 100 replications each, truth 0.800 at the read-out voxel, prevalence 1.

```
arm                              mean g    bias  mean se  sd of g  se/sd  cover  half/truth  report
12 studies,  0 images, tau 0.0    1.306  +0.506    0.179    0.104   1.72   0.10       0.44    1.00
12 studies,  0 images, per-study  1.314  +0.514    0.193    0.133   1.45   0.16       0.47    1.00
12 studies,  2 images, tau 0.0    1.048  +0.248    0.162    0.101   1.61   0.67       0.40    1.00
12 studies,  6 images, tau 0.0    0.879  +0.079    0.092    0.082   1.12   0.85       0.23    1.00
12 studies, 12 images, tau 0.0    0.810  +0.010    0.056    0.056   1.00   0.96       0.14    1.00
12 studies,  0 images, tau 0.3    1.486  +0.686    0.413    0.355   1.16   0.58       1.01    0.87
12 studies,  6 images, tau 0.3    0.890  +0.090    0.164    0.127   1.29   0.91       0.40    0.87
12 studies, 12 images, tau 0.3    0.817  +0.017    0.113    0.091   1.24   0.97       0.28    0.87
24 studies,  0 images, tau 0.0    1.289  +0.489    0.116    0.086   1.34   0.00       0.28    1.00
24 studies,  0 images, per-study  1.296  +0.496    0.133    0.108   1.24   0.01       0.33    1.00
24 studies,  6 images, tau 0.0    0.976  +0.176    0.093    0.073   1.27   0.56       0.23    1.00
24 studies, 24 images, tau 0.0    0.807  +0.007    0.039    0.036   1.07   0.94       0.10    1.00
```

Weight-share fit over all six `tau = 0` arms, two of them at `f = 1` from different study counts:

```
  f=0.00  measured  0.506  model  0.506
  f=0.17  measured  0.248  model  0.243
  f=0.25  measured  0.176  model  0.181
  f=0.50  measured  0.079  model  0.079
  f=1.00  measured  0.010  model  0.000     (12 studies)
  f=1.00  measured  0.007  model  0.000     (24 studies)
```

`b0 = 0.506`, `r = 5.40`, max residual 0.010. The two `f = 1` rows come from 12 and 24 studies
and give 0.010 and 0.007 -- an independent check that the study *count* does not enter, only the
share, which is what the model asserts and what the 6-of-12 versus 6-of-24 contrast first showed.

Projected onto shares real collections actually have:

```
  2 images of 20 studies (f=0.10)  bias 0.317  = 40% of a 0.80 effect
  3 images of 20 studies (f=0.15)  bias 0.259  = 32%
  2 images of 10 studies (f=0.20)  bias 0.215  = 27%
  5 images of 20 studies (f=0.25)  bias 0.181  = 23%
```

That is the number the "at least 2 images" scoping needs to be read against. It buys the
*calibration machinery* -- without a donor there is no scale at all -- but on these projections it
leaves the magnitude 25-40% high, and the interval will not cover.

Two further readings of the table.

**The two `per-study` rows are a matched pair and both say the same thing.** +0.506 -> +0.514 at
12 studies and +0.489 -> +0.496 at 24. The correction is not marginal, it is null, twice.

**Heterogeneity buys coverage by width in exactly the way the metrics-that-lie note warns.** At 12
coordinate-only studies, tau 0.3 raises coverage from 0.10 to 0.58 while raising `half/truth` from
0.44 to 1.01. The interval becomes as wide as the effect. Anyone reading coverage alone would
record heterogeneity as an improvement.

## E3 — the power-spread prediction, confirmed (answered)

Section 18's prediction, recorded before the run: the prevalence slope on the truth should rise
as the spread of study power grows, while the product's slope changes least. True `mu = 0.6`,
24 studies, prevalence swept 0.25 to 1.00, 24 replications per cell.

```
roster        fitted prevalence at true 0.25/0.50/0.75/1.00   slope    fitted g across the sweep
fixed         0.275  0.476  0.597  0.680                      +0.534   0.481 0.667 0.702 0.734
n varies      0.261  0.480  0.708  0.963                      +0.934   0.625 0.665 0.643 0.670
both vary     0.219  0.472  0.739  0.979                      +1.019   0.632 0.620 0.646 0.644
```

Confirmed, and the clearest evidence of it is the column I had not thought to predict.

**With a homogeneous roster, `g` is contaminated by the prevalence.** It should be 0.6 at every
cell, because the true magnitude does not change across the sweep. Under the fixed roster it runs
0.481 to 0.734 -- a 53% swing driven entirely by a parameter it is supposed to be separate from.
Once sample sizes vary it is flat at 0.62-0.67 across the whole sweep, and varying thresholds too
changes nothing further. That is the separation working, and it is a sharper diagnostic than the
prevalence slope because the target is a constant.

**The prevalence slope nearly doubles**, +0.534 to +1.019, and lands essentially on 1.

**The part of the prediction that was wrong**: I said the marginal's slope would stay roughly
constant, the product being identified all along. It went +0.733 to +1.118, a 53% improvement --
smaller than the prevalence's 91% but not "roughly constant". The product is better identified
than the factors under a homogeneous roster, not fully identified.

Practical consequence, and it is a cheerful one for once: real collections *do* have
heterogeneous sample sizes, so the regime that breaks the separation is not the common one. It is
now a checkable precondition rather than an unknown -- and one CBES could compute and report.

## E7 — the `g_marginal` cancellation is driven by reporting density (answered)

Crossing true magnitude with true prevalence in one fit, 24 coordinate-only studies, 40
replications. The target is the true marginal effect `pi * mu`.

```
site                 rep dens  true mu       g  g bias  true pi     pi  true marg  g_marg  marg bias
strong, universal        0.95     0.80   1.292  +0.492     1.00  0.997      0.800   1.289     +0.489
strong, uncommon         0.36     0.80   1.170  +0.370     0.40  0.557      0.320   0.640     +0.320
weak, universal          0.20     0.45   0.755  +0.305     1.00  0.847      0.450   0.645     +0.195
weak, uncommon           0.09     0.45   0.593  +0.143     0.40  0.619      0.180   0.422     +0.242
```

Confirmed as predicted: the absolute marginal bias tracks reporting density, worst at the
strong universal site (0.489 at density 0.95) and smallest at the sparse ones. At that site
`prevalence` reads 0.997, so there is nothing to offset `g`'s +0.492 and `g_marginal` inherits
all of it. The cancellation is not structural; it needs silence to work on.

One qualification I should state rather than let the headline stand alone: in *relative* terms
the ranking inverts -- 61% of the truth at the dense site against 134% at the sparse weak one --
so "worst where reporting is dense" is true of the absolute error and false of the relative one.
Which matters depends on whether a reader is comparing regions or quoting a number.

## The naive count beats the fitted prevalence at the only reading the docs endorse (answered)

Both estimators scored on the *same* simulated collections -- 24 coordinate studies, four sites
at true prevalences 0.25/0.50/0.75/1.00, 40 replications. The naive estimator is simply
`(studies with a focus within r mm) / (studies in the collection)`, which needs no likelihood,
no EM, no selection model and no threshold inference, and is biased low by construction.

```
effect 0.8   estimate at 0.25/0.50/0.75/1.00      mean|bias|  RMSE   rho    exact
  CBES prevalence     0.315 0.534 0.819 0.925          0.061  0.170  +0.820   42%
  naive, 10 mm        0.187 0.354 0.560 0.717          0.170  0.200  +0.954   80%
  naive, 15 mm        0.222 0.420 0.646 0.804          0.102  0.148  +0.951   80%

effect 0.5
  CBES prevalence     0.276 0.428 0.615 0.716          0.129  0.303  +0.570   22%
  naive, 15 mm        0.104 0.136 0.273 0.331          0.414  0.421  +0.877   60%

effect 0.4
  CBES prevalence     0.171 0.368 0.426 0.534          0.250  0.370  +0.542   20%
  naive, 15 mm        0.066 0.105 0.132 0.177          0.505  0.509  +0.619   25%
```

The result splits cleanly, and both halves matter.

**The censored likelihood does what it was built to do.** It corrects the naive count's downward
bias, and by a lot: mean absolute bias 0.061 against 0.102-0.170 at a strong effect, 0.129
against 0.414 at a moderate one, 0.250 against 0.505 at a weak one. Silence really is being read
as evidence, and reading it works. That is a genuine result for the model.

**And it loses, decisively, at the ordinal reading -- the only reading the documentation
endorses.** Rank correlation +0.820 against +0.954 at a strong effect; the four-site ordering is
exactly right in 42% of maps against 80%. At a moderate effect, 22% against 60%. Not one regime
where the fitted prevalence orders better. The likelihood buys bias correction with variance, and
the variance is what destroys the ranking.

So the docstring is in an awkward position of its own making. It tells the reader to ignore the
level and read the order, which is precisely the aspect where a count anyone could compute in
three lines does better. The defensible positions are (a) stand behind the level, which means
retracting "read it ordinally" and owning a calibration claim, or (b) emit the naive count too
and say plainly which to use for which purpose. Both are honest; the present combination is not.

On RMSE the two split by regime -- naive 15 mm wins at a strong effect (0.148 to 0.170) and CBES
wins at weak ones (0.303 to 0.421, 0.370 to 0.509) -- so RMSE alone would have made this look
like a draw and hidden the ordering gap entirely. Another entry for the metrics-that-lie list.

## Reported heights are not decoration; within a study they are actively harmful (answered)

16 coordinate-only studies, 30 replications, five sites. Locations, counts, study membership,
sample sizes and thresholds untouched; only the reported statistics altered.

```
height input        r(g, truth)  mean g at sites  mean pi  r(pi, truth)
as reported               0.436            0.916    0.844         0.481
study-flattened           0.515            0.923    0.888         0.469
all-flattened             0.421            0.906    0.897         0.494

paired against 'as reported':
  study-flattened   r changes +0.0796  (sd 0.1485, paired p 0.006)
  all-flattened     r changes -0.0148  (sd 0.1901, paired p 0.672)
```

Stronger than section 17 predicted. Replacing every focus in a study with that study's own mean
reported height -- destroying all within-study variation and nothing else -- **improves** the
correlation with the truth by +0.080, at paired p 0.006. And destroying height information
entirely costs nothing measurable (-0.015, p 0.67).

So the within-study spread of reported heights carries *negative* information. That is what one
should expect once the winner's curse is taken seriously: the differences between one study's
reported peaks are differences in how far each local maximum overshot its own threshold, which is
noise dressed as signal, and the estimator weights it as signal. The between-study level, which
survives flattening, carries what little there is.

This is the most direct evidence yet that the estimating equation is in the wrong channel, and it
comes with a cheap intervention: a `peak_bias`-like option that replaces each study's reported
magnitudes with their own mean would have improved the map in this bed. Worth trying on real
collections before proposing it.

`prevalence` is indifferent to all of this (0.481 / 0.469 / 0.494), which is consistent with it
being a count-driven quantity.

## Where the coordinates-only bias comes from, stage by stage (answered)

12 coordinate-only studies, 40 replications, truth 0.800 at the read-out voxel.

```
stage                                      value  cumulative    step
1. truth at the read-out voxel             0.800      +0.000
   truth where the foci actually landed    0.673      -0.127  -0.127   <- localisation
2. mean g implied by the reported peaks    0.993      +0.193  +0.320   <- winner's curse
3. pooled, selection_model='none'          1.318      +0.518  +0.324   <- weighting + kernel
4. pooled, shipped selection model         1.306      +0.506  -0.012   <- selection model
```

Two things here that I had wrong.

**The pooling step contributes as much as the winner's curse**, +0.324 against +0.320. Everything
written about this bias so far -- in the PR, in the docstring, in these notes -- attributes it to
peak-height selection. Half of it happens *after* the reported values are in hand, in the
inverse-variance weighting and the spatial kernel. That half is in code rather than in the
literature's reporting practice, so unlike the winner's curse it might be fixable. It is the most
promising lead in this whole program and I have no mechanism for it yet; candidates are that the
kernel up-weights foci that landed closer to the read-out voxel (which are the larger excursions,
since a bigger local noise peak pulls the maximum toward the true peak) and that the variance
used for weighting, `1/n + g^2/2n`, is computed from the *reported* magnitude. The second would
deflate rather than inflate, so it is probably the first. Logged as a task; it needs measuring,
not reasoning.

**The selection model does essentially nothing**: -0.012 on a +0.518 bias. Correctly signed --
a censoring correction should pull down -- and negligible in size. That matches the earlier
finding that `peak_bias='per-study'` moves nothing, and it means the zero-inflated censored
likelihood is earning its keep on `prevalence` (where it beats the naive count's bias by a factor
of three) and not on `g`.

**Localisation helps rather than hurts**, by -0.127: the reported foci land where the truth is
*lower* than at the read-out voxel, because the read-out voxel is the true peak and jitter can
only move away from it. So localisation error partially offsets the winner's curse. That is worth
knowing because the two are usually described together as though they compounded.

## Is the coordinates-only `g` one scale error or a distortion? (answered, and a prediction falsified)

Five sites at true g 0.800 / 0.700 / 0.600 / 0.500 / 0.400, 12 studies, 40 replications.

```
                   g / true g at the five sites            spread  within-map rho   exact
coordinates only   1.630 1.577 1.539 1.640 0.996            1.65x   +0.922 (0.137)   57%
2 of 12 imaged     1.301 1.260 1.185 1.072 1.013            1.28x   +0.957 (0.067)   65%
12 of 12 imaged    1.011 1.013 1.020 0.964 1.002            1.06x   +0.972 (0.045)   72%
```

**Mostly one constant, failing at the bottom of the range.** Across the four stronger sites the
ratio runs 1.539 to 1.640 -- a spread of 1.07x, near enough to a single multiplicative constant
that the "read it relatively" defence holds there. The 1.65x headline spread is produced entirely
by the weakest site, where the ratio collapses to 0.996. So the relative map is trustworthy over
the part of the range that is reported consistently and breaks where reporting gets sparse, which
is the same boundary `g_marginal` was already known to have ("calibrates the upper part of the
range").

**My prediction was wrong.** I predicted a U shape -- that a couple of images would make the
*relative* map worse, because images remove the coordinate bias unevenly (most where reporting is
sparse, least where it is dense), distorting a pattern that was at least uniformly wrong before.
The spread instead falls monotonically, 1.65x to 1.28x to 1.06x, and the within-map ordering
improves at every step. Uneven local dilution does not distort the pattern, so it also does not
explain the old NIDM puzzle where images+coordinates correlated worse with a held-out truth than
images alone. That puzzle is still open (task #37).

**`g` orders voxels better than `prevalence` does.** 57% exact ordering over five sites against
`prevalence`'s 50% over four (a strictly easier task), and rank correlation +0.922 against +0.862.
If one of the two maps is to be recommended for a relative reading, the evidence favours `g` --
which is the opposite of what the current documentation does, since it tells the reader to
distrust `g`'s magnitude and read `prevalence` ordinally.

> **Superseded.** Every number in this section was measured on the known-variance z convention.
> Re-run with a genuine t and a supplied threshold, the ordering is far worse for both maps and
> the margin nearly disappears: `g` exact ordering 22% over five sites against `prevalence`'s 19%
> over four. The recommendation above no longer follows from a 3-point gap on differently-sized
> tasks. The ratio-spread conclusion also changes; see the corrected section below.

# CORRECTION: my simulators reported a known-variance z, the estimator expects a t

This invalidates the magnitude of every bias number recorded above for `g`. Read this before
using any of them.

**What happened.** My beds generate a study's map as `g_true + noise/sqrt(n)` with unit-variance
noise, then multiply by `sqrt(n)` to get a statistic. That is a *normal* statistic with known
variance: `z = d * sqrt(n)` exactly. The estimator receives it as `Z` and, correctly for real
data, treats a reported z as a p-value-preserving image of a **t** on `n - 1` degrees of freedom:
it maps z back to t and then to d. In the far tail -- and every reported peak is in the far
tail -- that map is strongly expansive. At `z = 5.33`, `n = 30`: `z / sqrt(n) = 0.98` while the
estimator returns `g = 1.25`, 27% larger.

So my simulator and the estimator disagreed about what the reported number means, and the
estimator is the one that is right for real data: fMRI group maps are t-maps, and a published z
is nearly always a transformed t.

**Measured, at the same foci, with the same reporting pipeline:**

```
statistic convention                  foci  mean stat  truth there  mean g    bias
known-variance z (what my beds did)    430      5.389        0.672   1.302  +0.630
proper t, reported as T                428      5.420        0.670   0.966  +0.296
```

Less than half. Roughly 0.33 of the bias attributed to the estimator all session is a convention
error of mine.

**What this invalidates.** Every quantitative claim about the size of `g`'s bias:

- the coverage table's `bias` column, and therefore its coverage numbers -- a bias near +0.17
  (0.296 at the foci less 0.127 of localisation) against an `se` near 0.18 does not give 0.10
  coverage, it gives something respectable
- `b0 = 0.506` in the weight-share fit
- the stage decomposition, whose "conversion convexity +0.323" step *is* this error, not a
  property of the estimator (and the Jensen term inside it is only +0.043 -- the rest is the
  level shift of the z-to-t map)
- the ratios in the scale-error test, and the absolute biases in the marginal-cancellation test

**What survives**, because it is a comparison made under one convention rather than an absolute:

- `se/sd` near 1 in every coverage arm -- that is a statement about the standard error, and the
  standard error was not converted through anything
- the *functional form* of the weight-share model, `bias(f) = b0 (1-f) / ((1-f) + r f)`, and the
  finding that the image *share* rather than the donor count is what matters; `b0` and `r` must
  be re-measured
- the naive-count-versus-fitted-prevalence comparison: both estimators read the same collections,
  so the ordering gap (80% against 42%) is a like-for-like result
- flattening reported heights improving the map by +0.080 -- again a within-convention contrast
- the power-spread result (E3), for the same reason
- the occupancy-model window of detectability, which never went through CBES's conversion at all

**The rule I broke.** My own note says suspect the test before the theory. For most of a session I
had a large, stable, reproducible bias and went looking for mechanisms in the estimator --
inventing a "pooling step" contribution, then a "conversion convexity" one -- when the first
question should have been whether the number I was feeding in meant what the estimator thought it
meant. Task #30 in this program was *exactly this bug* in the field simulator, fixed earlier and
then reintroduced in a new bed. A convention mismatch is not a subtle failure mode here; it is
the recurring one, and it deserves a standing check rather than vigilance.

**And one genuine finding falls out of it.** The z-to-t map at a reported peak is a long
extrapolation into the t tail and is therefore very sensitive to the assumed degrees of freedom.
A 27% amplification at `n = 30` means an error in the effective df propagates strongly into `g`.
Papers do not always report the df behind a z-map, and software differs in what it puts there.
That is a real caveat about the estimator, it is independent of my bug, and it is worth measuring.

## One thing in the estimator that is independently corroborated: `peak_information`

Worth recording separately, because most of this program's findings are cautions and this one is
not.

`peak_information` compares the mean reported peak height against what peaks of pure noise would
average at the same threshold, and warns when the excess is small that "their magnitudes are
close to uninformative about the effect size". In the corrected coverage bed -- a *favourable*
case: peak g of 0.80, peak t around 4.4, every study reporting 3-4 clusters, statistics on the
right convention -- the excess comes out at **+0.17 to +0.23 z** across 112 fits, which is below
the threshold at which it warns. So it warns, on data where the effect is strong and universally
detected.

Independently, and by a route that shares nothing with it, `are_heights_decoration` destroyed the
within-study variation in reported heights and found the map got **better** by +0.080 (paired
p 0.006), while destroying height information entirely cost nothing measurable.

Two unrelated measurements, one conclusion: the reported heights carry almost nothing, and what
within-study variation they carry is noise the estimator weights as signal. The diagnostic is
telling the truth, it fires in the right regime, and it fires even when conditions are good --
which is the correct behaviour for a warning about an information channel that is closed by the
reporting practice rather than by the collection at hand.

It also sharpens section 17 rather than softening it. I had half expected the excess to look
healthier once the statistic convention was fixed, since the convention error was inflating the
*recovered* magnitudes. It does not: the diagnostic works on the z scale, where the convention
error never applied, so its verdict was never affected by my bug. That is a small piece of luck
and a reason to trust it.

## The weight-share model survives the convention correction, with the shape parameter intact

Refitted on the corrected (t-convention) arms:

```
  image fraction  measured   model   resid
            0.00     0.259   0.260  -0.001
            0.17     0.127   0.121  +0.006
            0.50     0.027   0.039  -0.012
            1.00    -0.018   0.000  -0.018
```

```
                   b0       r
z convention     0.506    5.40
t convention     0.260    5.74
```

This is the cleanest confirmation available that the *form* was right and the *constant* was
corrupted. `b0` halves, exactly as the convention error predicts -- it is the raw coordinate bias,
which is the quantity my simulator was inflating. `r` barely moves, 5.40 to 5.74: an image is
worth about five and a half coordinate studies for magnitude, in both beds. A wrong functional
form would not have separated that way; both parameters would have shifted to absorb the change.

Revised projections onto shares real collections have:

```
  2 images of 20 studies (f=0.10)  bias 0.159  = 20% of a 0.80 effect
  3 images of 20 studies (f=0.15)  bias 0.129  = 16%
  2 images of 10 studies (f=0.20)  bias 0.107  = 13%
  5 images of 20 studies (f=0.25)  bias 0.089  = 11%
```

Much less alarming than the 25-40% the uncorrected fit implied, and it changes the practical
reading of the scoping decision. "At least 2 images" now leaves the magnitude 13-20% high at
realistic collection sizes rather than 27-40%. That is a caveat rather than a disqualification,
and it is in the range where a documented warning is a reasonable response -- which the earlier
numbers were not.

The qualitative points stand unchanged: the image *share* rather than the donor count is what
matters, no correction inside the coordinate channel moves the bias (`peak_bias='per-study'`:
+0.259 to +0.255), and the channel is diluted rather than corrected.

## E5 — achievability as an output (answered)

The user's question was "how can I use this estimand to power my next study?" My earlier answer
was unsatisfying: power under zero inflation needs `pi` and `mu` separately, both are badly
estimated, and there is a ceiling at `pi + (1-pi)*alpha`. Correct, and not usable.

The reporting-probability framing (section 19) gives a usable answer, because the quantity a
planner wants *is* the estimand. The chance that a planned study of `N` subjects reporting at
threshold `u` produces a surviving cluster at this voxel is

    P(N) = pi * D(mu, N, u)

which is the same reporting probability, evaluated at the planned design rather than the
collection's. And "produces a surviving cluster near here" is exactly what a replication counts
as a hit, so this is not a proxy for the planning question -- it is the planning question.

```
   site (mu, pi)    N=20    N=30    N=50    N=80   N=150   N=400  ceiling
    (0.80, 1.00)    0.61    0.86    0.99    1.00    1.00    1.00     1.00
    (0.60, 0.75)    0.20    0.37    0.62    0.74    0.75    0.75     0.75
    (0.45, 1.00)    0.10    0.20    0.46    0.77    0.99    1.00     1.00
    (0.45, 0.40)    0.04    0.08    0.18    0.31    0.39    0.40     0.40
    (0.25, 0.75)    0.01    0.02    0.05    0.11    0.31    0.72     0.75
```

Sample size for a target, where reachable:

```
  mu=0.80, pi=1.00   50%: N = 17    80%: N = 27
  mu=0.60, pi=0.75   50%: N = 39    80%: unreachable (ceiling 75%)
  mu=0.45, pi=1.00   50%: N = 54    80%: N = 85
  mu=0.45, pi=0.40   50%: unreachable (ceiling 40%)
```

Three things worth saying about this as a deliverable.

**The ceiling is `pi`, and it is the first thing a planner needs.** However large the study, it
cannot report an effect that is not there in the population it samples. A site with prevalence
0.40 caps at a 40% replication rate and no sample size fixes it. That is a qualitatively
different answer from a power calculation, and giving it first would head off the commonest
mistake -- planning N against an effect size as though detection were the only obstacle.

**It reframes "failed replication".** A study that does not find a peak where a meta-analysis
said one was is not necessarily underpowered or wrong; at `(mu 0.45, pi 0.40)` it had a 92% chance
of missing at N = 30. A map of `P(N)` at the reader's own planned design turns that from a
judgement into an expectation.

**But it inherits the identification problem, and only partly.** `P` is identified where the
factors are not (this is the whole point of section 19 -- the product is what the data pin down),
so `P` at designs *inside* the collection's range of sample sizes is estimable to about 0.04. The
ceiling, however, is `pi` alone -- and that is the badly identified factor, only separable inside
the window of detectability. So the honest output is a graded one: `P` at designs the collection
supports, reported confidently; the ceiling reported only where the collection's spread of power
identifies `pi`, and refused elsewhere. That is exactly the diagnostic in task #49, arriving from
a second direction, which is some reason to think it is the right diagnostic.

## A prediction about threshold inference, recorded before the test

My own protocol says to hand the real height threshold to the estimator through metadata rather
than letting it be inferred from the smallest reported value. Every bed this session let it be
inferred -- `threshold="study-min"` is the default -- and that now matters more than it did
before, because section 22 showed `prevalence` is governed entirely by where the fitted magnitude
sits relative to the assumed cutoff.

The mechanism predicts a direction. Under cluster-extent reporting a reported focus is a
cluster's *maximum*, which sits well above the cluster-forming cut. So the smallest reported
value in a study is the smallest cluster maximum, not the threshold, and `study-min` -- even
with the order statistic undone -- should infer a cutoff **above** the true forming cut.

A cutoff inferred too high makes the fitted magnitude look closer to it, which makes more of the
silence explicable as censoring, which inflates the prevalence. So:

- `threshold="study-min"` should give a **higher** fitted prevalence than supplying the true
  cluster-forming cut as a float.
- In the coverage bed the true prevalence is 1, so an inflated prevalence is *closer* to the
  truth there -- the inference error and the truth happen to point the same way, which means the
  coverage bed cannot distinguish "right for the right reason" from "right by accident". The
  prevalence sweep bed, where the truth ranges over 0.25 to 1.00, can.
- `g` should move less, since it is identified mainly by the reported values rather than by the
  silence.

If instead `study-min` and the supplied cut agree closely, the order-statistic correction is
doing its job under extent-based reporting too, and one documented worry can be retired.

## Threshold inference is the whole story for `prevalence` (answered — and it retracts a caveat)

Prediction confirmed, and the consequence is much larger than the prediction.

20 coordinate-only studies, 30 replications, true magnitude 0.70, four sites with true
prevalences 0.25 / 0.50 / 0.75 / 1.00, cluster-extent reporting at a forming cut of z = 3.0902:

```
threshold setting           median cut   pi@0.25 pi@0.50 pi@0.75 pi@1.00   mean g
(the truth)                                 0.25    0.50    0.75    1.00     0.70
study-min (the default)          4.015     0.482   0.860   0.961   0.998    0.873
pooled-min                       3.532     0.320   0.588   0.809   0.992    0.911
the true forming cut             3.090     0.212   0.466   0.722   0.953    0.951
the library default 3.2905       3.291     0.248   0.501   0.752   0.972    0.935
```

**`study-min` infers a cut of 4.015 against a true 3.090** -- 0.93 z too high, exactly as
predicted, because under cluster-extent reporting the smallest reported value is the smallest
cluster *maximum*, not the threshold, and the order-statistic correction cannot know that.

**And with the right threshold, `prevalence` is nearly unbiased.** 0.212 / 0.466 / 0.722 / 0.953
against 0.25 / 0.50 / 0.75 / 1.00 -- errors of -0.04, -0.03, -0.03, -0.05, with no compression
worth naming. The default gives 0.482 / 0.860 / 0.961 / 0.998: inflated at every site and
squashed against 1 at the top.

So the documented pathology of `prevalence` -- "compressed toward the middle of the range", a
true 0.25 coming back as 0.49 to 0.60, "read it ordinally, not as a fraction" -- **is largely an
artefact of the default threshold inference, not of the censored likelihood.** That is a
retraction of a caveat I shipped to the docstring earlier today, whose measured numbers were all
taken at the default.

A fixed plausible constant beats inference. The library's own
`DEFAULT_REPORTING_THRESHOLD_Z = 3.2905` lands at 0.248 / 0.501 / 0.752 / 0.972, which is about
as good as supplying the truth -- because p < 0.001 is close to what studies actually use, and
being roughly right beats being precisely wrong.

**`g` moves the other way**, 0.873 at `study-min` against 0.951 at the true cut, on a truth of
0.70. Same mechanism, opposite sign: a higher assumed cutoff means more of the reported mass is
attributed to censoring, which pulls the magnitude down and pushes the prevalence up. So the two
outputs cannot both be optimised by one threshold choice, and the *product* should be the most
stable thing across the column -- 0.421, 0.536 (pooled-min, pi@0.50 x g), and so on. Worth
checking directly.

What this implies for earlier results, all of which were measured at the default:

- The within-map ordering numbers (50% exact at a strong effect, 12% at a weak one) are
  threshold-limited, not intrinsic. They must be re-measured with a supplied threshold before the
  docstring's ordinal caveat is believed.
- The naive-count comparison too. CBES's prevalence lost on ordering while winning on bias; with
  a correct threshold its bias nearly vanishes, so the ordering may well improve with it. That
  re-run is the next thing to do.
- The window-of-detectability mechanism (section 20) is unaffected -- it is about where the truth
  sits relative to the *real* threshold, and it is what explains why a mis-set threshold does so
  much damage.

### Re-measured on the correct convention, and the threshold finding is narrower than it looked

`prevalence_within_map` rerun with a genuine noncentral t, and with both threshold settings:

```
effect 0.8                        0.25    0.50    0.75    1.00   within-map rho   exact
  fraction reporting             0.180   0.383   0.573   0.729
  fitted, study-min              0.279   0.641   0.869   0.955   +0.762 (0.190)    19%
  fitted, fixed 3.2905           0.306   0.684   0.890   0.961   +0.762 (0.190)    19%

effect 0.4                        0.25    0.50    0.75    1.00
  fraction reporting             0.021   0.047   0.083   0.081
  fitted, study-min              0.191   0.342   0.396   0.464   +0.328 (0.602)     6%
  fitted, fixed 3.2905           0.224   0.402   0.468   0.540   +0.328 (0.602)     6%
```

Three corrections to what I said an hour ago.

**1. The ordering is worse than I shipped, not better.** On the correct convention the four-site
ranking is exactly right in **19%** of maps at a strong effect and **6%** at a weak one, against
the 50% and 12% I measured with the buggy statistic and put into the docstring. Rank correlation
+0.762 and +0.328 against +0.862 and +0.586. So the docstring caveat I shipped is right in
direction and too generous in degree; it needs the corrected numbers.

**2. The threshold setting does not touch the ordering at all.** +0.762 and 19% under both
settings, identical to three decimals. That is not a coincidence: changing the assumed cutoff
applies a roughly common inflation across sites, and a monotone transform preserves ranks. So the
threshold choice moves the *level* and leaves the *order* exactly where it was. Whatever is
wrong with the ordering is not the threshold.

**3. The threshold finding is narrower than the headline suggested.** In the cluster-extent bed
supplying the true cut brought prevalence to within 0.05 of the truth at every site. Here, with
*direct height* thresholding, `study-min` infers the cut correctly -- the smallest reported value
really is near the threshold when there are no clusters -- and a residual inflation remains
anyway: 0.684 against a true 0.50 with the fixed threshold. So:

> Threshold inference is badly wrong **under cluster-extent reporting**, where the smallest
> reported value is a cluster maximum, and fixing it largely fixes calibration in that regime.
> Under direct height thresholding the inference is already fine and the inflation that remains
> is the censored likelihood's own.

The two beds also differ in how much silence there is to explain -- reporting fractions of
0.18-0.73 here against much higher ones in the field bed -- and the inflation tracks that, which
is the section 22 mechanism again: more unexplained silence, more for the censoring term to
attribute, more inflation. That is the variable to isolate next, rather than declaring either bed
the representative one.

### The threshold trade-off across all three outputs, and why one choice still wins

Mean absolute error across the four sites (20 coordinate-only studies, 30 replications, true
magnitude 0.70 at every site, cluster-extent reporting at a forming cut of z = 3.0902):

```
threshold setting           median cut   prevalence      g   g_marginal
study-min (the default)          4.015        0.201  0.173        0.270
pooled-min                       3.532        0.056  0.211        0.166
the true forming cut             3.090        0.037  0.251        0.116
the library default 3.2905       3.291        0.008  0.235        0.132
```

No single choice is best for all three, so the trade-off is real -- but it is not symmetric, and
one choice still wins once the *reasons* are looked at rather than the numbers alone.

**`prevalence` and `g_marginal` both want a plausible fixed threshold.** The library's own
constant gives a mean absolute prevalence error of **0.008** -- 0.248 / 0.501 / 0.752 / 0.972
against 0.25 / 0.50 / 0.75 / 1.00 -- against 0.201 for the default inference. On `g_marginal` the
true cut and the fixed constant are the two best (0.116, 0.132) and the default is worst by a
factor of two.

**`g` prefers the default inference, for a bad reason.** `g` is biased *high* by peak selection,
and a cutoff assumed too high makes the censoring term attribute more of the reported mass to
truncation, which pulls the magnitude down. So `study-min`'s advantage on `g` (0.173 against
0.251) is an accidental partial cancellation of a different bias by a threshold error, not the
threshold being right. Relying on it means relying on two errors staying in proportion, which is
exactly the kind of coincidence this program has found breaking down whenever a regime changes.

`g` is also not constant across the four sites at any setting -- 0.866 / 0.803 / 0.876 / 0.947 at
`study-min`, 0.986 / 0.926 / 0.943 / 0.947 at the true cut, for a truth of 0.70 everywhere -- so
it remains contaminated by the prevalence whatever the threshold. The threshold choice is not
what fixes that; heterogeneous study power is (E3).

**Recommendation, with the evidence attached:** supply a plausible fixed threshold, or leave
`threshold` at the library constant, and treat `study-min` as appropriate only where tables came
from voxelwise-height thresholding -- the one regime where the smallest reported value really is
near the cut. On cluster-extent tables it is reliably too high by about 0.9 z, and it damages the
two outputs that are identified while flattering the one that is not.

## The complete corrected coverage table (final)

100 replications per arm, truth 0.800 at the read-out voxel, prevalence 1 at every site,
calibrated reporting regime, studies reporting a genuine t on `n - 1` degrees of freedom.

```
arm                              mean g    bias  mean se  sd of g  se/sd  cover  half/truth  report
12 studies,  0 images, tau 0.0    1.059  +0.259    0.164    0.075   2.19   0.72       0.40    1.00
12 studies,  0 images, per-study  1.055  +0.255    0.171    0.080   2.14   0.75       0.42    1.00
12 studies,  2 images, tau 0.0    0.927  +0.127    0.134    0.098   1.36   0.87       0.33    1.00
12 studies,  6 images, tau 0.0    0.827  +0.027    0.090    0.071   1.27   0.99       0.22    1.00
12 studies, 12 images, tau 0.0    0.782  -0.018    0.065    0.059   1.10   0.94       0.16    1.00
12 studies,  0 images, tau 0.3    1.102  +0.302    0.249    0.157   1.58   0.80       0.61    0.90
12 studies,  6 images, tau 0.3    0.812  +0.012    0.147    0.127   1.16   0.96       0.36    0.90
12 studies, 12 images, tau 0.3    0.785  -0.015    0.110    0.099   1.12   0.93       0.27    0.90
24 studies,  0 images, tau 0.0    1.048  +0.248    0.113    0.051   2.24   0.28       0.28    1.00
24 studies,  0 images, per-study  1.046  +0.246    0.118    0.059   2.01   0.35       0.29    1.00
24 studies,  6 images, tau 0.0    0.873  +0.073    0.080    0.058   1.38   0.89       0.20    1.00
24 studies, 24 images, tau 0.0    0.778  -0.022    0.045    0.041   1.11   0.97       0.11    1.00
```

Weight-share fit over all six `tau = 0` arms: `b0 = 0.260`, `r = 6.40`, max residual 0.022. `r`
sits around 6 rather than being tightly pinned -- it read 5.74 on four points and 6.40 on six --
so "an image is worth five or six coordinate studies for magnitude" is the honest statement. The
discriminating contrast holds a third time: the same six donors give +0.027 at `f = 0.50` and
+0.073 at `f = 0.25`.

**`se/RMSE` predicts coverage and `se/sd` does not**, cleanly across the table:

```
                       se/sd   se/RMSE   coverage
12 studies,  0 img      2.19      0.61       0.72
24 studies,  0 img      2.22      0.45       0.28
12 studies,  6 img      1.27      1.18       0.99
24 studies, 24 img      1.10      0.97       0.97
```

`se/sd` ranks the four arms 2.19 / 2.22 / 1.27 / 1.10 -- monotone in the wrong direction relative
to coverage. `se/RMSE` ranks them 0.61 / 0.45 / 1.18 / 0.97, in exactly coverage order. This is
the metrics-that-lie entry, now confirmed on the corrected data.

**Coverage still degrades with more studies** in the coordinates-only arms: 0.72 at 12 studies to
0.28 at 24, and 0.75 to 0.35 with `peak_bias='per-study'`. Milder than the uncorrected 0.10 to
0.00, and the same signature: a fixed bias with a shrinking interval. It is the one failure mode
that improves on every accuracy metric while the inference gets worse.

**The all-donor arms are now a clean reference.** 12 of 12: bias -0.018, `se/sd` 1.10, coverage
0.94. 24 of 24: -0.022, 1.11, 0.97. Both against nominal 0.95. The pooling, the observed-
information standard error and the read-out path are correct; every failure in this table is in
the coordinate channel.

### The naive-count comparison, re-run on the correct convention (confirmed)

Same collections, both estimators, studies now reporting a genuine noncentral t thresholded at
p < 0.001 on each study's own t scale (direct height thresholding, so `study-min` infers the cut
correctly here and the threshold finding does not apply):

```
effect 0.8    estimate at 0.25/0.50/0.75/1.00    mean|bias|   RMSE    rho    exact
  CBES prevalence     0.295 0.620 0.834 0.951         0.075  0.190  +0.804    40%
  naive, 15 mm        0.193 0.402 0.569 0.744         0.148  0.179  +0.957    82%

effect 0.5
  CBES prevalence     0.271 0.381 0.510 0.586         0.199  0.326  +0.523    18%
  naive, 15 mm        0.083 0.116 0.191 0.202         0.477  0.483  +0.646    35%

effect 0.4
  CBES prevalence     0.199 0.315 0.317 0.419         0.312  0.422  +0.259    12%
  naive, 15 mm        0.048 0.070 0.098 0.109         0.544  0.547  +0.445    25%
```

Unchanged and sharper. CBES wins on bias in every regime by a factor of two to two and a half
(0.075 against 0.148; 0.199 against 0.477; 0.312 against 0.544). The naive count wins on ordering
in every regime by roughly a factor of two (82% against 40%; 35% against 18%; 25% against 12%).

Two things this settles.

**The ordering gap is not a threshold problem and cannot be fixed by supplying one.** Ordering is
threshold-independent (a change of cutoff applies a roughly common inflation, preserving ranks),
and this bed's inference is already correct because the reporting is voxelwise. So the gap is
intrinsic: it is section 22's mechanism, a well-ordered count multiplied by a noisy correction.

**CBES's level is genuinely good where the effect is strong** -- mean absolute bias 0.075 on
prevalences spanning 0.25 to 1.00. That is a number worth standing behind, and it makes the
"read it ordinally, not as a fraction" instruction the weaker of the two available positions
rather than the safer one. The honest framing is the reverse of the current one: trust the level
in the strong-effect regime, and use a count if what you want is to rank regions.

## Is coordinates-only `g` one scale error? (corrected answer, and a falsification reversed)

Rerun with a genuine t on `n - 1` degrees of freedom, donor images derived from the same t through
the estimator's own conversion, and the reporting threshold supplied. Five sites at true g
0.800 / 0.700 / 0.600 / 0.500 / 0.400, 12 studies, 40 replications.

```
                   g / true g at the five sites            spread  within-map rho   exact
coordinates only   1.326 1.314 1.372 1.411 1.210            1.17x   +0.760 (0.338)   22%
2 of 12 imaged     1.130 1.108 1.044 0.905 0.828            1.36x   +0.942 (0.070)   52%
12 of 12 imaged    0.967 0.962 0.965 0.943 0.959            1.03x   +0.952 (0.071)   62%
```

Compare the uncorrected run, which read 1.65x / 1.28x / 1.06x and monotone.

**Coordinates-only `g` really is close to one scale error.** A spread of 1.17x across a
twofold range of true effect, with ratios between 1.21 and 1.41, means a single multiplicative
constant would repair most of the map. That is materially better than the 1.65x the buggy
convention suggested, and it is the strongest support the PR's "read `g_relative`, not `g`"
framing has had. The weakest site is no longer an outlier (1.21 against 1.33 at the strongest),
so the earlier caveat that the relative reading breaks at the bottom of the range does not
survive the correction either.

**And the U-shape I predicted and then recorded as falsified is real.** 1.17x with no images,
**1.36x with two**, 1.03x with all twelve -- non-monotone, worst in the middle, which is what I
wrote down in advance and then retracted on the strength of the uncorrected run. The mechanism
proposed then also holds: an image contributes at every voxel while a coordinate study
contributes only near its own foci, so the coordinate channel's *local* weight share is highest
where most studies reported -- at the strong sites. Those keep more of the coordinate bias, and
the ratio falls monotonically with the truth (1.13 at 0.80 down to 0.83 at 0.40).

**But the distortion coexists with better ordering, which is why I could not see it.** Exact
ordering goes 22% to 52% to 62% and rank correlation +0.76 to +0.94 to +0.95 -- monotone
improvement, right through the arm with the worst ratio spread. The tilt is monotone in the truth
and in the direction that *stretches* the map rather than compressing it, so it increases contrast
between strong and weak sites while making each ratio individually less accurate. A pattern
metric and a calibration metric therefore disagree about the 2-image arm, and I had only been
looking at one of them.

Two lessons, of different kinds. The scientific one: adding a couple of images buys a large gain
in ranking and a small loss in proportionality, and which matters depends on the reading. The
methodological one: a prediction "falsified" by a measurement that later turns out to have been
mis-instrumented is not falsified, and the retraction needs revisiting when the instrument is
fixed -- I had not gone back over the earlier falsifications after finding the convention bug, and
should have.

## Audit: which results rested on the known-variance z convention

Prompted by discovering that a "falsified" prediction was only falsified by a mis-instrumented
measurement. Rather than keep stumbling on these, here is the full list, checked by grepping for
beds that declare `kind: "Z"` while building a statistic as something times `sqrt(n)`.

Re-run on a genuine t, conclusions confirmed or corrected:

| bed | outcome |
| --- | --- |
| `interval_coverage` | redone; bias halves, coverage 0.72/0.28 rather than 0.10/0.00, `se/sd` and the weight-share *form* unchanged |
| `is_g_a_scale_error` | redone; ratio spread 1.17x not 1.65x, and the U-shape prediction reinstated |
| `prevalence_within_map` | redone; ordering *worse* (19% / 6% exact, not 50% / 12%) |
| `prevalence_vs_naive` | redone; conclusion unchanged and sharper |
| `are_heights_decoration` | redone; **conclusion reversed and retracted** |
| `decompose_bias` | its "conversion convexity" stage *is* this bug; the other stages stand |
| `identification_by_power_spread` (E3) | redone; confirmed more cleanly, plus a caveat that the bed has no winner's curse |
| `marginal_cancellation` (E7) | redone; confirmed and sharper |

Unaffected, and why:

- `threshold_sensitivity`, `calibrate_coverage_bed`, `heights_on_real_data` -- written after the
  fix, or on real data where the estimator's assumption is the correct one.
- `reporting_probability` and the occupancy/window-of-detectability work -- the exact detection
  model never passes through the estimator's conversion at all.
- Every earlier real-data result (held-out HCP, NIDM pain split-half, the NeuroVault paradigm
  sets) -- real collections report real t or z maps, so there is no mismatch to make.

**Audited clean: the error-rate beds never had it.** `within_analysis_null_rates` -- the bed
behind the PR's numbers -- uses NiMARE's own `create_effect_size_coordinate_studyset`, which emits
`observed_z = t_to_z(observed_d / scale, dof)`: a proper p-value-preserving image of a t on the
right degrees of freedom. `mixed_null` and its variants do it correctly by hand, as
`t_to_z(val * sqrt(n), n - 1)`. `few_images_null`, `approx_null_images` and `config_matrix` draw a
reported z directly from `rng.uniform(3.3, 5.0)`, and there is no true effect size in those beds
for a converted statistic to disagree with. So the PR's error rates stand.

Which sharpens the lesson rather than softening it: **`mixed_null` contains the correct idiom,
written in an earlier session.** I had the right pattern in this repository and wrote four new
beds today without reusing it. That is a regression in practice, not a gap in knowledge, and it
is why the convention now lives in one place (`reporting.study_t_field`) with a check every bed
runs (`reporting.assert_statistic_convention`) rather than in whichever line I happen to type.

The general point for the protocol: after finding an instrumentation bug, **go back over every
conclusion the instrument produced, including the ones that were negative.** A bug that inflates
an estimate also inflates the evidence against predictions that said the estimate would be
smaller, and those retractions need revisiting too. I found this one by accident rather than by
audit, which is the wrong way round.

## E3 and E7 re-run on the correct convention (both confirmed, one with a caveat I owe)

### E3 — power spread identifies the split

With a genuine noncentral t, per-study cuts derived from a reporting p on each study's own t
scale, and the threshold handed in rather than inferred. True `mu = 0.6` at every cell.

```
roster        fitted prevalence at true 0.25/0.50/0.75/1.00   slope   fitted g across the sweep
fixed         0.453  0.702  0.926  0.943                     +0.677   0.387 0.497 0.547 0.603
n varies      0.318  0.594  0.803  0.980                     +0.878   0.506 0.482 0.537 0.588
both vary     0.233  0.488  0.752  0.986                     +1.009   0.554 0.561 0.604 0.577
```

Confirmed, and more cleanly than on the buggy convention. The prevalence slope rises +0.677 to
+1.009 and lands on 1. The `g` column, which should be a constant 0.6, swings 0.387 to 0.603 under
the fixed roster -- a 56% range driven entirely by a parameter it is meant to be separate from --
and is flat at 0.55-0.60 once both sample sizes and thresholds vary. The marginal slope goes
+0.804 to +0.982, a 22% improvement against the prevalence's 49%, so my original prediction that
the product "stays roughly constant" while the factors improve is closer to right here than it
was on the uncorrected data (where it moved 53%).

**The caveat I owe on this bed.** With both varying, `g` comes back at 0.554-0.604 against a truth
of 0.600 -- apparently unbiased, which flatly contradicts the coverage bed's +0.26. The difference
is that this bed has **no winner's curse**: it draws one value at the site and jitters its
location, rather than taking the maximum of a smooth field over a surviving cluster. So E3
measures *identification* -- whether the prevalence/magnitude split can be recovered -- and not
selection bias. The two are separate and additive, and the honest statement is that heterogeneous
power stops `g` being contaminated by prevalence; it does nothing about the peak-height inflation,
which is what the coverage bed measures. I would have over-read this without the contrast.

### E7 — the `g_marginal` cancellation needs silence to work on

```
site                 rep dens  true mu      g  g bias  true pi     pi  true marg  g_marg  marg bias
strong, universal        0.93     0.80  1.062  +0.262     1.00  0.997      0.800   1.058     +0.258
strong, uncommon         0.39     0.80  1.042  +0.242     0.40  0.371      0.320   0.386     +0.066
weak, universal          0.28     0.45  0.618  +0.168     1.00  0.602      0.450   0.355     -0.095
weak, uncommon           0.12     0.45  0.647  +0.197     0.40  0.220      0.180   0.136     -0.044
```

Confirmed and sharper. At the densely-reported site the prevalence saturates at 0.997, so there is
nothing to offset `g`'s +0.262 and `g_marginal` inherits essentially all of it (+0.258). At the
other three the cancellation works well: +0.066, -0.095, -0.044 against truths of 0.32, 0.45 and
0.18. And the prevalence at the "strong, uncommon" site is accurate (0.371 against a true 0.40),
because that site sits inside the window of detectability.

So `g_marginal` is trustworthy where the effect is not reported by nearly everyone, and inherits
the full peak-height inflation where it is. Since "reported by nearly everyone" is what makes a
voxel interesting to a reader, that is a caveat to state in the docstring rather than a
reassurance.

## The caveat on the convergence comparison, written before the numbers arrive

The smoke run (2 splits) has `g_marginal` at rank correlation 0.553 and top-decile AUC 0.850
against ALE's 0.233 / 0.665 and MKDA's 0.375 / 0.698. That is the first strongly favourable
comparative result in this program, which is exactly when I should write the objections down
rather than after.

**The reference structurally favours the magnitude estimators.** The truth is an
inverse-variance pooling of held-out *images* -- a magnitude map, and specifically an estimate of
`pi * mu`, which earlier work established is the quantity an IBMA estimates and the one
`g_marginal` shares. ALE and MKDA estimate convergence, a different quantity. Scoring both
against a magnitude reference asks how well each approximates a magnitude, which is the question
one of them was built for. A convergence map that perfectly captured where studies agree would
still lose here wherever agreement and magnitude diverge.

So the fair reading, whatever the numbers say, is narrow: **if what a reader wants is a map of
where the effect is large, this comparison speaks to that; if what they want is a map of where
studies agree, ALE answers that and this comparison does not test it.** I should not report the
result as "CBES beats ALE" without that clause.

**Two things that are legitimate advantages rather than artefacts.** CBES uses the reported
heights and ALE/MKDA do not, and the heights carry real information (23a) -- that is a genuine
difference in inputs. And CBES models silence, which is information the convergence statistic
discards. Both are reasons to expect it to do better, and neither is a rigged comparison.

**One thing to check rather than assume.** The voxel set is CBES's own `n_studies > 0`, so all
estimators are scored on the same voxels but the set is chosen by one of them. If CBES's coverage
radius happens to exclude voxels where the convergence maps do well, that would flatter it.
Worth a sensitivity check on the union of coverage instead.

**And the size of the evidence.** One collection, 21 studies, ten splits of the same data. That
is a single source, which is what the retracted height-flattening claim also rested on before
the real-data check reversed it. A second collection -- the NeuroVault paradigm sets -- should
carry this before it goes near the PR.

## #38 — does a plain convergence map localise as well as CBES? (answered, and the answer splits)

NIDM pain, 21 studies, cluster-extent reporting with one max-statistic focus per cluster, 10
splits, held-out inverse-variance pooling as the reference, every estimator reading exactly the
same tables. Scored on localisation only, since a convergence statistic has no effect-size scale.

```
estimate              rank r vs truth   AUC top decile
CBES g_marginal                 0.425            0.763
CBES prevalence                 0.240            0.667
CBES g                          0.216            0.653
MKDA density / KDA              0.294            0.652
ALE                             0.203            0.643
```

The answer splits, and the split is the interesting part.

**CBES's magnitude map `g` does not localise better than a convergence statistic.** AUC 0.653
against MKDA's 0.652 and ALE's 0.643 -- level, on a reference the magnitude estimator is aimed at
and the convergence ones are not. All the censored-likelihood machinery, the selection model and
the effect-size conversion buy nothing here over a smoothing kernel.

**`g_marginal` does, by a clear margin.** AUC 0.763, a gain of 0.11 over the best convergence arm
and 0.11 over `g` itself. So what earns its keep is specifically the *product* of magnitude and
prevalence, not the magnitude machinery alone -- which is the same conclusion arrived at from the
estimand side (`g_marginal` shares an IBMA's estimand; `g` has no external reference) and from
the mechanism side (`g_marginal` is magnitude times an empirical reporting probability).

**`prevalence` alone is also level with convergence**, 0.667 against 0.652. That is unsurprising
-- it is the count channel expressed differently -- and it is a useful sanity check that the
comparison is measuring what it claims.

With the caveat pre-committed before these numbers: the reference is an image-based estimate of
`pi * mu`, which is exactly what `g_marginal` estimates and what ALE and MKDA do not, so some of
that 0.11 is definitional rather than earned. The defensible claim is the narrow one:

> **If what a reader wants is a map of where the effect is large, `g_marginal` does better than a
> convergence statistic on a held-out image reference, and `g` does not. If what they want is a
> map of where studies agree, ALE answers that question and this comparison does not test it.**

Still owed before this goes near the PR: a second collection (the NeuroVault paradigm sets), and
a sensitivity check on the voxel set, which is currently CBES's own coverage. A paired test across
splits is running, since the splits share studies and an unpaired comparison would be swamped by
the between-split variance.

### Paired across splits, which is the right test since the halves share studies

```
estimate              rank r   (sd)     AUC    (sd)
CBES g_marginal        0.425  0.094   0.763  0.097
MKDA density / KDA     0.294  0.083   0.652  0.056
CBES prevalence        0.240  0.114   0.667  0.057
CBES g                 0.216  0.086   0.653  0.083
ALE                    0.203  0.045   0.643  0.030

paired against MKDA density, the strongest convergence arm:
  CBES g           rank r  -0.078  (sd 0.146, p 0.128)     AUC  +0.001  (sd 0.092, p 0.980)
  CBES g_marginal  rank r  +0.131  (sd 0.097, p 0.002)     AUC  +0.111  (sd 0.068, p 0.001)
  CBES prevalence  rank r  -0.054  (sd 0.093, p 0.097)     AUC  +0.016  (sd 0.039, p 0.241)
  ALE              rank r  -0.091  (sd 0.077, p 0.005)     AUC  -0.009  (sd 0.036, p 0.433)
```

The split is statistically clean, not a matter of eyeballing means.

- **`g` against the best convergence arm: AUC +0.001, p 0.980.** Indistinguishable. Not "slightly
  worse" or "about the same" -- the point estimate of the difference is one thousandth of an AUC
  unit. The magnitude machinery buys nothing over a smoothing kernel on this measure.
- **`g_marginal`: +0.111 AUC at p 0.001 and +0.131 rank correlation at p 0.002.** Significant on
  both, with a paired sd of 0.068 on the AUC difference, so ten splits are enough to see it.
- **`prevalence`: +0.016 AUC, p 0.241.** Level, as expected for a count-channel quantity scored
  against convergence statistics.
- ALE sits slightly below MKDA on rank correlation (-0.091, p 0.005) and level on AUC.

So the answer to #38, with the pre-committed caveat attached: the product earns its keep and
neither factor does, on a reference that is itself the product. That is three independent routes
to the same conclusion -- estimand, mechanism, and now localisation -- and it is the strongest
case in this program for treating `g_marginal` as the deliverable and `g` as a diagnostic.

### The voxel set reverses it, and the pre-committed worry was the right one

```
                    AUC on CBES's coverage    AUC on the whole mask
CBES g_marginal              0.763                    0.644
CBES prevalence              0.667                    0.641
CBES g                       0.653                    0.641
MKDA density / KDA           0.652                    0.669
ALE                          0.643                    0.753

paired against MKDA, whole mask:
  CBES g           AUC  -0.028  (sd 0.005, p 0.000)
  CBES g_marginal  AUC  -0.025  (sd 0.005, p 0.000)
  ALE              AUC  +0.084  (sd 0.015, p 0.000)
```

Complete reversal. On CBES's own support `g_marginal` wins by +0.111; over the whole brain
**ALE is the best arm and all three CBES maps are the worst**, every difference at p < 0.001.
The concern I wrote down before the numbers -- that scoring on CBES's coverage could flatter it --
is exactly what happened.

**The decisive number, which separates artefact from limitation.** CBES assigns exactly 0 outside
its coverage, and a tie block takes mid-ranks and contributes nothing to an AUC, so the whole-mask
penalty could have been pure scoring artefact. Measured over the same ten splits:

```
  share of the brain CBES covers at all:         0.092
  share of the truth's top decile CBES covers:   0.341
  enrichment of the top decile in covered voxels: 3.70x
```

So it is not an artefact. CBES's coverage is genuinely aimed at signal -- a 3.7x enrichment is
substantial -- but **two-thirds of the truth's strongest voxels lie where CBES returns zero.** It
is not declining to rank noise; it is declining to estimate most of where the effect is.

Part of that is a genuine mismatch of objects rather than a fault: an inverse-variance pooling of
z maps is a smooth, spatially extensive field, while a coordinate table is a few dozen points, and
a 20 mm sphere around each covers 9% of the brain. A reader wanting the extent of an effect will
not get it from coordinates at any threshold. But the consequence for the comparison is
unavoidable: ALE's kernel assigns a decaying value everywhere and so ranks the uncovered
two-thirds, badly but better than a tie.

**Corrected answer to #38**, replacing what I wrote an hour ago:

> Within the region CBES estimates, `g_marginal` localises a held-out image reference
> significantly better than any convergence statistic (+0.111 AUC, p 0.001) and `g` is
> indistinguishable from MKDA (+0.001, p 0.980). But CBES estimates only 9% of the brain and
> 34% of the truth's top decile, and over the whole brain ALE is the best arm while every CBES map
> is the worst. So whether the machinery earns its keep depends on whether one counts the
> two-thirds of the signal CBES declines to estimate -- and for a reader who wants to know the
> extent of an effect rather than its peaks, that two-thirds is the answer they came for.

`coverage_radius` is the lever this identifies. It is currently `2 * fwhm`; a wider radius would
cover more of the truth at the cost of claiming estimates further from any evidence. The
trade-off is now measurable, which it was not before, and it deserves a sweep.

## The second collection owed to #38 is blocked by data, not by effort

Two things checked, and both are worth recording rather than quietly dropping.

**The NeuroVault "animal" collection cannot support a coordinate-based test at all.** Its 11
movie-watching z maps peak at:

```
  study   max |z|   cluster-forming cut is 3.09
   8836      3.89
   8838      0.70
   8854      1.50
   8891      0.56
   8892      1.58
   8893      3.24
   8894      1.23
   8895      1.22
   8956      1.36
   8962      3.45
   9000      2.62
```

Seven of eleven peak below 2. Under realistic reporting -- a p < 0.001 forming cut with a
family-wise extent test -- **every one of the eleven reports nothing**, which the extraction duly
returns. So this collection is not a weak test of a coordinate method; it is not a test of one,
because there would be no published table to meta-analyse. Worth knowing because
`second_collection.py` used a bare uncorrected `U = 3.2905` on the same data, which would have
yielded foci from only the three studies peaking above 3.29 -- so that earlier experiment was
running on three studies, not eleven.

**The harvested 258-study corpus is strong enough but is not grouped into usable collections.**
98% of its studies have a peak clearing the forming cut and 80% clear a plausible corrected
height, so strength is not the problem. Paradigm membership is:

```
  None / Other                43
  rest eyes open              10
  rest eyes closed             8
  go/no-go task                7
  episodic recall              6
  2nd-order rule acquisition    6
  lexical decision task         6
  monetary incentive delay      6
```

The two rest-state sets have no task contrast to meta-analyse. The largest genuine task paradigm
is seven studies, which split in half gives three supplying coordinates -- too thin to carry the
#38 comparison, whose effect sizes are on the order of 0.1 in AUC.

**So the honest position:** every comparative and magnitude claim about CBES rests on the 21-study
NIDM pain collection, because it is the only assembled collection with enough studies sharing a
paradigm whose maps are strong enough to yield coordinate tables under realistic reporting. That
is a limitation of the validation rather than of the method, it is not fixable by more analysis of
what is on disk, and it should be stated plainly wherever those claims appear. Curating a second
collection -- grouping the 43 unlabelled studies by contrast name, or harvesting a new paradigm
with 15+ studies -- is the work that would lift it, and it is data work rather than modelling.

## The default kernel width is what cost CBES the convergence comparison

Sweeping the right parameter. NIDM pain, 8 splits, held-out image reference, `g_marginal` scored.

```
                brain covered  top decile covered  AUC on covered  AUC whole mask
fwhm  6 mm             0.013               0.060           0.677           0.527
fwhm 10 mm             0.090               0.336           0.757           0.642   <- default
fwhm 16 mm             0.263               0.636           0.798           0.745
fwhm 24 mm             0.608               0.879           0.780           0.782
MKDA density           1.000               1.000              --           0.667
ALE (earlier run)      1.000               1.000              --           0.753
```

**The whole-mask deficit was a kernel-width artefact, not a property of the method.** At the
default 10 mm, CBES covers 9% of the brain, captures a third of the truth's top decile, and loses
to MKDA over the whole mask (0.642 against 0.667). At 16 mm it covers 26% and 64%, and *beats*
MKDA on the whole mask (0.745). At 24 mm it covers 61% and 88%, and beats ALE too (0.782 against
0.753). So the earlier finding -- "over the whole brain ALE is the best arm and every CBES map is
the worst" -- holds only at the default kernel.

**And it is not a trade of precision for extent.** The `AUC on covered` column improves as well,
0.757 at 10 mm to 0.798 at 16 mm. A wider kernel makes the estimate better *where it already
existed* and also extends it. Only at 24 mm does the covered-AUC turn down slightly (0.780),
which is where the trade finally appears. On this collection the sweet spot for both columns is
around 16 mm and the default is well below it.

**The qualification that matters, and it is a real one.** The optimum depends on how smooth the
reference is. The truth here is an inverse-variance pooling of z maps -- a smooth, spatially
extensive field -- and a wide kernel is rewarded for matching its extent. Against a punctate
truth a narrow kernel would win. So the honest claim is not "16 mm is right" but:

> On a reference of the kind an image-based meta-analysis produces, CBES's default 10 mm kernel is
> narrow enough to cost it both accuracy and the comparison against convergence estimators, and
> widening it to 16 mm improves both the estimate and its extent. Whether 16 mm is right in
> general depends on the extent of the effects being pooled, which no single collection settles.

For reference, the convergence estimators' own kernels are in this range: ALE uses an
N-dependent Gaussian of roughly 10-12 mm FWHM and MKDA a 10-15 mm sphere. A 16 mm kernel is on
the wide side of the literature's practice but not outside it.

**The silence radius, swept separately since it changes the estimates and not their extent:**

```
radius  8 mm    AUC on covered 0.688
radius 20 mm    AUC on covered 0.757   <- default
radius 30 mm    AUC on covered 0.767
```

Monotone and still rising at 30 mm, consistent with the docstring's own note that `prevalence`
increases with this radius at every true value. The default is close to the best of these but not
at it; the gain from 20 to 30 mm is 0.010, against 0.041 from widening the kernel from 10 to
16 mm. So the kernel is the parameter that matters and the silence radius is second order.

## The kernel width trades the map against the interval

Directly on the uncertainty question, and it complicates the kernel recommendation. Coordinates
only, 12 studies, truth 0.800, prevalence 1:

```
                        bias  mean se  sd of g  se/sd  coverage  half/truth
fwhm 10 mm (default)  +0.233    0.163    0.041   3.99      1.00        0.40
fwhm 16 mm            +0.232    0.109    0.038   2.85      0.17        0.27
fwhm 24 mm            +0.245    0.085    0.033   2.56      0.00        0.21
```

**Widening the kernel does not touch the bias and shrinks the interval, so coverage collapses.**
The bias is flat to within noise because it is a property of the reported heights, which the
kernel does not alter. But a wider kernel pools more foci into each voxel, so the standard error
falls by half, and an unchanged bias against a halved interval takes coverage from 1.00 to 0.00.

Same signature as adding more studies -- a fixed bias with a shrinking interval -- and it means
the two deliverables want opposite settings:

| what you are using | what the kernel should be |
| --- | --- |
| the map: where the effect is and how extensive | **wider**; 16 mm beat 10 mm on both accuracy and extent on real data |
| the interval: whether a voxel's magnitude is pinned | **narrower**; only the wide interval at 10 mm covers, and it covers by being wide |

Neither is good on its own. The 10 mm interval covers at 1.00 against a nominal 0.95 with `se/sd`
of 3.99 -- not calibrated, merely wide enough to contain a +0.23 bias, which is the
metrics-that-lie pattern again. Widening the kernel strips away that accidental protection and
exposes the bias, which is arguably the more honest state even though the coverage number looks
worse.

So the recommendation on task #56 has to be split rather than stated flatly: a wider kernel is
better for the map and worse for the interval, and the reason the interval looks better at 10 mm
is that a bias it cannot see is hidden by a width it did not earn. What would actually fix the
interval is removing the bias, which needs images -- and at 6 of 12 images the bias is +0.013 and
coverage 0.99 at the default kernel already.

## A gap in the coverage table I should have caught: it never used the recommended configuration

Every arm of the coverage table ran with `peak_bias=None`. The docstring's recommendation for a
collection with both kinds of study is `peak_bias="per-study"` with `peak_bias_scale="images"`,
and `do_coordinates_help.py` even carries a note from an earlier session saying that judging
mixing without it "tests a configuration the code tells you not to use".

The two settings work by different mechanisms and that is why it matters:

- **`peak_bias=None` (what I measured).** Images enter as extra unbiased contributions alongside
  the coordinates. The coordinate values are untouched, so the bias falls only because the
  coordinate channel's share of the pooling weight falls. That is exactly the weight-share model,
  `bias(f) = b0 (1-f) / ((1-f) + r f)`, and it fitted to a maximum residual of 0.022 -- because it
  is the right model *for this configuration*.
- **`peak_bias="per-study"`, `peak_bias_scale="images"` (the recommendation).** The images are
  additionally used to read off the scale constant that the coordinate arm's own values are then
  divided by. The coordinate contributions are *corrected* rather than merely outvoted.

So the headline I have been reporting -- that the coordinate channel is "diluted, never
corrected" -- may be a property of the configuration I chose rather than of the method. If the
calibrated setting corrects it, then at a low image fraction it should do much better than
dilution predicts, which is precisely the regime that matters: two images among twenty studies is
where the weight-share model says +0.16 of bias survives.

Recorded before the numbers: if calibration works, the calibrated arms should beat their
`peak_bias=None` counterparts at 2 of 12 images by well more than the seed noise, and the
weight-share fit should not describe them. If they come out the same, the scale calibration is not
doing anything the dilution does not already do, and the "diluted, never corrected" framing stands
for both settings.

Either way this is a gap in how I set the bed up, not a discovery: the recommended configuration
was documented and I did not use it. The habit it argues for is to run the configuration the
docstring recommends *first*, and treat anything else as the variant.

## Note for a future session: the `benchmark` check is noisy

`benchmark` failed on commit 546b48a with:

```
| +  | 53.3±10ms  | 66.5±0.5ms | 1.25 | bench_cbma.TimeCBMA.time_mkdachi2_studyset |
PERFORMANCE DECREASED.
```

It was noise and a single re-run came back green. Three things made that the right diagnosis
rather than a guess, and they are the checks to repeat rather than re-derive:

1. **The benchmark is of `MKDAChi2`, which the PR cannot reach.** The only shared file the branch
   touches is `nimare/meta/utils.py`, and that diff is purely additive -- four new functions and
   one import reformatted across lines, with no existing function modified. `MKDAChi2` imports
   from `cbma.base`, `cbma.utils` and `meta.kernel`, none of whose used functions changed.
2. **The base measurement's own error bar spans the difference.** 53.3 ± 10 ms is a 19%
   coefficient of variation; 53.3 + 10 = 63.3 against an "after" of 66.5 ± 0.5. The PR-side
   measurement is the tight one.
3. **No other benchmark moved**, across about 96 of them.

Also worth knowing: the suite contains no CBES benchmark, so the PR's own new code is not timed
at all. Adding one would be a genuine improvement and is deliberately not being done here -- it
widens the PR beyond what was asked -- but it is the reason a benchmark failure on this PR is
*a priori* unlikely to be about the PR.

### And the prediction resolves against me: the coordinate channel IS corrected

Truth 0.800, 12 studies, 40 replications, the same bed:

```
configuration                                   mean g    bias  mean se  cover
2 of 12: peak_bias=None (what the table used)    0.904  +0.104    0.141   0.93
2 of 12: per-study, scale read off the images    0.753  -0.047    0.109   1.00
```

At the image fraction that matters -- two of twelve, `f = 0.17` -- the recommended configuration
cuts the bias from **+0.104 to -0.047**, less than half the magnitude and now slightly *under* the
truth rather than over it. So the coordinate values really are being corrected, not merely
outvoted.

**This retracts a headline I reported twice.** "The coordinate channel is diluted, never
corrected" is false of the estimator; it is true only of `peak_bias=None`. And the weight-share
model that fitted so beautifully --

```
    bias(f) = b0 (1 - f) / ((1 - f) + r f),   b0 = 0.260, r = 6.40, max residual 0.022
```

-- is a correct and well-specified model of the configuration the documentation tells users not
to rely on. Its excellent fit was never evidence that dilution is the mechanism available; it was
evidence that dilution is the mechanism *when you switch the correction off*, which is what I had
done in all twelve arms.

The projections that hung off it go too: "2 images of 20 studies leaves the magnitude 20% high"
was computed from `b0` and `r`, and on this evidence the calibrated setting would do materially
better. They need re-deriving from calibrated arms before being quoted anywhere.

What survives is narrower and still useful: with the correction switched off, dilution is the
whole story and the image *share* rather than the donor count governs it. That is a real finding
about the mechanism of `peak_bias=None`, and it explains why `peak_bias='per-study'` alone -- with
`peak_bias_scale` left at its default 1.0 -- moved nothing (+0.259 to +0.255): the per-study
factors correct the part of the bias that varies between studies, and without a scale read off
images there is nothing to fix the common part, exactly as the docstring says.

**The lesson, which is the same one as this morning in a new place.** I characterised a method for
a whole session in a configuration its own documentation warns against, having read that warning
earlier in the session and written it into the notes. The habit is not "read the docs" -- I had --
it is: **make the documented default the first arm of every comparison, and label any other
setting as the variant it is.**

### The complete calibrated-versus-dilution comparison

```
truth 0.800, 12 studies, 40 replications      mean g    bias  mean se  cover
2 of 12: peak_bias=None                        0.904  +0.104    0.141   0.93
2 of 12: per-study, scale read off the images   0.753  -0.047    0.109   1.00
6 of 12: peak_bias=None                        0.811  +0.011    0.092   1.00
6 of 12: per-study, scale read off the images   0.759  -0.041    0.084   1.00
```

The full picture is better than the one I drew from the 2-image rows alone, and better for the
estimator.

**The calibration pins the scale to about -5% of the truth, independent of the donor count.**
-0.047 at two donors and -0.041 at six: a small, stable over-correction that does not care how
many images there are. That is the property that matters for a real collection, where the image
share is whatever the literature happened to provide and is usually small.

**Dilution's accuracy depends entirely on the share, so it wins only where you scarcely need the
coordinates.** +0.104 at two images against +0.011 at six. At six of twelve it is the more
accurate of the two -- but a collection with half its studies imaged is one where an IBMA on the
images is available anyway.

So the two settings are not ranked; they have different failure modes, and the recommended one
has the better-behaved failure. Stated for the docstring rather than for me:

> With the scale read off image donors, the magnitude comes back within about 5% of the truth
> whether two studies or six supply images, erring slightly low. With `peak_bias=None` it depends
> on the share of studies imaged -- 13% high at two of twelve, 1% high at six -- so the
> configuration that needs fewest images is also the one that is insensitive to how many there
> are.

Coverage is 1.00 in three of the four cells and 0.93 in the fourth, all against a nominal 0.95,
so the intervals are conservative throughout -- and at these biases they are conservative for the
right reason rather than by being wide enough to hide an error, which is what the coordinates-only
arms were doing.

Correction to my own retraction, written an hour ago: I said the calibrated setting "more than
halves the bias", which is true at two donors and false at six, where it is worse. The accurate
statement is that it removes the *dependence on the donor count* at the cost of a small constant
over-correction.

## The definitive coverage table: documented configuration primary

The bed was restructured so the configuration the docstring recommends is the primary arm set and
`peak_bias=None` is kept alongside on every imaged arm as a labelled variant. Before reading the
new numbers, the shared arms reproduce the old ones exactly, which is the check that the
restructure changed nothing but the selection:

```
arm                                  old table   new table
12 studies,  0 images (per-study)       +0.255      +0.255
12 studies,  2 images, peak_bias=None   +0.127      +0.127
```

So the two runs are directly comparable and the difference between the calibrated and uncalibrated
arms is the configuration rather than any drift in the bed.

```
arm                                   mean g    bias  mean se  sd of g  se/sd  cover  half/truth
12 studies,  0 images                  1.055  +0.255    0.171    0.080   2.14   0.75       0.42
12 studies,  0 images @ fwhm 16        1.048  +0.248    0.109    0.064   1.70   0.28       0.27
12 studies,  0 images @ fwhm 24        1.049  +0.249    0.086    0.057   1.50   0.05       0.21
12 studies,  2 images, calibrated      0.762  -0.038    0.107    0.081   1.32   0.99       0.26
12 studies,  2 images, peak_bias=None  0.927  +0.127    0.134    0.098   1.36   0.87       0.33
```

(remaining arms still running)

The 2-image pair is the one that matters, because two donors among a dozen studies is the regime a
real collection sits in. The calibrated arm is **four times closer to the truth** (0.038 against
0.127) on a smaller interval (0.107 against 0.134) with better coverage (0.99 against 0.87). There
is no axis on which the uncalibrated setting is preferable there, which is what makes it the right
primary and makes a whole session's worth of numbers taken in the other setting the variant.

### A crash found by running the documented configuration

Re-running the coverage table with `peak_bias_scale="images"` primary hit an unhandled error on
the first all-donor arm:

```
12 studies, 12 images, calibrated  ->  ValueError: No study contributed any in-mask voxels.
```

Reproduced minimally: 2 of 6 and 5 of 6 studies imaged calibrate fine (scale 0.677 and 0.476),
6 of 6 raises. The cause is structural rather than incidental. The scale is read off the donors by
comparing, at shared voxels, what the coordinate-only fit says against what each donor's image
says; when every study is a donor the roster less the donors is empty, so `_statistic` gets an
empty focus table and no image studies, and `_accumulate` raises from three frames down with a
message that names neither the cause nor the configuration.

Fixed (commit c8590c0) by returning 1.0 with `scale_source_ = "unset"`, which is the right answer
rather than a fallback: a donor's own peaks are dropped in favour of its image, so with every
study imaged there are no coordinate values for a scale to act on.

Worth noting how it was found. Two deliberate tests already exercise the all-donor collection
shape, and the coverage work had already identified it as the best-behaved configuration
measured -- but nothing had ever combined that shape with the documented `peak_bias_scale`
setting, because every arm of the table used `peak_bias=None`. The bug had been sitting behind
exactly the configuration gap that made the magnitude numbers describe a variant. Running the
documented default first would have found it hours earlier, which is the second time today that
habit would have paid.

### The all-donor fix verified behaviourally

The two all-donor arms now agree to every printed digit across 100 replications:

```
12 studies, 12 images, calibrated       0.782  -0.018    0.065    0.059   1.10   0.94
12 studies, 12 images, peak_bias=None   0.782  -0.018    0.065    0.059   1.10   0.94
```

That is the check the fix needed rather than merely a plausible-looking number. With every study
supplying an image the calibration has nothing to act on -- each donor's own peaks are dropped in
favour of its image, so no coordinate value exists for a scale to multiply -- and the fix returns
1.0 with `scale_source_ = "unset"`. A scale of 1.0 makes the `per-study` path arithmetically
identical to `peak_bias=None`, so the two rows *must* coincide, and they do, exactly. If they had
merely been close it would have meant the fallback was doing something.

It also confirms the arm is still a usable reference: bias -0.018, `se/sd` 1.10, coverage 0.94
against a nominal 0.95, which is the row that shows the pooling and the observed-information
standard error are correct. The fix did not cost that.

### And this retracts "the interval has no validated operating point"

I wrote earlier, in the success-criteria reflection, that "the interval never reaches nominal
coverage in any genuinely coordinate-based configuration" and therefore "the interval on `g` has
no validated operating point". The 12-study calibrated arms say otherwise:

```
arm                                   bias  mean se  se/sd  coverage
12 studies,  2 images, calibrated   -0.038    0.107   1.32      0.99
12 studies,  6 images, calibrated   -0.026    0.083   1.26      0.98
12 studies, 12 images (no coords)   -0.018    0.065   1.10      0.94
12 studies,  0 images               +0.255    0.171   2.14      0.75
```

The 2-image row is genuinely coordinate-based: ten of twelve studies speak only through
coordinate tables, and the two donors are supplying a scale rather than carrying the estimate.
It covers 0.99 against a nominal 0.95 with a bias of -0.038 on a truth of 0.800 -- under 5%.

So the corrected claim, which is narrower than "it works" and much narrower than what I said:

> With the documented configuration and at least two image donors, the interval on `g` covers --
> 0.94 to 0.99 against a nominal 0.95 across image fractions from 0.17 to 1.00 -- erring
> conservative, with `se/sd` from 1.10 to 1.32. Coordinates-only it does not: 0.75 at twelve
> studies and 0.35 at twenty-four, and worse as studies accumulate.

Two honest qualifications on the positive half. The coverage is *conservative* rather than
calibrated: `se/sd` of 1.32 means the interval is a third wider than the estimator's own
variability, so it covers partly by being generous. And it is one bed, one truth value, one
reporting regime.

But the shape of the claim changes. "No validated operating point" said the quantity was beyond
rescue; what is actually true is that the operating point requires the documented configuration
and a couple of donors, and I had been measuring the configuration without them. That is the third
conclusion this session that was a property of my setup rather than of the estimator, and the
count is itself the finding: when a result is negative, the first question should be whether the
harness is configured the way a user would configure it.

### The all-image floor, and two hypotheses eliminated

Every calibrated row of the coverage table is read against the all-image arm, which sits at
-0.018 with twelve studies and -0.022 with twenty-four on a truth of 0.800. Before quoting the
calibrated rows as "the coordinate channel's residual" it was worth knowing what that floor is.

**It is not the bed's own conversion.** The donor image is built from the study's t through
`peak_stat_to_hedges_g`, the same conversion the estimator uses, so a Jensen term was the first
suspect. Measured without CBES in the loop at all -- 400 draws at each of n = 20, 30, 40, 60,
reading the donor's own g at the read-out voxel:

```
n= 20  mean donor g +0.8031  bias +0.0031
n= 30  mean donor g +0.8034  bias +0.0034
n= 40  mean donor g +0.8034  bias +0.0034
n= 60  mean donor g +0.8032  bias +0.0032
```

Unbiased to +0.003, and *upward*. So the estimator introduces the shortfall.

**It is not the spatial kernel.** The read-out voxel is a local maximum of the truth, so any
spatial averaging pulls it down, which made peak attenuation the obvious candidate. Swept over
four kernel widths on the all-image arm:

```
kernel   mean g     bias  mean se      sd  cover
     4    0.762   -0.038    0.064   0.047   0.96
     6    0.762   -0.038    0.064   0.047   0.96
    10    0.762   -0.038    0.064   0.047   0.96
    16    0.762   -0.038    0.064   0.047   0.96
```

Identical to four decimals at every width, which is the point: with every study donating an
image the kernel has no foci to spread, so it is inert and the hypothesis was untestable this
way rather than merely wrong. This is the second time today the "confirm the parameter moved the
thing you are attributing to it" rule has caught a story before it was written down -- the first
was sweeping `coverage_radius` and getting 0.089 at every radius from 8 to 45 mm.

It also corrects a reading I had already made. At 24 replications the arm gives -0.038 and I
briefly took that as "narrowing the kernel makes it worse than the -0.018 in the table". The
table's -0.018 is 100 replications; 0.047/sqrt(24) = 0.0096, so the two differ by 2 standard
errors of sampling and nothing else. Comparing a 24-rep number against a 100-rep number is the
same class of error as comparing configurations, just cheaper to make.

**What is left is the weighting.** Hedges' variance is `1/n + g^2 / (2(n-1))`, a function of the
*observed* effect, so a voxel that drew high gets a larger variance and less inverse-variance
weight than one that drew low -- a downward bias by construction. The order is right: at n = 30
and g = 0.8 the g^2 term is 0.011 against 1/n = 0.033, so a third of the weight varies with the
square of a noisy quantity. `experiments/image_floor.py` substitutes a draw-independent variance
(`1/n` everywhere) and changes nothing else. If that removes the floor the mechanism is not a
property of this bed: real `g_var` maps carry the g^2 term, so any image pooling that uses them
inherits it.

### The weight-share model transfers out of sample, and says an image is worth three to five studies

The dilution model `bias(f) = b0 (1-f) / ((1-f) + r f)` was fitted on the twelve-study rows and
then evaluated at the twenty-four-study image fractions, which the fit never saw. Both channels
are measured against the all-image arm, since that floor is what an image-only fit lands on.

```
peak_bias=None   b0 = +0.273 (measured from the coordinates-only arm, not fitted), r = 4.57
   2 of 12  (fitted)   f=0.167   measured +0.127   predicted +0.125   err +0.002
   6 of 12  (fitted)   f=0.500   measured +0.027   predicted +0.031   err -0.004
  12 of 12  (fitted)   f=1.000   measured -0.018   predicted -0.018   err +0.000
   2 of 24  (HELD OUT) f=0.083   measured +0.163   predicted +0.171   err -0.008
   6 of 24  (HELD OUT) f=0.250   measured +0.073   predicted +0.086   err -0.013

calibrated       b0 = -0.032 (fitted; no coordinates-only arm exists), r = 3.00
   2 of 24  (HELD OUT) f=0.083   measured -0.064   predicted -0.047   err -0.017
   6 of 24  (HELD OUT) f=0.250   measured -0.049   predicted -0.038   err -0.011
  24 of 24  (HELD OUT) f=1.000   measured -0.022   predicted -0.022   err +0.000
```

Two things worth keeping.

**`r` is between 3 and 4.6, and it is the same across study count.** One study that shares its
map carries three to five studies' worth of pooling weight. That is the number to quote when
someone asks what an image is worth in a mixed collection, and it is far more useful than an
image *fraction*, because it converts: two images among twenty-two coordinate tables carry about
as much weight as nine coordinate tables, which is why 2 of 24 still moves the estimate a long
way.

**The calibrated channel's own residual is `b0 = -0.032`, against +0.273 uncorrected.** So
reading the scale off the images removes about 88% of the coordinate channel's bias and
overshoots slightly. That is the honest one-line summary of what `peak_bias_scale="images"`
buys, and it is a statement about the *correction* rather than about dilution.

I first read the held-out misses as the calibrated configuration behaving qualitatively
differently, and wrote a task saying it overshoots "nearly twice" what dilution predicts, blaming
the voxel set the scale is read over. Holding the image fraction fixed instead shows most of it
is not image-related at all:

```
                             12 studies  24 studies   shift
coordinates only, fwhm 10        +0.255      +0.246  -0.009
coordinates only, fwhm 16        +0.248      +0.243  -0.005
coordinates only, fwhm 24        +0.249      +0.244  -0.005
all images, calibrated           -0.018      -0.022  -0.004
```

About -0.005 of uniform downward drift on doubling the studies, with no images and no
calibration in play. After crediting dilution and that drift, roughly -0.008 is left that is
specific to the calibration -- real, but a few standard errors at 100 replications, not a
different regime. The lesson is narrow and repeatable: when a model misses out of sample, check
whether the misses are specific to the mechanism you are about to blame, by finding a row where
that mechanism is switched off.

### Coverage is nothing but the bias-to-width ratio, and that is bad news about large collections

Sixteen arms of the calibrated table span coverage from 0.00 to 0.99. All of it is predicted by
a shifted normal -- the estimate sitting `b` away from the truth with spread `sd`, judged against
a half-width of `1.96 se`:

```
coverage = Phi((1.96 se - b) / sd) - Phi((-1.96 se - b) / sd)

                           b/se  pred   actual
12 st,  0 img              1.49  0.84     0.75
12 st,  0 img fwhm 24      2.90  0.08     0.05
12 st,  2 img cal         -0.36  0.98     0.99
12 st, 12 img cal         -0.28  0.96     0.94
24 st,  0 img              2.08  0.40     0.35
24 st,  2 img cal         -0.83  0.94     0.91
24 st,  2 img None         1.60  0.70     0.58
24 st, 24 img cal         -0.49  0.94     0.97

mean absolute error 0.034, maximum 0.124 over all sixteen arms
```

There is no residual pathology to explain. The interval fails exactly when, and exactly as much
as, the bias-to-width ratio says it should. Every earlier attempt to characterise the interval as
"too narrow" or "conservative" was describing `b/se` in words.

The consequence is not comfortable. `se` falls roughly as `1/sqrt(studies)` and the bias does not
fall at all, so **coverage decreases as a collection grows**, in every configuration measured:

```
coordinates only       0.75 at 12 studies  ->  0.35 at 24
2 images, calibrated   0.99 at 12 studies  ->  0.91 at 24
6 images, calibrated   0.98 at 12 studies  ->  0.92 at 24
```

So the usual reassurance is inverted here, and the calibrated configuration's good coverage may
be a *small-collection* property rather than an operating point. The residual bias the
calibration leaves is about 0.03 on a truth of 0.800 -- some 4% -- and the interval is honest
only while `se` stays comfortably larger than that. Extrapolating the fitted model, 48 studies at
the same weight share as the 2-of-24 arm lands near 0.85 and 96 near 0.75.

That is an extrapolation, so it is being tested rather than asserted:
`experiments/does_more_studies_hurt.py` holds the image fraction at 1/12 and grows the collection
through 12, 24 and 48 studies. If coverage tracks the prediction, what the PR can claim about the
interval is conditional on collection size and has to say so. If coverage holds up instead, the
normal-shift model breaks between 24 and 48 studies -- which matters just as much, because it is
the model every other coverage statement here is being read through.

The reframing this suggests for #57 is worth stating separately. The guidance should not be about
an image *fraction* at all: with `r` between 3 and 4.6, what matters is the weight share
`r*n_img / (r*n_img + n_coord)`, and that falls when coordinate-only studies are added. Two
images among ten coordinate tables hold about 29% of the weight; the same two among twenty-two
hold 15%. **Adding coordinate-only studies to a mixed collection makes the magnitude estimate
worse**, not better -- it dilutes the only thing correcting the peak-height inflation -- while
simultaneously narrowing the interval. Both effects push coverage the same way.

### The all-image floor is inverse-variance weighting, and it is not CBES's to fix

The floor is the `g^2` term in Hedges' variance. Substituting a draw-independent variance (`1/n`
for every voxel of a donor, nothing else changed, paired on seeds) removes it:

```
studies  variance   mean g     bias     sem   mean se      sd  cover
     12    hedges    0.774   -0.026  0.0072    0.065   0.056   0.95
     12      flat    0.798   -0.002  0.0078    0.062   0.060   0.92
     24    hedges    0.776   -0.024  0.0051    0.045   0.039   0.98
     24      flat    0.801   +0.001  0.0053    0.043   0.041   0.97
```

Both draw-independent rows are unbiased to within 0.3 standard errors; both Hedges rows sit at
-0.024 to -0.026 and are identical across study count, which is the signature of a bias in each
weight rather than a small-sample artefact.

And it reproduces with no brain, no estimator and no images -- 200,000 replications of scalar
inverse-variance pooling, `g_i ~ N(mu, 1/n_i)` with `n_i ~ U(20,40)`:

```
12 studies, mu 0.4:  Hedges weights -0.0111  (-2.8%)   draw-independent -0.0002
12 studies, mu 0.8:  Hedges weights -0.0183  (-2.3%)   draw-independent -0.0000
12 studies, mu 1.2:  Hedges weights -0.0214  (-1.8%)   draw-independent -0.0000
24 studies, mu 0.8:  Hedges weights -0.0190  (-2.4%)   draw-independent +0.0001
```

The mechanism is textbook: `var(g) = 1/n + g^2 / (2(n-1))` is computed from the *observed* effect,
so a study that drew high gets a larger variance, less weight, and the pooled estimate is pulled
toward zero. It is a bias in each weight rather than a small-sample artefact, so it does **not**
shrink with more studies -- -0.0183 at twelve and -0.0190 at twenty-four.

`_accumulate` weights every contribution by `a = weights / var_g`, so this is not confined to the
image channel: `peak_stat_to_hedges_g` computes `var_g = bias^2 * ((n1+n2)/(n1 n2) + d^2/(2(n1+n2)))`
from the observed peak statistic too. For coordinates the variance is not merely a weight -- it is
the sigma^2 the censored likelihood is written against -- so the same substitution is not available
there without a two-step scheme, and reported peaks have a hugely inflated observed `d`, which
makes the size of the distortion in that channel a separate question rather than the same one.

**CBES cannot correct this, and that is the useful conclusion.** A donor's `g_var` map is an
opaque input. A map produced by converting a t map carries the `g^2` term; a map produced by a
proper mixed model at each voxel does not. The estimator has no way to tell which it was handed,
so substituting a pooled estimate into a variance formula it only assumes would be wrong as often
as right. The fix belongs where `g_var` is created, or in the user's knowledge of what they
supplied.

Practical consequence, and the reason this is worth recording rather than filing: about 2-3% of
the all-image arm's downward bias is this, not CBES. It also means the -0.032 residual attributed
above to the calibrated coordinate channel is measured against a floor that has a known and
separable cause, so the corrected channel's own residual is nearer -0.01 than -0.032 once the
weighting bias is removed from both sides.

### Coverage erosion with collection size, tested at a size never otherwise measured

The prediction from the bias-to-width model was that coverage falls as a collection grows, since
`se` shrinks as roughly `1/sqrt(studies)` and the bias does not shrink at all. The table's own
evidence for it was mixed: the coordinates-only rows are clean (share fixed at zero, bias constant
at +0.25, coverage 0.75 -> 0.35), but the two-donor rows are not, because between them the donors'
share of the pooling weight *halved*, so dilution is mixed in with the study count.

Grown with the weight share held fixed at one donor per twelve studies, 40 replications per row:

```
studies  images   mean g     bias  mean se      sd  se/sd  cover  predicted
     12       1    0.746   -0.054    0.120   0.095   1.27   0.97       0.97
     24       2    0.738   -0.062    0.081   0.056   1.44   0.97       0.96
     48       4    0.725   -0.075    0.054   0.035   1.53   0.72       0.81
```

The erosion is real and this is the first time it has been seen at forty-eight studies. But it
takes hold much later than the twelve-to-twenty-four step suggested -- flat across that doubling,
then falling sharply -- for a reason worth keeping: **`se/sd` rises with study count** (1.27, 1.44,
1.53), so the reported interval shrinks more slowly than the estimator's actual spread and partly
offsets the constant bias. The model's own prediction at forty-eight is 0.81 against 0.72
measured, about 1.3 standard errors at this replication count.

Two smaller things fall out. The bias *grows* slowly with study count at fixed share (-0.054,
-0.062, -0.075), which the dilution model does not explain and which is the same residual drift as
#58. And the extrapolation I wrote down before running this -- "48 studies near 0.85 and 96 near
0.75" -- was roughly right at 48 by luck: it assumed `se/sd` constant, which it is not, and got a
similar answer because the bias also grew.

So the claim survives its own test and the mechanism I first gave for it was incomplete. That is
the second time today that stating a claim, naming the measurement that would break it, and
running that measurement changed the claim rather than confirming it.

### The relative map has no external truth, and my test of it was the voxel-set mistake again

`g_relative = g / P95(|g|)` is built so that a common multiplicative inflation cancels exactly. If
the peak-height bias were one constant `c`, then `g = c*mu` and the ratio is free of it. That made
it the obvious candidate for an interval that covers coordinates-only, where the interval on `g`
covers 0.75 and 0.35. Measured against `mu / P95(|mu|)`, with `se_relative = se / P95(|g|)`:

```
arm                        truth     mean     bias    se/sd  cover
12 studies, 0 images       4.572    0.872   -3.700     2.11   0.00
24 studies, 0 images       4.572    0.893   -3.679     1.43   0.00
12 studies, 2 images       4.572    2.868   -1.704     1.29   0.03
```

**This is my test being wrong, not the estimator.** The fitted `P95(|g|)` came out 1.212
coordinates-only against the truth's `P95(|mu|)` of 0.175 -- a factor of seven. The truth field is
a handful of Gaussian blobs in a 30^3 volume and therefore mostly zero, so its 95th percentile
sits in the blobs' skirts; the fitted map is kernel-smoothed, censored and inflated over its whole
covered support, so its 95th percentile sits somewhere else entirely. Two quantiles of
differently *shaped* distributions do not cancel a common factor, whatever the factor is.

That is exactly the failure recorded as task #38 -- "it depends entirely on the voxel set" -- and I
made it in a script whose own docstring warns against it and then picks the wrong voxel set. The
honest conclusion is narrower and more useful than either a pass or a fail here:

> **`g_relative` has no external truth to be scored against.** Its normalizer is a property of the
> estimator's own map, so "does `g_relative` cover" is not a well-posed question without also
> fixing what the truth's normalizer is computed over -- and any choice of that is a choice about
> the answer.

One thing the arms do say cleanly: `se/sd` for the relative map (2.11, 1.43, 1.29) matches `se/sd`
for `g` in the same fits (2.14, 2.01, 1.32). Dividing by the normalizer scales the estimate and its
error together, so relative precision is unchanged -- which is what `se_relative = se / P95(|g|)`
assumes and is the one part of the proposal that survives.

The question itself is still worth answering, and it does not need a normalizer.
`experiments/does_the_scale_cancel.py` asks it as a **ratio between two voxels**: two sites of the
bed differ only in magnitude, 0.800 and 0.600, so the truth is 1.333 and no voxel set enters. If
`g(v1)/g(v2)` covers coordinates-only, the unidentified scale really is the whole problem and a
ratio-valued output is worth having; if it does not, "coordinates identify the pattern but not the
scale" is too generous a summary.

### The whole coverage table scored the wrong interval

The class docstring tells callers to build the interval from `se` against a **t on `dof`**, not
against a normal, and gives a reason: a normal interval was measured at 73-94% of nominal against
91-97% for the t. Every arm of the coverage table -- all three versions of it -- scored
`g +/- 1.96 se`.

That is not a rounding difference. `dof` is `n_eff - 1` with `n_eff` a Kish effective sample size
over the studies reaching a voxel, so it is far below the study count: measured directly on a
twelve-study fit with two image donors it is **4.46**, where the t critical value is 2.666 against
1.960. The interval the documentation recommends is **36% wider** than the one I measured.

Propagating each arm's own bias and spread through a t, without refitting anything:

```
arm                   measured  model@1.96  |   t,dof=5  t,dof=11  t,dof=23
12 st,  0 img             0.75        0.84  |      0.99      0.94      0.89
12 st,  2 img cal         0.99        0.98  |      1.00      0.99      0.99
24 st,  0 img             0.35        0.40  |      0.83      0.59      0.49
24 st,  2 img cal         0.91        0.94  |      0.99      0.97      0.95
24 st, 24 img cal         0.97        0.94  |      0.99      0.97      0.96
```

At a dof near 5, which is what the one direct measurement suggests, the coordinates-only arms go
from 0.75 and 0.35 to something near 0.99 and 0.83. **That would overturn the table's central
conclusion** -- that the interval is usable with donors and not without them -- and it is the
conclusion I had just finished writing into the docstring.

So the emphasis was misplaced in a specific way. The bias measurements stand: +0.255 coordinates
only against -0.038 calibrated is a fact about the point estimate and no interval changes it. What
does not stand unqualified is "coordinates alone leave a bias no interval width can absorb",
because the interval the documentation actually recommends is a third wider than the one that
sentence was measured against.

The bed now returns `dof` and reports coverage under both intervals side by side, `cov(z)` and
`cov(t)`, with the critical value taken per replication rather than at the mean dof since it is
nonlinear in it. Re-running. The docstring carries the caveat in the meantime rather than waiting
for the table.

The general lesson is sharper than "check your critical value". The bed was built to measure
whether the documented interval covers, and it measured a *different* interval than the one
documented three paragraphs above the table it produced -- the same class of error as running
every arm at `peak_bias=None` when the docstring recommends `peak_bias_scale="images"`. Both times
the fix is the habit already written down: **make the documented configuration the first thing you
measure**, and that includes the documented way of reading the output, not only the documented
settings.

### The max-statistic guard is well-targeted and slightly too permissive

Measured on the cell that exposed the familywise defect and on a control, 20 studies, global null,
N 10-1000, 200 permutations, 60 simulations each:

```
foci/study  guard fired  voxel FWE  rejected|kept  n kept
         2        0.783      0.050          0.231      13
         6        0.017      0.050          0.051      59
```

Two readings, and they point the same way.

**The two-foci cell is mitigated but not controlled.** The overall 0.050 -- against 0.150 to 0.180
before the guard -- is not error control: a withheld fit cannot reject, so the aggregate is just
`(1 - 0.783) * 0.231`. Among the fits the guard *passes*, which is the population a user is in when
they get a p-value at all, rejection is 0.231 against nominal 0.05, exact binomial p = 0.0245 with
a 95% lower bound of 0.066. The mechanism is right and the threshold is in the wrong place.

**The six-foci control is clean, and that is what makes tightening cheap.** The guard fires on 1
fit in 60, and the rate is nominal both overall (0.050) and among the kept (0.051, n = 59). So the
guard is not refusing legitimate collections, and there is a lot of headroom: a threshold strict
enough to catch the marginal two-foci fits would cost almost nothing here.

What remains is *which* threshold. The guard refuses only when a null is both sparse (under 20
distinct attained maxima of 200) and narrow (coefficient of variation under 0.05). Distinct counts
run about six throughout the two-foci regime, so the thirteen kept fits are almost certainly
sparse-but-not-narrow -- kept on the docstring's reasoning that "a coarse but wide distribution
still separates the observed value from the bulk". `experiments/which_limb_lets_them_through.py`
cross-tabulates both statistics against rejection so the answer comes from the joint distribution
rather than from tuning a constant to a thirteen-fit measurement.

### The excess interval width is the censoring term. The tau2 story was wrong.

`se/sd` runs 1.1 to 2.1 across the coverage table, and under the documented `t` interval that is
what makes every default-kernel arm cover regardless of configuration. I had recorded
DerSimonian-Laird truncation as the leading cause -- a truncated estimator of a quantity whose
true value is zero has a positive mean, so the fit would be charging the interval for
heterogeneity that is not there -- and written down the number it had to hit: fitted `tau2` near
0.12 when the truth is 0.

A two-fit smoke test settled it, on a coordinates-only twelve-study fit:

```
configuration        g       se   fitted tau2
baseline          1.093    0.189      0.0000
tau2 none         1.093    0.189      0.0000
selection none    1.093    0.100      0.0000
```

Fitted `tau2` is **exactly zero**, so there is no truncation bias to speak of, and switching DL
off changes the `se` by nothing at all. Target 0.12, value 0.0000. The hypothesis is dead, and it
died for three seconds of compute rather than the eighty-replication run it was written for --
which is the habit working, not a wasted prediction.

What moves the `se` is the censoring term. `selection_model="none"` nearly halves it, 0.189 to
0.100, which against the coordinates-only spread of 0.080 takes `se/sd` from about 2.4 to about
1.25. Most of the excess, though not all.

So the excess width is the observed information of the censored mixture, which carries uncertainty
about which component each observation came from. That is genuine uncertainty about a latent
quantity rather than an arithmetic error, and two readings remain open with different
consequences:

* **The `se` is right and `sd` is the wrong comparison.** The observed information is the error for
  the parameter mu *in the selection model*; the replication-to-replication spread of the point
  estimate is a different quantity. On this reading nothing needs fixing and the lesson is only
  that coverage must never be read as validating the uncertainty.
* **The information is overstated.** Under a correctly specified model the two should agree
  asymptotically, so a factor of 2.4 at twelve studies is evidence of misspecification or of the
  profiling step double-counting.

One piece of evidence already bears on this and points somewhere useful. The oracle bed -- the
brute-force MLE against the same censored mixture with *known* variances -- covered 94.5% to
98.4%, close to nominal. So the information is about right where the model is true, which makes
the excess in this bed more likely to be misspecification of the **reporting process** than of the
likelihood. That is a different and more interesting defect than a variance bug, and it is
testable: the bed's reporting is a cluster-forming cut with a fixed assumed extent, while the
estimator's censoring term assumes height thresholding.

### Why `se/sd > 1` is anomalous in one direction only

I had left two readings of `se/sd` open: that the observed information is the right error for `mu`
in the selection model and the replication spread is simply a different quantity, or that the
information is overstated. A variance decomposition settles which way the anomaly points.

The censored likelihood conditions on the **positions** of the reported foci and models their
magnitudes; it does not model where the peaks fell. So `se` estimates a variance conditional on
the design -- which studies contributed, at which distances. The replication `sd` is marginal over
designs, because every replication draws fresh studies, fresh sample sizes, fresh noise and
therefore fresh focus positions. By the usual decomposition,

```
sd^2  =  E[ Var(g_hat | design) ]  +  Var( E[g_hat | design] )   >=   E[ Var(g_hat | design) ]
```

so a well-calibrated conditional `se` should come out **smaller** than the marginal `sd`, by
however much the design varies. `se/sd` ought to sit at or below 1.

It is 1.1 to 2.4. The anomaly is therefore not a yardstick mismatch that could excuse either
direction -- the mismatch that does exist pushes the ratio the other way, so the measured excess
understates how overstated the `se` is. That is evidence for the second reading, and it is
independent of anything about `tau2` or the critical value.

Which leaves a sharper question than "is the se right". The oracle bed covers 94.5% to 98.4% with
known variances, so the information is about right where the model holds. The excess appears when
the reporting process is realistic. `experiments/is_it_the_reporting_scheme.py` changes only how
foci are extracted -- cluster-forming cut versus voxelwise FDR versus Bonferroni -- while holding
the estimator fixed, which is the one comparison that separates a bad model of reporting from a
bad likelihood.

### A second, separable reason the interval is too wide: `dof` counts fewer studies than `se` uses

`n_eff` is Kish's `(sum w)^2 / sum w^2` over the kernel weights, and **only studies whose foci
reach a voxel contribute a `w`**. The censored likelihood also uses the studies that reported
nothing nearby -- their silence bounds `mu`, which is the entire point of the selection model. So
the `se` draws on more studies than the `dof` counts, and referring that `se` to a `t` on
`n_eff - 1` charges the interval for a smaller sample than it actually used.

Measured on a twelve-study coordinates-only fit:

```
             at the read-out voxel   median over covered voxels
n_studies                   10.00                         2.00
n_eff                        4.76                         1.14
dof                          3.76                         1.00
```

At the read-out voxel the likelihood is using ten studies (reports plus silences) and the
interval is referred to 3.76 degrees of freedom. `t(0.975, 3.76) = 2.87` against
`t(0.975, 9) = 2.26` -- **27% of extra width from the `dof` choice alone**, before anything about
the `se` itself. Over the map as a whole it is far worse: a median `dof` of 1.00 gives a critical
value of 12.71.

So the over-wide interval has two separable causes, and they compound:

1. the `se` is about twice the estimator's own spread, located in the censoring term
   (`selection_model="none"` nearly halves it; `tau2` is exactly zero and irrelevant);
2. the `dof` counts only the studies that spoke, while the `se` uses the studies that stayed
   silent too.

The second is an internal inconsistency rather than a modelling question -- the same quantity is
being credited with information in the numerator and denied it in the denominator. But what the
right `dof` is for a censored-likelihood observed information is a real statistical question, not
a typo to patch: candidates are the censoring roster `n_studies`, a Kish count weighting silent
studies by their censoring information, or abandoning the `t` reference altogether in favour of a
profile-likelihood interval, which needs no `dof` at all. That last is the principled answer and
the most work.

Recording it rather than changing it, because it moves every interval the estimator reports.

### The roster `dof` is not a fix, and that closes the calibration question

Referring the `se` to the censoring roster instead of `n_eff` is the obvious remedy for the
inconsistency above, and it is worth about 15-20% of width. It does not change which arms cover.
Projecting the roster critical value onto each arm's own bias and spread:

```
arm                              dof(n_eff)  cov(t)  dof(roster)  cov(roster)
12 st, 0 img, fwhm 10                   4.5    0.98          9.0         0.95
12 st, 0 img, fwhm 16                   7.8    0.52          9.0         0.49
12 st, 0 img, fwhm 24                   9.7    0.13          9.0         0.17
12 st, 2 img calibrated                 4.4    1.00          9.0         0.99
12 st, 2 img peak_bias=None             4.4    0.99          9.0         0.96
24 st, 0 img, fwhm 10                   9.1    0.56         18.9         0.51
24 st, 0 img, fwhm 24                  20.2    0.00         18.9         0.00
24 st, 2 img calibrated                 8.6    0.96         18.9         0.96
24 st, 2 img peak_bias=None             8.6    0.74         18.9         0.77
```

The same spread, 0.00 to 0.99, either way. **The bias decides which arms cover; the `dof` decides
only how wide they are.** So the `dof` question is a consistency defect worth fixing on its own
terms and not a route to a calibrated interval.

And there is nothing unexplained left in the calibration. The bias-to-width model predicts
`cov(t)` across all sixteen arms to a mean absolute error of **0.021** -- better than the 0.034 it
achieves on `cov(z)`, so the median `dof` is a sufficient summary of the per-replication variation
in the critical value. Everything the coverage column does is accounted for by two numbers, the
bias and the width, and the only open questions are why the `se` is about twice the estimator's
spread (#60) and what to do about the bias.

That is worth stating plainly because it retires a line of enquiry. Three separate hypotheses
about the interval have now been tested and closed -- `tau2` truncation (falsified, fitted `tau2`
is exactly zero), the critical value (real, 36% of width, does not change which arms cover), and
the `dof` reference (real, 15-20% of width, does not change which arms cover). The interval is not
mis-scaled. It is wide because the `se` is large, and it misses because the point estimate is
biased, and those are the only two things left.

### What the interval problems do *not* touch: the p-values

Worth checking rather than assuming, because it bounds how much the `se` findings cost. `z` is
`g / se` (effectsize.py:3107) and `p` compares `|z|` against a permutation null of `|z|`, with
each permutation re-fitting and so getting its own `g` and its own `se` (3438-3443).

A permutation p-value is valid for whatever statistic it is computed on, provided the statistic is
computed the same way under permutation. It does not require the `se` to be calibrated. So an
`se` that is twice the estimator's spread inflates the observed `z` and every null `z` alike, and
the significance map is unaffected by everything in the three sections above. The docstring
already says this and says it correctly -- "the p-values are unaffected either way, they come from
the permutation null, not from referring `z` to any distribution".

The separation is worth keeping in view when reading the rest of these notes. The uncertainty
findings bear on the **interval on `g`**, which is wide and sometimes misses. They do not bear on
which voxels are called significant, except through the max-statistic guard, which is a separate
defect with its own measurements. A user who reads the thresholded map and ignores `se` is not
affected by any of it.

### A prediction to be wrong about: the reporting-scheme test may explain the *bias*, not just the width

`is_it_the_reporting_scheme.py` was written to ask whether the excess interval width is
misspecification of the reporting process. It also prints the bias per scheme, and on reflection
that column may be the more important one. Writing the prediction down first.

`peak_bias="per-study"` computes each study's `rho_k` from the random-field peak-height
distribution at that study's own threshold and sample size. That is an *absolute* correction, not
a relative one -- so on the face of it the scale should already be pinned, and the docstring's
claim that the common scale is unidentified needs a reason. Measured, `per-study` moves the bias
from +0.259 to +0.255: it does essentially nothing.

The reason, I think, is that a reported focus is not the object the theory describes. The theory
is about a local maximum of a smooth statistic field exceeding a **height** threshold. What the
bed reports -- and what papers report -- is the maximum within a *cluster* that survived a
forming cut and an extent criterion, in a field whose smoothness the estimator does not know. So
`rho_k` is wrong by a factor that is common across studies because it depends on the things not
modelled (the smoothness, the cluster criterion), which is exactly the signature of "an
unidentified common scale".

If that is right, then under **voxelwise FDR or Bonferroni** -- where a reported peak really is a
local maximum clearing a height threshold -- the theoretical `rho_k` should be close to correct
and the coordinates-only bias should fall well below +0.255, without any image donor. And the
consequence would be large: the peak-height inflation would not be irreducible, it would be a
consequence of assuming height thresholding when papers report clusters, and the fix would be to
model cluster reporting rather than to require images.

If instead the bias stays near +0.255 under every scheme, the inflation is not about the reporting
event and the unidentified scale is genuinely unidentified from coordinates -- which is what the
docstring currently claims and what four measured-and-rejected remedies already support.

Stating it in advance because the last three predictions recorded this way -- the weight-share
model transferring, coverage eroding with collection size, `tau2` truncation explaining the
width -- came out two right and one wrong, and the wrong one was worth more than either of the
others.
