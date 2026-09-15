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
