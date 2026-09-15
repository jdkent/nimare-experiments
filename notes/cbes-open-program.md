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
| `identification_by_power_spread` (E3) | re-running |
| `marginal_cancellation` (E7) | re-running |

Unaffected, and why:

- `threshold_sensitivity`, `calibrate_coverage_bed`, `heights_on_real_data` -- written after the
  fix, or on real data where the estimator's assumption is the correct one.
- `reporting_probability` and the occupancy/window-of-detectability work -- the exact detection
  model never passes through the estimator's conversion at all.
- Every earlier real-data result (held-out HCP, NIDM pain split-half, the NeuroVault paradigm
  sets) -- real collections report real t or z maps, so there is no mismatch to make.

**Still to audit, not yet done:** the false-positive-rate beds (`mixed_null`, `mixed_null2`,
`mixed_null3`, `few_images_null`, `approx_null_images`, `config_matrix`, `sim_validate`,
`threshold_correction_fit`). These measure error rates rather than magnitudes, and a permutation
p-value is invariant to a monotone per-study transform of the values -- but the conversion is
*not* a common transform, it depends on each study's `n`, so a roster with heterogeneous sample
sizes could in principle see its null shifted. The error rates were established in earlier
sessions and are load-bearing for the PR, so this needs checking rather than assuming. Logged as
a task.

The general point for the protocol: after finding an instrumentation bug, **go back over every
conclusion the instrument produced, including the ones that were negative.** A bug that inflates
an estimate also inflates the evidence against predictions that said the estimate would be
smaller, and those retractions need revisiting too. I found this one by accident rather than by
audit, which is the wrong way round.
