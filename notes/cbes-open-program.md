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

### That prediction was mis-specified, and the smoke test says something better

Three fits, before the full run:

```
cluster, infers u    g 1.093   bias +0.293   foci/study  3.58
cluster, told u      g 1.049   bias +0.249   foci/study  3.58
fdr,     told u      g 1.040   bias +0.240   foci/study 47.92
```

**The prediction is wrong, and it was wrong because the test does not vary what I said it
varied.** `scheme` sets the multiplicity correction and therefore the *threshold*; the selection
event is set by `focus`, and `focus="max"` takes the maximum statistic within each surviving
cluster under *every* scheme. So FDR reports far more clusters at a lower cut, but each reported
focus is still a spatial maximum. I claimed the test would show whether the reporting event
matches the theory, and it cannot: the event is the same in all three arms.

What the arms do say is worth more than what I predicted. Going from 3.58 to **47.92 foci per
study** -- a thirteenfold change in how much of the map is tabulated, at a much lower threshold --
moves the bias from +0.249 to +0.240. **Table density and threshold level barely touch it.** The
inflation is not about how selectively peaks are reported.

Which leaves spatial maximum selection as the whole of it, and that is consistent: the maximum of
a smooth field over a blob exceeds the value at the blob's centre, whatever cut admitted the blob.
It also explains why `peak_bias="per-study"` does nothing. Its `rho_k` comes from the random-field
peak-height distribution, which depends on the field's **smoothness** -- and the estimator is
never told the smoothness, so `rho_k` is wrong by a factor common to every study. That is exactly
the shape of "an unidentified common scale".

Telling the estimator the true threshold is worth +0.293 to +0.249, about 0.044 -- real, and the
one part of the story `threshold="study-min"` is responsible for. It bears on #53 and it is much
smaller than the 0.20 of prevalence accuracy that decision turns on.

The constructive line this opens: **papers report smoothness.** If the residual is a
smoothness-dependent factor then a smoothness metadata field could pin the scale with no image
donor at all, which is the thing four other remedies failed to do. Testable by sweeping the bed's
own smoothness and checking whether the bias tracks it the way random-field theory says it should.
Recorded as a direction rather than a result.

### The complete table, and the one row that makes the point

All 21 arms, scored against both intervals. Raw log in
`results/coverage_both_intervals_21arm.log`.

```
arm                                     bias  se/sd   dof cov(z) cov(t)  width_t
12 studies,  0 images                 +0.255   2.14   4.5   0.75   0.98     0.57
12 studies,  0 images @ fwhm 16       +0.248   1.70   7.8   0.28   0.52     0.32
12 studies,  0 images @ fwhm 24       +0.249   1.50   9.7   0.05   0.13     0.24
12 studies,  2 images, calibrated     -0.038   1.32   4.4   0.99   1.00     0.36
12 studies,  2 images, peak_bias=None +0.127   1.36   4.4   0.87   0.99     0.45
12 studies,  6 images, calibrated     -0.026   1.26   7.0   0.98   1.00     0.25
12 studies,  6 images, peak_bias=None +0.027   1.27   7.0   0.99   1.00     0.27
12 studies, 12 images, calibrated     -0.018   1.10  11.0   0.94   0.97     0.18
24 studies,  0 images                 +0.246   2.01   9.1   0.35   0.56     0.33
24 studies,  0 images @ fwhm 16       +0.243   1.79  16.2   0.03   0.04     0.20
24 studies,  0 images @ fwhm 24       +0.244   1.59  20.2   0.00   0.00     0.16
24 studies,  2 images, calibrated     -0.064   1.35   8.6   0.91   0.96     0.22
24 studies,  2 images, peak_bias=None +0.163   1.48   8.6   0.58   0.74     0.29
24 studies,  6 images, calibrated     -0.049   1.41  10.8   0.92   0.96     0.19
24 studies,  6 images, peak_bias=None +0.073   1.38  10.8   0.89   0.93     0.22
24 studies, 24 images, calibrated     -0.022   1.11  23.0   0.97   0.97     0.12
12 studies,  0 images, tau 0.3        +0.268   1.70   3.6   0.84   0.94     0.91
12 studies,  6 images, tau 0.3        -0.048   1.06   6.7   0.92   0.94     0.43
12 studies, 12 images, tau 0.3        -0.015   1.12  11.0   0.93   0.96     0.30
```

**The row to read first is the second from last.** Coordinates-only with real heterogeneity
covers **0.94** -- nominal, to the eye -- with a half-width of **0.91 of the effect**. On a truth
of 0.800 the interval runs from about 0.34 to 1.80. It covers because it admits almost any
magnitude. Coverage alone cannot tell that row from `24 studies, 24 images` at 0.97 and a width
of 0.12, and no amount of care about the critical value or the degrees of freedom changes which
of those two a user should trust.

That is the whole lesson of the last two hours, and it cost three superseded versions of this
table to learn: the first had a mis-calibrated reporting regime, the second scored the wrong
`peak_bias` configuration, the third scored the wrong interval. Each time the fix was to measure
what the documentation describes -- the regime papers actually report in, the settings the
docstring recommends, and finally the interval the docstring tells callers to build. And each
time the correction changed a conclusion I had already written down.

The bias-to-width model predicts `cov(t)` across all twenty-one arms to a mean absolute error of
**0.024**, against 0.036 for `cov(z)`, so the median `dof` is a sufficient summary and there is no
behaviour left unaccounted for. Two numbers describe the entire coverage column.

### The over-wide `se` is masking the bias, and under heterogeneity two errors cancel

Sixty replications, twelve studies, three configurations so the two candidate causes separate:

```
donors  true tau  configuration  fitted tau2  mean se      sd  se/sd     bias  cover
     0      0.00       baseline       0.0039    0.176   0.077   2.28   +0.260   0.83
     0      0.00      tau2 none       0.0000    0.171   0.074   2.32   +0.258   0.80
     0      0.00 selection none       0.0039    0.109   0.081   1.35   +0.272   0.17
     0      0.30       baseline       0.0146    0.261   0.156   1.67   +0.250   0.90
     0      0.30      tau2 none       0.0000    0.232   0.163   1.43   +0.260   0.80
     0      0.30 selection none       0.0146    0.128   0.137   0.93   +0.333   0.20
     6      0.00       baseline       0.0054    0.084   0.063   1.33   -0.031   1.00
     6      0.00      tau2 none       0.0000    0.082   0.063   1.29   -0.034   1.00
     6      0.00 selection none       0.0058    0.075   0.067   1.12   -0.042   0.93
```

**`tau2` is not the cause, and at real heterogeneity it is underestimated sixfold.** Fitted
`tau2` is 0.0146 against a true 0.0900 at `tau = 0.3`. My truncation hypothesis had the sign
backwards: DerSimonian-Laird here is not spuriously positive, it is badly low.

**The censoring term is the cause, and removing it lands `se/sd` at or below 1.** The `se` drops
38% coordinates-only, 51% under heterogeneity, 11% with six donors, and `se/sd` goes 2.28 → 1.35,
1.67 → **0.93**, 1.33 → 1.12. The 0.93 is the telling one: below 1 is where the variance
decomposition says a well-calibrated *conditional* `se` should sit, so the censoring term is the
anomaly and the rest of the machinery is about right.

**But the inflated `se` is what holds the coverage up.** Switch the censoring term off and
coordinates-only coverage collapses from 0.83 to **0.17**, and from 0.90 to 0.20 under
heterogeneity, because the bias is unchanged (+0.260 → +0.272, +0.250 → +0.333) while the interval
narrows by half. That is not a reason to keep the inflation. It is the opposite: **an honest
interval around a biased point estimate ought to miss**, and the over-wide `se` is masking a bias
the user should be told about. Fixing it would make the estimator report plainly that
coordinates-only `g` is not trustworthy, instead of covering by accident. Less flattering, more
honest.

**And under heterogeneity two errors cancel.** The censoring term inflates the `se`,
DerSimonian-Laird underestimates `tau2` sixfold, and the net is coverage of 0.90 at 1.96 and 0.94
at the documented `t`. So the near-nominal coverage in the heterogeneity arms is not evidence the
model is right; it is two errors of opposite sign, and either one fixed alone makes coverage
worse. That is worth knowing before anyone quotes the τ arms as reassurance -- I nearly did.

The three donor-plus-heterogeneity rows complete the picture and add one thing:

```
donors  true tau  configuration  fitted tau2  mean se      sd  se/sd     bias  cover
     6      0.30       baseline       0.0604    0.150   0.132   1.13   -0.063   0.92
     6      0.30      tau2 none       0.0000    0.109   0.140   0.78   -0.032   0.82
     6      0.30 selection none       0.0557    0.110   0.134   0.82   -0.044   0.88
```

**The `tau2` underestimation is largely a coordinates-only problem.** With six donors
DerSimonian-Laird recovers 0.0604 of a true 0.0900 -- 67%, against 16% coordinates-only. Images
carry direct variance information and the estimate improves accordingly, which is another thing
donors buy that was not on the list.

The `tau2 none` row at real heterogeneity is a sanity check the bed passes: ignoring genuine
heterogeneity gives `se/sd` 0.78 and coverage 0.82, too narrow in exactly the way it should be.

One residual: with the censoring term off, coordinates-only `se/sd` is still 1.35, while with
donors it is 1.12 and under heterogeneity 0.93 or 0.82. So the residual is specific to the
coordinates-only homogeneous configuration and shrinks when either the donors or the heterogeneity
change.

One candidate was ruled out by reading the code rather than simulating: `_apply_peak_bias` scales
`g` by `rho_k` and `var_g` by `rho_k**2`, consistently, so the correction does not put the estimate
and its variance on different scales.

What is left is that **the kernel weights are random and correlated with the values.** A
coordinate study contributes at the weight of whichever of its foci is nearest the voxel, and both
which focus that is and how far it lands are functions of that study's noise; the pooled variance
treats the weights as known constants. A study that drew a high peak close to the voxel contributes
a large value *at a large weight*, and the formula does not know the weight was chosen partly by
the noise that set the value. The dilution pattern fits: weight-1 image contributions dominate the
weighted sum (1.35 to 1.12), and between-study variance dominates the weighting noise (1.35 to
0.93). `experiments/are_the_weights_the_residual.py` places the foci at fixed positions so the
weights are identical in every replication and only the magnitudes vary -- a deliberate break with
realism, stated as such, for a mechanism this bed cannot otherwise isolate.

## Redesigning the estimator around silence

The instruction: rebuild around the silence of coordinates and remove the now-dead code
aggressively. What that means concretely, worked out before deleting anything.

**The premise.** Coordinate tables contribute *presence and absence only*. Magnitudes come from
images. Measured on three collections (NIDM pain, HCP MOTOR_LH, HCP EMOTION_FACES), that has the
lowest rmse in all three (paired p <= 0.0021) and the best-calibrated level in all three.

**What dies, and why it is genuinely dead rather than merely unused:**

* The **kernel** -- `fwhm`, `_kernel_support`, `_study_voxel_weights`, `_focus_geometry`,
  `kernel_min_weight`. Checked by grep: these are reached only from `_accumulate`, i.e. only to
  spread coordinate *magnitudes* over voxels. Silence geometry is `coverage_radius`, which stays.
  This retires #56 by deletion.
* All of **`peak_bias`** -- the factors, the scale, the calibration from donors, the
  `scale_source_`/`n_scale_donors_`/`scale_interval_` plumbing, `_scale_is_pinned`,
  `peak_stat_to_hedges_g`, `null_peak_overshoot`, `peak_information`, `null_peak_mean_g`. The
  entire winner's-curse correction exists to repair reported peak heights. No heights, no
  correction. This also retires the all-donor calibration crash and the two-donor gate.
* **Threshold inference** -- `infer_threshold_from_minimum`, `_expected_min_peak`, `"study-min"`,
  `"pooled-min"`. All three infer the cut from the smallest reported *statistic*. The threshold
  is still needed (silence is only informative against a cut) but must now be supplied or
  defaulted. This retires #53 by deletion.
* **The peak statistic itself** -- `stat_column`, `_resolve_stat_column`, `_reported_z`. A
  consequence worth naming: the estimator no longer needs a reported statistic per focus, so it
  can consume the tables most of the coordinate literature actually publishes.
* **`g_relative` and `g_absolute`** and `_relative_g`. They exist because the coordinate scale was
  unidentified. With magnitudes from images only, `g` is on the effect-size scale by
  construction and there is nothing to relativise.
* **`_permute_magnitudes`** -- there are no coordinate magnitudes to permute.

**What survives, checked rather than assumed.** The censoring machinery entire
(`_censoring_terms`, `_ReportingPairs`, `_SilentPairs`, `_coverage_entries`,
`_apply_selection_model`, `_working_sets`, `_update_prevalence`, `_fit_chunk`,
`_observed_information`), image loading and pooling, `tau2`, `prevalence`, `coverage_radius`.

**The null survives, which I had expected to be the blocker.** `_permute_image_values` reassigns
each image's values among its own voxels -- many states even with one image, unlike a sign flip --
and tests the same hypothesis the magnitude shuffle did: within a study, effect size is unrelated
to location. The coordinate silence pattern is *fixed* across permutations, so it contributes the
same structure to the observed statistic and to every null draw, and cancels. The relocation null
is not needed, which matters because it was already measured and rejected: not conditional on
multiplicity, and permuting whole rows broke the multiplicity invariant.

**One decision is forced rather than optional.** `n_eff` is Kish over the pooled weights, and with
coordinate magnitudes gone only images carry a weight, all equal to 1 -- so `n_eff = k` and
`dof = k - 1`, which is **0 at one image** and makes the documented interval `nan`. #62 stops
being a consistency question. Taking the roster: the likelihood genuinely uses every study's
report or silence, and with equal image weights the Kish count adds nothing over a plain count.

### A near-miss worth recording

`null_effect_variance` is the one surviving caller of the conversion being deleted, so its
`d = 0` value had to be inlined. I wrote the usual approximation, `J**2 / N` for one sample, and
checked it against the pre-surgery module rather than assuming: **0.04608 against 0.05150 at
N = 20, wrong by 12%.** `d_to_g` uses the *exact* variance,
`(N-1)(1 + N d**2) h**2 / (N (N-3)) - d**2`, which at `d = 0` is `(N-1) h**2 / (N (N-3))`, while
the two-sample branch does use the approximate form. The two designs do not share a formula and
assuming they did would have silently mis-scaled every silence in the model.

Both designs now match the old path to 0.000e+00. The habit that caught it is the one already in
`PROTOCOL.md`: name the check that would show the change is wrong, and run it before moving on.

## The redesign around the coordinate channel's silence -- and the limb that was missing

Shipped. `nimare/meta/cbma/effectsize.py` 3784 -> ~2600 lines. Removed: `fwhm`, the spatial
kernel, `peak_bias`, `peak_bias_scale`, `stat_column`, `kernel_min_weight`, `use_images`,
`g_relative`, `g_absolute`, `scale_interval_`, `peak_information_`, threshold inference
(`study-min`/`pooled-min`), `permute-magnitudes`, and ~20 helper functions. Images are now
**required**; a coordinates-only collection is refused rather than returned as zeros.

### Two defects found by measuring rather than by reading

**1. The cutoff was left on the z scale, which took the whole channel inert.** The z->g
conversion of each study's reporting threshold lived inside `_apply_peak_bias`, which the
surgery deleted. A z of 3.29 is about 18 sampling standard deviations at N = 30, so
`prob_silent_null` saturated at 1 and the censoring term said nothing at any `mu`. Caught
because a test asserting `threshold=2.5` differs from `threshold=5.0` returned values
identical to eight decimals. Reinstated as `reporting_cutoff_to_g`, whose docstring carries
the assumed-dof sensitivity table (3.30 -> 0.653 at df=29 against 0.604 at df=1000).

**2. Dropping the *report* limb of the indicator biases `mu` down past the threshold.** The
silent pairs were the only evidence about the reporting indicator, so the model read the
observed silence fraction against a denominator excluding every study that reported.

The three arms, field simulator, truth known exactly, 6 seeds, true g = 0.5 at the focus,
study cutoffs ~0.6 g, stratified because 9204 of 9261 voxels have truth < 0.05:

| arm | rmse truth<0.05 | rmse 0.05-0.25 | rmse truth>=0.25 | bias truth>=0.25 | g at focus |
|---|---|---|---|---|---|
| images only            | 0.120 | 0.095 | 0.112 | +0.038 | 0.544 |
| silence only           | 0.108 | 0.076 | 0.091 | -0.060 | 0.352 |
| **+ report at the named voxel** | 0.114 | 0.076 | **0.074** | **-0.042** | **0.444** |
| + report over the 20 mm sphere | 0.457 | 0.343 | 0.119 | +0.097 | 0.457 |

**The asymmetry is the finding.** A silence is a statement about a neighbourhood -- nothing
within the radius cleared the cut. A report is a statement about *one voxel*, because a
reported peak is a local maximum selected for being large and displaced from where the effect
is. Asserting it across the sphere is catastrophic (rmse 0.457). Asserting it at the named
voxel only is the design.

**A non-spatial testbed cannot see this.** A one-voxel grid-search likelihood said "restore
the limb" unconditionally -- mean mu 0.351 -> 0.524 for a true 0.500, rmse 0.192 -> 0.129 --
and has no radius to get wrong. Logged against the calibration rule in CLAUDE.md: the 0-D bed
cannot speak to anything that turns on the extent over which an event is asserted.

### Real data, NIDM pain, 21 studies, 8 paired splits, held-out half as truth

Coordinates extracted the way papers produce them (cluster-forming cut, whole clusters,
one focus per cluster).

| estimate | r | rank r | AUC | bias | at top | rmse |
|---|---|---|---|---|---|---|
| images only, pooled     | +0.621 | +0.484 | 0.893 | +0.136 | +0.137 | 0.269 |
| CBES, silence off       | +0.633 | +0.499 | 0.899 | +0.142 | +0.155 | 0.272 |
| **CBES `g`**            | +0.579 | +0.486 | 0.889 | **+0.075** | **-0.065** | **0.210** |
| CBES `g_marginal`       | +0.536 | +0.458 | 0.861 | -0.004 | -0.226 | 0.191 |

rmse -22% (paired p = 0.0003), bias -45% (p < 0.0001), and the +0.137 overestimate at the
strongest voxels becomes a slight under. rank r +0.002 (p = 0.90) and AUC -0.004 (p = 0.32)
are unmoved; Pearson r costs 0.042 (p = 0.017). **The channel corrects the level and leaves
the ordering alone.** That is the same split as before -- images win the pattern, the
censoring wins the level -- but the level is now much better rather than slightly worse.

### jdkent's design point, confirmed

"The shrinkage shouldn't be to 0, but underneath the estimated threshold." Correct, and it is
what the model does once the cutoff is on the right scale: `mu` is pushed down only as far as
the silences carry it, toward the cut rather than toward nothing. The pull toward zero comes
through `pi` -- `g_marginal` is the better map where nothing was reported (bias +0.056 against
+0.095 for the images alone) and the worse one where an effect exists (-0.096).

### Also shipped

`create_effect_size_coordinate_studyset(n_image_studies=k, image_dir=...)`. Without it there
is no way to simulate a valid CBES input at all, since images are now required. Writes
`g`/`g_var` for the first k studies from the same field draw that produced their peaks, so the
two channels are consistent. Requires `simulate_field=True`.

Test file rewritten from scratch: 105 tests -> 64, all passing, organised around what each
channel does rather than around the deleted parameters.

## jdkent's threshold idea, in its safe form: the minimum as a bound, not an estimate

Shipped as `clamp_threshold=True`.

Estimating the cut from the smallest reported value was measured and rejected -- it undoes an
order statistic it cannot identify and overshot a cluster-forming cut by about 1 z. But there is
a strictly weaker use: anything a study reported *cleared* its cut, so

    c_k <= min_j |z_kj|

is a hard inequality. So it can **clamp** an assumed constant rather than replace it. It can only
move a cutoff down, only for a study whose own table contradicts the assumption, and never below
the truth.

Three regimes, field simulator, 8 seeds, truth known exactly, assumption 3.29 z:

| regime | cutoffs moved | rmse quiet | rmse middle | rmse effect |
|---|---|---|---|---|
| studies at 2.4/2.8/3.29/3.8 z | 9 of 20 | 0.1236 -> **0.1094** (p=0.0001) | 0.0817 -> **0.0734** (p=0.015) | 0.0727 -> 0.0743 (p=0.79) |
| every study at 3.29 z | 0 of 20 | bit-identical | bit-identical | bit-identical |
| every study at 4.5 z (thin tables) | 0 of 20 | bit-identical | bit-identical | bit-identical |

A strict no-op in both regimes where it should not act, including the thin-table case (a paper
reporting only its strongest peaks, where the minimum sits far above the cut). Safe on by default.

**Unexplained, recorded rather than relied on.** In the first regime the clamp *beats the oracle*
(0.1094 against 0.1123 quiet, 0.0734 against 0.0749 middle), and in the third the true thresholds
are worse than the too-low assumption (0.1278 against 0.1059 quiet; bias at the effect +0.0154
against -0.0653). A cutoff slightly below the truth is compensating for something. Likely cause:
the model treats a report as `|g| >= c` while a reported peak is `|g| >= c` **and** a local
maximum -- a strictly smaller event -- so P(report) is overstated and a low c offsets it. Testable
by replacing the exceedance probability with the RFT peak-height survival; not done.

### Test-harness lesson, third this session

The "heights never reach the estimate" test scaled reported statistics by `3z + 7`. That is not
monotone in `|z|`: it pulls negative peaks toward zero and *lowers* the smallest reported
magnitude, which is now exactly the quantity the estimator reads. The test failed, the code was
right. Use a pure scaling when the thing under test is a function of `|z|`.

## The interval, re-measured after the fixes: the estimate improved and the interval got worse

Task #67. Every previous se/sd figure was taken while the censoring term was inert (the cutoff
was on the z scale), so none of them were measuring the coordinate channel at all.

Field simulator, truth known exactly, 30 replications per arm, `threshold="reporting_threshold"`:

| arm | bias quiet | se/sd quiet | bias effect | se/sd effect | cov(t) effect | width effect |
|---|---|---|---|---|---|---|
| 20 studies, 1 image | +0.114 | 2.71 | -0.081 | 1.55 | 0.98 | 0.99 |
| 20 studies, 2 images | +0.093 | 3.17 | -0.082 | 1.78 | 0.99 | 0.92 |
| 20 studies, 5 images | +0.076 | 3.67 | -0.039 | 1.00 | 1.00 | 0.59 |
| 20 studies, 20 images | +0.045 | 3.74 | -0.021 | 1.64 | 0.99 | 0.35 |
| 2 images, silence off | +0.101 | 2.05 | -0.034 | 1.29 | 1.00 | 5.71 |
| 2 images, tau 0.3 | +0.093 | 3.17 | -0.088 | 2.15 | 0.98 | 1.62 |

**se/sd is now 1.55 to 3.74, worse than the 1.1 to 2.1 recorded before.** And the excess is
*located*, not inferred: switching the silence off drops it from 1.78 to 1.29 where the effect
is. The same indicator that fixes the magnitude inflates the error. Coverage is 0.97 to 1.00
everywhere and therefore says nothing -- at two images it covers with a half-width of 0.92 of the
effect.

The `silence off` width (5.71) is not comparable: with no censoring roster `dof` falls back to the
Kish count over image weights, which at two images is 1, and t(1) has a critical value of 12.71.
That is a correct statement about two studies, not a wider interval for the same information.

Bias also drifts the other way across the strata: **positive where the truth is near zero, negative
where the effect is**, at every image count. The absolute-value floor explains the first; the
second is over-shrinkage that more images reduce (-0.081 at one image to -0.021 at twenty) but the
silence does not.

Next: #66, the peak-height survival. If P(report) is overstated because a reported peak is a local
maximum rather than any exceedance, that is a candidate for both the negative bias at the effect
and for why a too-low cutoff outperforms the true one.

## Two attempts at the residual -0.039 bias, both falsified. The shipped choice is a joint optimum.

Where the effect is largest, `g` comes back at -0.039 of bias. Two mechanisms proposed, both
tested, both dead -- and both in a way that *confirms* the shipped configuration.

**1. "A report is over-credited, because the model treats it as any exceedance."** A reported
peak is `|g| >= c` **and** a local maximum, a strictly smaller event whose conditional
probability falls with mu, so the true P(report | mu) should rise more slowly than the plain
exceedance. Probed by raising the report limb's probability to a power `alpha`, which flattens
its mu-sensitivity without changing its range:

| alpha | 0.1 | 0.3 | 0.5 | 0.7 | **1.0** | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|---|---|
| bias, effect | -0.061 | -0.055 | -0.049 | -0.045 | **-0.039** | -0.024 | -0.008 | +0.011 | +0.024 |
| rmse, effect | 0.091 | 0.083 | 0.077 | 0.073 | **0.070** | 0.076 | 0.098 | 0.134 | 0.170 |

Sign backwards from the prediction: down-weighting makes the bias *worse*, monotonically. And
`alpha = 1` -- the principled value, already shipped -- is the rmse optimum. `alpha ~ 2` zeroes
the bias and costs 40% of rmse, so the knob trades bias for variance and buys nothing.

**2. "Displacement discards the studies that detected the region."** A peak sits where the noise
helped, so at the true focus a study that plainly detected the region names a *neighbouring*
voxel, read as sign 0 rather than as a detection -- leaving an indicator sample enriched for
genuine failures. Remedy: two radii, 20 mm for the silence (a reporting extent) and something
small for the report (a localisation error).

| report radius | rmse quiet | rmse middle | rmse effect | bias effect |
|---|---|---|---|---|
| **named voxel** | **0.1129** | **0.0737** | **0.0697** | **-0.0392** |
| 4 mm | 0.1337 | 0.0803 | 0.1131 | +0.0939 |
| 6 mm | 0.1585 | 0.1636 | 0.1224 | +0.1023 |
| 8 mm | 0.1818 | 0.2598 | 0.1257 | +0.1056 |
| 12 mm | 0.2786 | 0.3736 | 0.1257 | +0.1056 |
| 20 mm | 0.4534 | 0.3753 | 0.1264 | +0.1061 |

**There is no middle ground: the penalty starts at the first ring.** The bias flips from -0.039
to +0.094 at 4 mm -- one ring overshoots by more than the original undershoot, because a 4 mm
sphere asserts the indicator at seven voxels instead of one. Paired p < 0.0001 in the quiet
stratum at every radius.

The two probes agree, and that is the useful part: the report limb is correctly weighted at
`alpha = 1` and correctly extended at 0 mm, so **the shipped configuration is the joint optimum
of a two-parameter family** rather than an arbitrary choice. The mechanism in (2) is real; the
likelihood is simply far more sensitive to asserting `|g| >= c` where it is false than to
discarding an observation.

So the -0.039 is not in the report limb. Remaining candidates, none tested: Hedges' variance as
an inverse-variance weight (a known 2-3% downward pull, too small on its own), tau2 estimated
about the naive mean and so biased low (measured to move `g` from 0.883 to 0.822 at a true 0.8
when alternated), and the stratum definition itself -- 7 voxels at truth >= 0.25.

## The compression is gone. That was the old design's headline failure.

Never measured for the redesign until now, and it is the single biggest improvement.

Three-bin stratification cannot answer this -- it mixes the slope with the floor, and the top bin
held 7 voxels. So: four well-separated foci at true g of 0.2, 0.4, 0.6, 0.8 inside *one* map,
8 collections, 20 studies, 2 images, regressing estimate on truth over the 198 signal voxels.

| estimate | slope | intercept | @0.2 | @0.4 | @0.6 | @0.8 | recovered range |
|---|---|---|---|---|---|---|---|
| images only | 0.872 | +0.044 | 0.242 | 0.446 | 0.586 | 0.764 | 3.2-fold |
| **CBES g** | 0.729 | +0.050 | 0.188 | 0.345 | 0.605 | **0.798** | **4.2-fold** |
| CBES g_marginal | 0.750 | +0.017 | 0.151 | 0.293 | 0.589 | 0.765 | 5.1-fold |

`g` recovers **4.2-fold for a true 4-fold**. The old design gave 1.2-fold for 11-fold. An unknown
overall scale would leave the ratio alone, so the old compression was a real defect; it is now
gone, and if anything slightly over-spread.

And the *shape* of the improvement is the right one. `g` beats images-only at the strong foci
(0.798 for a true 0.800, 0.605 for 0.600) and sits further below at the weak ones (0.188 for
0.200 against images-only's inflated 0.242). A focus whose true effect is 0.2 against a cutoff
near 0.6 g was reported mostly by luck and *should* be shrunk. But that is also the entire source
of the slope of 0.73, so **the weak end is where to expect over-correction, not the strong end**.

### This retracts the -0.039 "bias where the effect is largest"

The three-bin measurement that drove #66 binned truth >= 0.25 into one column spanning 0.25 to
0.50, averaging voxels whose estimate is slightly high with voxels whose estimate is low and
reporting the mixture as a bias. At the focus itself `g` is 0.798 for a true 0.800 -- essentially
unbiased. So the quantity two experiments were chasing was partly an artefact of the binning.
Both falsifications still stand on their own terms (alpha = 1 and 0 mm are still the joint
optimum), but the motivating number was weaker than stated.

Slope and intercept are the honest summary for this estimator. Recorded in the docstring, which
previously claimed the coordinate channel "corrects the bias rather than the compression" -- now
measurably wrong and replaced.

### Correction: the weak end is where the correction *works*, not where it over-corrects

I shipped "the weak end is where to expect over-correction" off the absolute values. The relative
errors say the opposite:

| truth | images | err% | CBES g | err% | g_marginal | err% |
|---|---|---|---|---|---|---|
| 0.2 | 0.242 | **+21.0** | 0.188 | **-6.0** | 0.151 | -24.5 |
| 0.4 | 0.446 | +11.5 | 0.345 | -13.8 | 0.293 | -26.8 |
| 0.6 | 0.586 | -2.3 | 0.605 | +0.8 | 0.589 | -1.8 |
| 0.8 | 0.764 | -4.5 | 0.798 | -0.3 | 0.765 | -4.4 |
| **mean \|err%\|** | | **9.8** | | **5.2** | | **14.4** |

At the *weakest* focus CBES is more accurate than the images by a factor of three (-6% against
+21%), which is precisely what the selection correction exists for: a focus at a true 0.2 against
a cutoff near 0.6 g is reported mostly by luck, so pooling only what got reported reads it as far
too strong. The only focus where `g` is worse is 0.4 (-13.8% against +11.5%), comparable in size
and opposite in sign.

**So the slope of 0.729 is not the foci being compressed.** It comes from the blob skirts, where
the truth runs 0.05 to 0.2 and both arms are dominated by the floor that reading a map as `|g|`
imposes. The per-focus column, not the slope, is what says what a peak's magnitude is worth.

`g_marginal` is the one map that *is* worse at the weak end (-25%, -27%), because `prevalence`
falls toward its floor exactly where few studies reported. Narrowed in the docstring to "prefer
it only when comparing against an image-based reference", where it is the matching estimand.

Third correction this session from reading absolute rather than relative error. Worth a habit:
for a quantity spanning a range, report relative error per stratum, never a pooled absolute one.

### Second correction, with an n this time: 24 seeds and standard errors

The per-focus table moved twice because I shipped it off 6 and 8 seeds, where the SE on a value
of 0.2 to 0.8 is 0.02 to 0.03 -- about 10% relative, which is the size of the effects I was
reporting. At 24 seeds:

| estimate | slope | @0.2 | @0.4 | @0.6 | @0.8 | range | mean \|err\| |
|---|---|---|---|---|---|---|---|
| images only | 0.911 | 0.230+-.025 (+15%) | 0.419+-.029 (+5%) | 0.635+-.027 (+6%) | 0.792+-.025 (-1%) | 3.44x | **6.7%** |
| **CBES g** | 0.759 | 0.194+-.019 (-3%) | **0.329+-.022 (-18%)** | 0.634+-.016 (+6%) | 0.825+-.017 (+3%) | **4.26x** | 7.5% |
| g_marginal | 0.769 | 0.151+-.018 (-25%) | 0.267+-.024 (-33%) | 0.576+-.015 (-4%) | 0.780+-.019 (-3%) | 5.18x | 16.1% |

**What survives:** the range. 4.26-fold for a true 4-fold against 3.44-fold for the images alone.
The compression that was the old design's headline failure is genuinely gone.

**What does not:** "5.2% against 9.8%". Mean absolute relative error is 7.5% for `g` against
**6.7%** for the images -- a wash, slightly the wrong way. The arms trade errors focus by focus
rather than one dominating, so a single averaged number was never going to be the summary.

**The one clear defect: -18% at the 0.4 focus, about 3 SE, where the images are +5%.** That is the
middle of the detection window (1 of 17 studies reported), so exactly where the likelihood is most
sensitive to the reporting model. Replicated at 24 seeds. Unexplained.

### The alpha probe, re-run in the regime the first one could not reach

The first sweep's middle stratum was *exactly* invariant to alpha (0.0737 at every value) -- a dead
test, not a null result: in a single-focus bed that stratum is the blob skirt, where nothing is
reported, so there were no `sign = -1` pairs to act on. Re-run on the four-foci bed, which does
span the window (0/18, 1/16, 4/12, 12/6 reported/silent at truths 0.2/0.4/0.6/0.8):

| alpha | g@0.2 | g@0.4 | g@0.6 | g@0.8 |
|---|---|---|---|---|
| **1.0** | 0.153 (-23%) | **0.357 (-11%)** | **0.606 (+1%)** | **0.803 (+0%)** |
| 0.7 | 0.153 (-23%) | 0.336 (-16%) | 0.585 (-3%) | 0.745 (-7%) |
| 0.5 | 0.153 (-23%) | 0.324 (-19%) | 0.596 (-1%) | 0.708 (-11%) |
| 0.3 | 0.153 (-23%) | 0.312 (-22%) | 0.554 (-8%) | 0.660 (-18%) |

alpha = 1 is best at every focus where the limb is live, so #66 is now closed on a live test rather
than a dead one. And `g@0.2` is exactly invariant at 0.153 -- which is the structural fact worth
keeping: **below the detection window the report limb does not exist**, nothing is reported and
silence is nearly flat in mu, so `g` there is the image estimate and nothing the coordinate
channel does can change it.

### Habit to keep

Three corrections this session came from a measurement, not a model: absolute error where relative
was wanted, a bin that mixed two behaviours, and an n too small for the effect being claimed. All
three produced confident wrong statements that survived until re-measured. Report n and SE beside
any number that goes into a docstring.

## The 0.4 deficit located, and #66 was closed on a probe that could not represent its hypothesis

Instrumented the four-foci bed to report, per focus, the observed reporting rate against the rate
the model computes at the true mu, and the mu that the observed count alone would imply. 16
collections.

| truth | reported | silent | observed rate | model's rate | ratio | fitted mu | mu from count alone |
|---|---|---|---|---|---|---|---|
| 0.2 | 0.12 | 17.31 | 0.007 | 0.007 | 1.00 | 0.211 (+6%) | -0.144 |
| 0.4 | 0.88 | 16.50 | **0.050** | **0.089** | **1.78** | 0.337 (-16%) | 0.161 |
| 0.6 | 6.50 | 10.94 | 0.373 | 0.402 | 1.08 | 0.618 (+3%) | 0.586 |
| 0.8 | 12.56 | 5.19 | 0.708 | 0.797 | 1.13 | 0.850 (+6%) | 0.754 |

**The defect is the report limb's functional form.** The model uses `P(|g| >= c | mu)`, but a
paper reports a voxel only if it cleared `c` **and** was a local maximum. So the model overstates
the reporting probability, by a factor that peaks at **1.78 in the middle of the window** and
falls to 1.08-1.13 above it and 1.00 below it (where nothing is reported at all). An overstated
P(report) against an under-observed count drags mu down, and it drags hardest exactly where the
factor is largest: the 0.4 focus.

**Two things I had wrong.**

*"The indicator channel dominates: 17 observations against 2 image values."* No. At 0.4 the count
alone implies mu = 0.161 and the images imply about 0.4; the fit lands at 0.337. The channels are
genuinely combined, and the images carry most of the weight. Below the window the count implies a
*nonsensical* mu = -0.144 and the fit is +6% accurate, which is the same fact in a starker form.

*"#66 falsified: the report limb is correctly specified."* The hypothesis -- that a reported peak
is a local maximum and not any exceedance, so P(report | mu) is overstated -- is now **confirmed
and quantified**. What I falsified was my *probe* for it: raising the limb's probability to a power
`alpha` rescales the log-probability by a constant factor in mu, whereas the real correction is
mu-dependent (1.00, 1.78, 1.08, 1.13 across the four truths, non-monotone). A one-parameter
uniform flattening cannot represent a non-monotone mu-dependence, so the sweep was never a test of
the hypothesis. Reopened with a quantified target: make `P(report | mu)` reproduce 0.007 / 0.050 /
0.373 / 0.708 rather than 0.007 / 0.089 / 0.402 / 0.797.

The radius sweeps stay falsified on their own terms -- those did test what they claimed.

## SDM-PSI vs the redesigned CBES: SDM wins the pattern by a wide margin. This is a design problem.

21-study NIDM pain. SDM's coefficient is the one already on disk -- 50 imputations over
coordinates built at U = 3.2905. CBES gets the *same* extraction: 2 of the studies as g/g_var
images and the other 19 as coordinate tables, because the redesign requires at least one image.
Reference: the inverse-variance mean of the 19 images CBES never saw, so CBES is not scored
against its own input; SDM is mildly flattered (2 of its 21 tables came from reference studies),
which is the conservative direction. 228,483 voxels.

| estimate | r | rank r | AUC | mag ratio | \|est\| top | \|ref\| top |
|---|---|---|---|---|---|---|
| **SDM-PSI coeff (21 tables)** | **+0.689** | **+0.627** | **0.877** | 0.41 | 0.210 | 0.462 |
| images only (2 images) | +0.407 | +0.283 | 0.751 | 0.63 | 0.311 | 0.462 |
| CBES g (2 img + 19 tab) | +0.358 | +0.287 | 0.755 | **0.62** | 0.314 | 0.462 |
| CBES g_marginal | +0.390 | +0.289 | 0.755 | 0.57 | 0.291 | 0.462 |

**SDM wins the pattern decisively** -- r +0.689 against +0.358, rank r +0.627 against +0.287,
AUC 0.877 against 0.755. Not close, and the 2-of-21 flattery cannot account for it.

**CBES wins the magnitude**: 0.62 of the reference at the top quartile against SDM's 0.41. SDM is
2.4x too low there, CBES 1.6x.

### What this says about the design

The redesign discarded the coordinate channel's *spatial* information along with its heights.
SDM's imputation uses each peak's location and height to build a whole per-study map, and that is
where its pattern advantage comes from. CBES's silence channel corrects the level and does
essentially nothing for the pattern -- on the split-half bed that showed up as rank r +0.002
(p = 0.90) and AUC -0.004 (p = 0.32), and here it shows up as what that costs.

Worse: **CBES's r (+0.358) is *below* images-only (+0.407).** The silence channel slightly hurts
the pattern, which matches the split-half result (r -0.042, paired p = 0.017). So with 2 images
the pattern is set by 2 studies' noise and the indicators do not rescue it.

The central open question for the design is now sharp: **can the coordinate *locations*
contribute to the pattern without their heights?** In principle the indicator carries location
information -- -1 where studies reported, +1 where all were silent -- but it enters only through a
censoring term that nudges a magnitude the images already set. SDM, with no images at all, builds
its pattern entirely from locations and gets nearly twice the correlation.

Two honest readings, and they are not exclusive:

1. CBES answers a different question (mu, the effect among studies that have one, on a real
   effect-size scale) and is better at it -- the magnitude column says so, and SDM's imputed
   coefficient is 2.4x low.
2. For anyone who wants a map of *where*, 21 coordinate tables through SDM beat 2 images plus 19
   silences through CBES, and it is not marginal.

Next: does the gap close as the image count rises? If CBES needs 5+ images to match SDM's
pattern, the honest recommendation is narrower than the docstring currently implies.

## Retraction: "SDM beats the redesign on pattern, and this is a design problem"

jdkent caught the reference. It is wrong, and wrong in SDM's favour.

CBES got 19 *coordinate tables*, and those tables were extracted from **the same 19 images that
formed the reference**. So "the 19 images CBES never saw" was false -- it saw thresholded
summaries of exactly those studies. SDM got 21 tables, including all 19 reference studies, and
builds its entire pattern by imputing effect mass at their peaks -- which sit at the reference's
own high-|g| voxels. So a large part of r = +0.689 is SDM recovering the reference's noise, while
CBES's magnitude came from 2 non-reference images and its tables entered only as weak indicators.
Not symmetric flattery. The comparison cannot be salvaged without re-running SDM on a subsample,
and there is no subset of studies the existing SDM run did not see.

One sub-conclusion survives, because it does not involve SDM: images-only (+0.407) was the single
arm with a clean reference and it still beat CBES g (+0.358).

### The design question, answered on a reference with no confound

`validate_redesign.py` already does it right: the work half supplies both the coordinates *and*
the images, the held-out half supplies the reference, and the two halves are disjoint studies.
21-study pain, 8 paired splits.

| estimate | r | rank r | AUC | bias | err at top | rmse |
|---|---|---|---|---|---|---|
| images only (pooled) | +0.621 | +0.484 | 0.893 | +0.136 | +0.137 | 0.269 |
| CBES, silence off | +0.633 | +0.499 | 0.899 | +0.142 | +0.155 | 0.272 |
| **shipped CBES g** | +0.579 | +0.486 | 0.889 | **+0.075** | **-0.065** | **0.210** |
| shipped CBES g_marginal | +0.536 | +0.458 | 0.861 | -0.004 | -0.226 | 0.191 |

**The coordinate channel buys a large improvement in level for a small loss in pattern.** rmse
-23% and bias -47% against silence-off, with the +0.155 overestimate at the strongest voxels
becoming -0.065. The cost is Pearson r (-0.054, paired p = 0.0026) and essentially nothing on the
scale-invariant metrics: rank r -0.013 (p = 0.03), AUC -0.010 (p = 0.05).

So "the redesign threw away the coordinate channel's spatial information" was reading a confound.
On a clean reference the pattern cost is small and the level gain is large. That is a defensible
trade and it is what the design was chosen for.

**Still open, and now correctly scoped:** whether a method that uses peak *locations* to build a
pattern (SDM) beats CBES on pattern when both are held to the same disjoint split. Needs SDM run
per split. Split-half rather than leave-one-out: one study's map is too noisy a reference to test
the magnitude claim, which is where CBES wins, and 8 SDM runs beat 21.

Blocker: `sdm_parse pp` works; `sdm_parse mean` exits 0 without producing `analysis_MyMean`, and
`MyMean=mean` exits 2. No pdftotext for the tutorial, and the argument names are not in the
binary's strings. The 50-imputation coefficient from the earlier full-collection run is intact.

## The `benchmark` CI check: two failures, both spurious, demonstrated not asserted

`bench_cbma.TimeCBMA.time_mkdachi2_studyset` was flagged twice -- 1.61x on `a011f45` and 1.54x on
`ecddf5c`. Neither commit can cause it: the first changes only `test_meta_effectsize.py`, the
second only a docstring. And the PR's changes to files that *are* on MKDAChi2's path are inert --
`meta/utils.py` is purely additive (new `_padded_flat_to_masked`, `_gpd_*`, `_max_statistic_maps`,
none called by MKDAChi2) and `studyset/requirements.py` adds one key (`"sum": np.sum`) to a dict
literal inside a cached branch.

Reproduced locally, same benchmark body, 8 samples after a warm-up:

| revision | median |
|---|---|
| f22be2c (PR base) | 48.1 ms |
| origin/main | 48.7 ms |
| PR head | **46.4 ms** |

All 46-49 ms, with the head slightly *faster*. CI claimed 42.8 -> 66.1 ms. The 66 ms exists only
in CI.

The mechanism is runner variance, and the CI logs show it on the **base** side of the same two
runs: `time_mkdadensity_dense` 506+-5ms vs 637+-3ms, `time_mkdachi2_dense` 797+-3ms vs 1.05+-0s,
`time_kda` 50.6 vs 56.8 ms. Swings of 25-30% on identical code, because asv builds and benchmarks
the base and the PR in separate virtualenvs at different moments on a shared runner. A 1.5x flag
on a 42 ms benchmark is within that.

`benchmark` then passed on the head with no change from me, which is the confirmation. No comment
posted on the PR: the check is green there, so there is nothing to stand down from -- but the
local numbers are recorded here in case it recurs.

## HCP, reference from held-out subjects: the correction does harm where prevalence is truly 1

The cleanest reference available. MOTOR_LH, 786 subjects: 480 cut into 16 synthetic studies of
30, coordinates extracted the way papers do (cluster-forming, whole clusters, one focus per
cluster, 8 mm apart, no cap -- 14.4 peaks per table); the remaining **306 subjects, used to make
no coordinate at all**, give the truth. The selection producing the peaks is therefore statistically
independent of the quantity compared against, which nothing before this was.

| estimate | r | rank r | AUC | mag ratio | \|est\| top | \|ref\| top |
|---|---|---|---|---|---|---|
| images only (2 images) | **+0.845** | +0.576 | **0.973** | **0.85** | 0.394 | 0.446 |
| CBES g | +0.830 | +0.575 | 0.967 | 0.63 | 0.290 | 0.446 |
| CBES g_marginal | +0.785 | +0.564 | 0.943 | 0.54 | 0.259 | 0.446 |

**Pooling the two images alone beats CBES on every metric, magnitude included.** That is the
opposite of the pain split-half, where CBES cut rmse 23% and bias 47%.

### Why, and it is not a bug

**HCP has prevalence exactly 1 by construction** -- every synthetic study draws from the same
population, so there is no between-study absence for the zero-inflated mixture to find. And the
images here are unthresholded maps of 30 subjects from that same population, so pooling two of
them is already nearly unbiased (0.85, the shortfall being Hedges' variance weighting and noise).
There is nothing for a selection correction to correct.

Measured directly: **fitted prevalence is 0.664 against a true 1.0**, rising to 0.929 at the
strongest decile. Where the effect is strong enough that most studies report it, pi is right;
where it is weaker, the model reads "failed to clear its threshold" as "has no effect". That is
the mixture's identifiability problem -- both explanations fit a silence -- and here every unit of
mass it puts on the null component is pure error.

Worse, it is not a clean pi/mu trade. If pi = 0.664 and the marginal is 0.85, then mu should read
1.28; it reads 0.63. The silences push mu down through the censoring term *and* pi down through
the mixture, so `g_marginal = pi * mu` is shrunk twice and comes back worst of all (0.54).

### What this means for the recommendation

The two beds bracket the regime:

  * **pain, 21 real studies** -- genuinely different paradigms and populations, so prevalence < 1.
    The correction helps: rmse -23%, bias -47%, the top-voxel overestimate eliminated.
  * **HCP synthetic studies** -- one population, prevalence 1. The correction only hurts:
    magnitude 0.63 against 0.85 for the images alone.

A user cannot easily tell which regime a real collection is in, and the estimator currently offers
no diagnostic for it beyond `prevalence` itself -- which is exactly the quantity that is wrong
when it matters. That is the most important open problem on the design, above the report-limb
functional form.

Note this is a stronger statement than the docstring's existing prevalence caveats. Those say
`prevalence` is compressed and should be read ordinally. This says the compression **propagates
into the magnitude** and can make `g` worse than doing nothing, in a regime that is not exotic:
any collection of similar studies of the same effect.

### Correction to the HCP diagnosis: the prevalence is the symptom, not the cause

I wrote that the silences "push mu down through the censoring term *and* pi down through the
mixture, so g_marginal is shrunk twice". Tested it by adding a plain-Tobit option (prevalence
pinned at 1, which is the truth on this bed) -- the machinery was already there, since
`_fit_chunk` sets `pi = 1` whenever the model is not zero-inflated, and only the option list
forbade it.

**It makes the magnitude worse, not better:** 0.60 of the HCP reference against the
zero-inflated 0.63, and 0.411 against 0.422 for a true 0.5 on the field simulator. Removing the
"this study has no effect" escape forces every silence to be explained by a small mu, so mu falls
further.

So the censoring term **over-shrinks mu whatever the prevalence does**, and a fitted pi below 1 is
the model partly *absorbing* that over-shrinkage rather than compounding it. Which makes the HCP
failure the same defect as the overstated `P(report)` (#66), from another direction: too high a
reporting probability against an under-observed count drags mu down, and on a prevalence-1 bed
there is no genuine absence for pi to absorb it into.

Option reverted rather than shipped -- it is never the right choice, and "pin the prevalence"
is the obvious-looking fix that someone else would otherwise reach for. The measurement is in the
docstring so it does not get re-tried.

**This consolidates the open problems.** What looked like three separate defects -- the 18%
deficit at the 0.4 focus, the se/sd of 1.55 to 3.74, and the HCP magnitude at 0.63 -- all point at
one thing: `P(report | mu) = P(|g| >= c | mu)` is too high because a paper reports a local
maximum, not any exceedance. Fixing that is now the single highest-value change to the estimator.

## SDM-PSI vs CBES on HCP, refereed by subjects neither saw: the pain result reverses completely

MOTOR_LH. 480 subjects cut into 16 synthetic studies of 30; coordinates extracted the way papers
produce them (14.4 peaks per table); **306 subjects used to make no coordinate at all** give the
truth. SDM gets all 16 tables (config + `pp` + `mi`, 87 s of MLE); CBES gets 2 of the studies as
g/g_var images and the other 14 as tables; `images only` is those 2 maps pooled.

| estimate | r | rank r | AUC | mag ratio | \|est\| top | \|ref\| top |
|---|---|---|---|---|---|---|
| SDM-PSI coeff (16 tables) | +0.563 | +0.454 | 0.845 | **0.21** | 0.096 | 0.446 |
| images only (2 images) | **+0.845** | **+0.576** | **0.973** | **0.85** | 0.394 | 0.446 |
| CBES g (2 img + 14 tab) | +0.830 | +0.575 | 0.967 | 0.63 | 0.290 | 0.446 |
| CBES g_marginal | +0.785 | +0.564 | 0.943 | 0.54 | 0.259 | 0.446 |

**The confounded pain comparison had SDM at r +0.689 against CBES's +0.358. On a clean reference
it is +0.563 against +0.830 -- reversed.** So the retraction was right, and the clean number is
more favourable to CBES on pattern than I had allowed for. SDM's magnitude is 4.8x low (0.21);
CBES's 1.6x (0.63).

### But the honest reading is not "CBES beats SDM"

The inputs differ by what each method requires. CBES has **two unthresholded maps of 30 subjects
each**; SDM has sixteen thresholded coordinate tables. And `images only` -- those same two maps,
no coordinates at all -- beats both on every metric. So what this bed shows is:

  * two shared maps carry essentially everything,
  * sixteen coordinate tables through SDM are far behind two shared maps,
  * and the coordinate channel adds nothing to the pattern on top of the maps while costing
    magnitude (0.63 against 0.85).

The last point is the HCP prevalence-1 regime already recorded above, not a new fact.

### What would make this a fair method comparison

Give SDM and CBES the same number of *studies* and let each use them as its model requires, with
the image count swept: at 1, 2, 5, 10 images CBES's advantage should shrink toward SDM's as the
maps thin out, and the crossing point is the practically useful number. One SDM run per
configuration is ~2-3 minutes of `mi` plus `pp` on this bed, so a sweep is affordable now that
the CLI works.

### The SDM CLI recipe, since it cost hours to find

`pp` does **not** write the configuration `mi` needs -- the GUI writes `sdmpsi_params.xml`, and
without it `mi` aborts with "Couldn't load configuration file. Please check that you've already
run preprocessing". The file is 50 lines and the only dataset-specific field is `nStudies`
(`nImputs` and `VoxelsMask` may be left at 0). So:

  1. write `sdm_table.txt` (study, n1, t_thr) plus `<study>.spm_mni.txt` of `x,y,z,t`, or
     `<study>.no_peaks.txt` for a study that reported nothing
  2. write `sdmpsi_params.xml`
  3. `sdm_parse pp`   (~15 min on 21 studies at 2 mm; faster here)
  4. `sdm_parse mi`   -> `analysis_MyMean/mi/coeff.nii.gz`, `tau2.nii.gz`, `mle_report.txt`

Also: arguments are **not** space-separated positionals -- `sdm_parse pp gray_matter 1 gray_matter 2`
joins them into a template name and looks for `cor_gray_matter_1_gray_matter_2_2mm.nii.gz`. Bare
`pp` defaults correctly, which is what to use.

### With parity -- SDM given the same two images -- it closes half the gap, and the conclusion sharpens

jdkent's point. ES-SDM takes a mixture by design ("combines reported peak coordinates and
statistical parametric maps"), so the earlier arm was unfair: SDM had 16 thresholded tables while
CBES had 2 unthresholded maps plus 14 tables. Supplying the same two studies to SDM as t maps
(`<study>.nii.gz`; it is strict about the header -- intent must be "t test" with dof, and the
xform codes must be flagged MNI, or it refuses the study silently):

| estimate | r | rank r | AUC | mag ratio | \|est\| top |
|---|---|---|---|---|---|
| SDM-PSI, 16 tables (unfair) | +0.563 | +0.454 | 0.845 | 0.21 | 0.096 |
| **SDM-PSI, 2 img + 14 tab** | **+0.722** | **+0.558** | **0.927** | **0.43** | 0.210 |
| images only, 2 img | **+0.845** | **+0.576** | **0.973** | **0.85** | 0.394 |
| CBES g, 2 img + 14 tab | +0.830 | +0.575 | 0.967 | 0.63 | 0.290 |
| CBES g_marginal | +0.785 | +0.564 | 0.943 | 0.54 | 0.259 |

**Giving SDM the maps closed about half its gap** -- r +0.563 to +0.722, magnitude 0.21 to 0.43 --
so a large part of the earlier difference was input asymmetry rather than method.

Two conclusions, and the second matters more:

1. **On identical inputs CBES beats SDM-PSI on every metric**, pattern and magnitude: r +0.830
   against +0.722, AUC 0.967 against 0.927, magnitude 0.63 against 0.43. (rank r +0.575 against
   +0.558 is inside single-split noise; the r and magnitude gaps are not.)

2. **Both are beaten by pooling the two maps and ignoring all 14 coordinate tables.** Two
   independent methods, given the same data, both do worse with the tables than without them.
   That is much stronger evidence for the prevalence-1 finding than CBES alone could be: it is
   not a quirk of this estimator's censoring term, it is that thresholded coordinate tables carry
   little that unthresholded maps of the same studies do not already carry -- and both methods
   pay for reading them.

Caveat: one split. The r and magnitude differences are large enough to trust; the rank
correlations are not separated.

Contamination checked rather than assumed: CBES's donor g/g_var maps were being written into the
same directory `sdm_parse pp` scans for `<study>.nii.gz`. Moved to a `cbes/` subdirectory, and the
SDM input directory verified to hold exactly `s00.nii.gz`, `s01.nii.gz` and `pp`'s own
`sdm_mask.nii.gz` before the number was believed.

## Retraction: prevalence is NOT what decides whether the correction helps

I shipped "the two beds bracket the regime" -- pain at prevalence < 1 where the correction helps,
HCP at prevalence 1 where it hurts. Built the regime as a dial to check it, and it is wrong.

`does_prevalence_decide_it.py`: `n_studies` studies of 30 subjects, of which `round(pi * n)` are
**MOTOR_LH** studies that carry the effect and the rest are **EMOTION_FACES** studies that do not
-- real subjects, real noise, real reported peaks, just in the wrong places, so they are silent at
the motor voxels for the right reason. Held-out MOTOR subjects give both estimands exactly:
`mu` = their g, `pi*mu` = prevalence x that. 6 splits.

| true pi | reported pi | g / mu | g_marg / pi*mu | images / pi*mu |
|---|---|---|---|---|
| 1.00 | 0.918 | **0.68** | 0.57 | 0.91 |
| 0.75 | 0.714 | **0.63** | 0.60 | 1.08 |
| 0.50 | 0.590 | **0.56** | 0.67 | 1.35 |
| 0.25 | 0.540 | **0.53** | 1.15 | 2.50 |

**`g` under-estimates `mu` at every prevalence and degrades monotonically as prevalence falls.**
It does not improve where the correction supposedly has something to correct. So prevalence is
not the variable.

**Nor is the statistic.** The pain bed reported whole-map rmse and bias; this bed reported a ratio
at the truth's top quartile, and I suspected the absolute-value floor. Measured both on the dial:
rmse `g` 0.130 against images 0.129 at pi = 1, and 0.281 against 0.241 at pi = 0.5. The two arms
tie or the images win on the *pain statistic* too, where on pain CBES won 0.210 against 0.269.

So two candidate explanations tested and rejected, and the difference between the beds is
genuinely unexplained. They differ in at least three ways at once: real studies against synthetic,
real between-study heterogeneity against essentially none, and a 19-study-level reference against
a 300-subject pooled one. The next dial to build is tau2, since that is the one I can add to the
synthetic bed directly.

**One thing this bed did improve on:** the prevalence estimate itself. Against a *designed* pi on
real data it reads 0.918, 0.714, 0.590, 0.540 for true 1.00, 0.75, 0.50, 0.25 -- so it tracks well
down to about 0.5 and then floors near 0.54. That is better than the simulator's compression
(a true 0.25 read 0.49 to 0.60 there) and is now the number in the docstring, since a designed
prevalence on real subjects beats a simulator for this purpose.

## Four candidate explanations tested. Only the *reference construction* moves it the right way.

Why does the correction help on the 21-study pain collection (rmse 0.210 against 0.269) and not on
HCP with held-out subjects? Dials built and measured, all on the MOTOR/EMOTION prevalence bed:

| dial | result | verdict |
|---|---|---|
| prevalence 1.00 / 0.75 / 0.50 / 0.25 | g/mu 0.68, 0.63, 0.56, 0.53 | **rejected** -- degrades monotonically as pi falls |
| statistic (whole-map rmse vs top-quartile ratio) | rmse 0.130 vs img 0.129 at pi=1; 0.281 vs 0.241 at 0.5 | **rejected** -- ties or loses on the pain statistic too |
| between-study tau = 0.0 / 0.3 / 0.6 | rmse g 0.126/0.135/0.159, img 0.126/0.123/0.140 | **rejected** -- no crossing at any tau |
| reference: pooled subjects vs study-level | rmse g 0.126 vs img 0.126, then **0.122 vs 0.125** | **moves it** |

**The mechanism, and it is a caution about the pain number.** The pain reference was an
inverse-variance mean of 19 study-level g maps. Hedges' variance `1/n + g^2/(2n)` is a function of
the *observed* effect, so a study that drew high gets less weight and the pooled reference is
itself pulled **down**. CBES is also biased down. An estimator biased low scores better against a
reference biased low. Reconstructing the HCP reference the same way -- held-out subjects cut into
synthetic reference studies whose g maps are pooled by inverse variance -- takes CBES from tied
(0.126 vs 0.126) to winning (0.122 vs 0.125), and shifts every ratio up (g/mu 0.69 to 0.71,
images 0.92 to 0.95).

**But it accounts for only about 2 points of the pain gap's 22.** So the direction of the
artefact is demonstrated and the magnitude is not: most of the pain result is still unexplained.
The remaining untested difference is real studies against synthetic ones -- different scanners,
paradigms, preprocessing, sample sizes from 14 to 43 -- which cannot be dialled on this bed.

### Two bugs in the bed, both mine, both caught by the numbers

* `shape = EFFECT[used_eff].mean(0)` shadowed the module-level mask `shape` that `report_peaks`
  reads, giving "maximum supported dimension for an ndarray is currently 64, found 29398".
* Adding `delta_k * pooled_mean` to a study's subjects divides a 480-subject mean by a
  30-subject sd -- an unbounded ratio where the latter is tiny, which produced `g/mu` of 4.4e7
  and a "reported prevalence" of 0.094. Fixed by scaling the study's *own* mean, so its g scales
  by exactly `(1 + delta_k)` and tau is between-study spread in the magnitude and nothing else.

Fourth and fifth times this session that a confident-looking number was a defect in the harness.

## The report-probability defect: mechanism confirmed by intervention, scalar fix ruled out

The diagnosis was that the likelihood uses `P(|g| >= c)` where a paper reports a voxel only if it
cleared `c` **and** was a local maximum, overstating the reporting rate by 1.00, 1.78, 1.08, 1.13
at true g of 0.2, 0.4, 0.6, 0.8. Implicated in four deficits: the 0.4-focus shortfall, `se/sd` of
1.55 to 3.74, the HCP magnitude at 0.63, and the prevalence absorbing heterogeneity.

**The fix is not a reweighting** -- `alpha` was tried twice and made things worse. It is that both
limbs must be complementary probabilities of the *same* event. With `E = P(|g| >= c | mu)` and `q`
the chance that a study exceeding here actually names this voxel:

    reported    q E          instead of   E
    silent      1 - q E      instead of   1 - E

Still a proper likelihood. And the arithmetic says where the benefit comes from:

    report limb   score = (qE)'/(qE) = E'/E      -- q cancels entirely
    silent limb   score = -q E' / (1 - q E)      -- scaled by roughly q

So `q` leaves the report limb alone and weakens the *silence* pull by about `q`, which is exactly
the over-shrinkage diagnosed. No smoothness needed to state it.

Measured, 6 collections:

| q | g@0.2 | g@0.4 | g@0.6 | g@0.8 | mean \|err\| |
|---|---|---|---|---|---|
| **1.00 (shipped)** | -23% | **-11%** | +1% | +0% | **8.9%** |
| 0.90 | -22% | -9% | +4% | +4% | 9.4% |
| 0.80 | -20% | -6% | +7% | +8% | 10.3% |
| 0.70 | -19% | -4% | +10% | +12% | 11.0% |
| 0.56 | -16% | **+1%** | +16% | +16% | 12.2% |

**The mechanism is confirmed.** At q = 0.56 -- the reciprocal of the 1.78 measured at the 0.4
focus -- that focus closes exactly, -11% to +1%. Acting on `P(report)` moves mu in the predicted
direction by the predicted amount, which no amount of reasoning had established.

**But no constant q helps overall.** Every value below 1 is worse on mean error, monotonically,
because the required correction is mu-dependent and non-monotone (1.00, 1.78, 1.08, 1.13): a
scalar lifts the weak foci and overshoots the strong ones. So the shipped q = 1 is the best
constant -- the third time this session the shipped configuration has turned out to be the
optimum of the tractable family (the others: alpha = 1 on the report limb, 0 mm report radius).

### What this makes actionable

The fix needs `q(mu)` -- the probability that an exceeding voxel is a local maximum -- which is a
function of the field's **smoothness**. And unlike cluster extent, which jdkent correctly said
papers do not report reliably, **estimated smoothness (FWHM) is routinely reported**: SPM and FSL
both print it, and it appears in methods sections. So the concrete route is a `smoothness`
metadata field feeding an RFT expected-maxima density, which turns this from "needs a quantity we
cannot have" into "needs a field papers already publish". Links to #63.

### Error bars on the 1.78, and a correction to the route I recommended

Put a Poisson CI on the reporting-rate ratio (24 collections, so the counts are the total across
them), because "overstated by 1.78" was a ratio of two small counts with no uncertainty attached:

| truth | reports | q = obs/model | 95% CI | reading |
|---|---|---|---|---|
| 0.2 | 3 | 1.03 | [0.00, 2.20] | **unmeasured** -- my "1.00" was meaningless |
| 0.4 | 21 | **0.58** | **[0.33, 0.84]** | excludes 1: the over-statement is real |
| 0.6 | 163 | 1.00 | [0.85, 1.15] | consistent with no over-statement |
| 0.8 | 299 | 0.89 | [0.79, 0.99] | excludes 1, mildly |

**The 1.78 survives** (1/0.58 = 1.72, CI [1.19, 3.0]), so the defect is real. Two things I stated
with more confidence than the data carry: the value at 0.2 is not measured at all, and **q is not
demonstrably non-monotone** -- the 0.6 and 0.8 intervals overlap heavily, so the shape is "well
below 1 at 0.4, at or just below 1 above the window", which is monotone-increasing-then-flat.

**And that overturns the route I recommended.** I said the fix needs `q(mu)` from reported
smoothness via an RFT expected-maxima density. Checked against the measurement, that has the
*sign backwards*. The RFT clump argument is `q ~ 1/clump size ~ u^3` with `u = (c - mu)/sigma`, so
it predicts q **falling** as mu rises (u goes +1.10, 0.00, -1.10 across the three foci). The data
show q **rising** (0.58, 1.00, 0.89).

The reason is that the RFT density describes a *zero-mean* field, and these are *signal peaks*. At
a strong focus the blob's own curvature makes that voxel the local maximum, so q approaches 1; at
a marginal focus the noise decides which of several exceeding voxels is the max, so q falls. The
governing quantity is the signal's curvature relative to the noise smoothness -- **not** the null
expected-maxima density, and not recoverable from a coordinate table.

So "supply smoothness and compute q from RFT" is withdrawn. What remains established: the defect
is real at the window, acting on `P(report)` closes it (q = 0.56 takes the 0.4 focus from -11% to
+1%), no constant q helps overall, and the correct `q` depends on a signal curvature this model
has no way to observe. Which makes this harder than I claimed an hour ago, not easier.

Also worth recording because it invalidated a hypothesis before I tested it: I suspected the
cluster-forming extent requirement was the mechanism. It cannot be, on this bed -- the field
simulator's rule is `(magnitude == maximum_filter(magnitude, size=3)) & (magnitude >= threshold)`,
every suprathreshold local maximum and no extent test at all. Read the generator, not the
assumption.

## Shipped: `coordinate_share`, the diagnostic all the caveats needed

Every warning accumulated this session is about *when* to trust the coordinate channel, and none
of them is checkable against a collection in hand. What is checkable is whether the channel is
even acting at a voxel -- and the quantity was already being computed and discarded.
`_observed_information` accumulates the images' contribution to `I_mu_mu` and the indicators'
contribution separately before summing, so the ratio is free.

  0  the images carry the estimate alone; the tables changed nothing here, so none of the
     coordinate caveats apply at this voxel
  1  the indicators carry it, and all of them do

On a 20-study collection with two donors: range 0.02 to 1.00, **median 0.10, and 0.79 at the
focus**. So on a typical map the images carry the estimate almost everywhere and the coordinates
take over exactly where studies reported. That is the stratification the whole design rests on --
the 54%-vs-6% bias reduction measured at the start -- now readable per voxel rather than only in
aggregate across a corpus.

Emitted only under `selection_model="zero-inflated"`; with the selection off the tables are inert
and the share is identically zero, so a map would be noise.

### And `coordinate_share` does NOT predict where the interval fails -- tested, falsified

The obvious next inference from the new map: since the `se` excess was located in the censoring
term, `se/sd` should be worse where the share is high. Tested on the four-foci bed, 16
replications, per-voxel spread across them:

| coordinate_share | voxels | se/sd | mean se | sd of g |
|---|---|---|---|---|
| [0.00, 0.05) | 123 | 2.93 | 0.183 | 0.063 |
| [0.05, 0.10) | 537 | 3.21 | 0.192 | 0.060 |
| [0.10, 0.25) | 10530 | 3.23 | 0.187 | 0.058 |
| [0.25, 0.50) | 4344 | 2.97 | 0.223 | 0.075 |
| [0.50, 1.01) | 91 | **2.28** | 0.372 | 0.163 |

Top decile of share 2.86, bottom decile 3.24. **Flat, and slightly better where the coordinates
act most** -- the opposite of the prediction.

So the over-statement is roughly uniform over the map (se/sd ~ 3 everywhere) and does not follow
the channel that produces it. Which also means the earlier reading -- "the excess is located in
the censoring term, since `selection_model="none"` halves the se" -- was conflating things:
switching the selection model off changes the estimator entirely, including the `dof` fallback,
so that comparison was never a within-fit localisation.

The docstring now says explicitly that the share bounds *which caveats apply* and says nothing
about the interval's width, so nobody infers the rule I just failed to establish. And a uniform
3x over-statement points at something structural -- the conditional-vs-marginal mismatch, or an
incomplete Louis missing-information correction -- rather than anything per-voxel.

## Retraction: every `se/sd` on record was computed on `sd(|g|)`, and so was every "bias" in the quiet stratum

Suspect the test before the theory, again -- and this one had been shipped.

`g` is a signed inverse-variance mean. The beds all did

    gs.append(np.abs(res.get_map("g", return_type="array").ravel()))

before taking a spread across replications. At a voxel whose truth is zero, `sd(|g|)` is about
0.6 of `sd(g)`, and `mean(|g|)` is about 0.8 `sd(g)` above zero rather than at it. So the
absolute value did two things at once: it shrank the denominator of `se/sd` and it manufactured a
positive bias out of nothing. Both landed in the class docstring.

**The calibration arm that should have caught it.** Give all 20 studies images and switch the
selection model off, and the fit is a textbook local inverse-variance random-effects
meta-analysis, where `se/sd` must be near 1:

| arm | stratum | bias | sd(g) | sd(&#124;g&#124;) | se | se/sd | se/sd(&#124;g&#124;) | cov |
|---|---|---|---|---|---|---|---|---|
| A textbook IVW (20 images, no selection) | null | -0.001 | 0.0378 | 0.0226 | 0.0433 | **1.15** | 1.92 | 0.98 |
| A | g=0.2 | -0.008 | 0.0372 | 0.0372 | 0.0435 | 1.17 | 1.17 | 1.00 |
| A | g=0.6 | -0.025 | 0.0426 | 0.0426 | 0.0460 | 1.08 | 1.08 | 1.00 |
| B images only (2 images, no selection) | null | -0.001 | 0.1225 | 0.0738 | 0.1503 | 1.23 | 2.04 | 1.00 |
| B | g=0.6 | -0.025 | 0.1317 | 0.1317 | 0.1555 | 1.18 | 1.18 | 1.00 |
| C shipped (2 images, zero-inflated) | null | -0.001 | 0.1120 | 0.0644 | 0.2462 | **2.20** | 3.82 | 0.98 |
| C | g=0.2 | -0.084 | 0.1041 | 0.0868 | 0.2036 | 1.96 | 2.35 | 1.00 |
| C | g=0.6 | -0.016 | 0.0836 | 0.0836 | 0.1156 | 1.38 | 1.38 | 1.00 |

Arm A reads 1.08 to 1.17, so the harness is sound. Note the `sd(|g|)` column: at the null
stratum it is 1.7x the signed one in every arm, and at the foci it is identical, because there
the estimate is far enough from zero that the absolute value is a no-op. That is exactly the
signature -- and it explains why the docstring's *quiet* column ran 2.71 to 3.74 while its
*effect* column ran 1.55 to 1.97.

**What the corrected numbers are.** Re-running the shipped interval table with signed `g`, 40
replications: quiet bias goes from +0.114 to **-0.000** at one image and from +0.093 to
**+0.001** at two; quiet `se/sd` from 2.71 to **1.51** and from 3.17 to **1.80**; the effect
stratum moves little, 1.55 to 1.39 and 1.78 to 1.66. So:

  * **The estimator is not biased upward where there is no effect.** That claim was an artifact
    of the absolute value and has to come out of the docstring.
  * **`se/sd` is 1.3 to 1.9, not 1.55 to 3.74.** Still conservative, still worth fixing, but not
    the three-fold failure on record.
  * **Some of the conservatism is not the censoring term at all.** A textbook all-image IVW fit
    reads 1.08 to 1.17 in this bed, so the baseline is already above 1 (DL `tau2` on the se, and
    the conditional-vs-marginal distinction). The shipped fit's remaining gap is 1.6 to 1.9.

**And the share result reverses.** Re-run at 48 replications on signed `g`:

| coordinate_share | voxels | se/sd | mean se | sd of g | bias |
|---|---|---|---|---|---|
| [0.00, 0.05) | 25 | 1.87 | 0.1935 | 0.1036 | +0.0029 |
| [0.05, 0.10) | 445 | 1.82 | 0.2036 | 0.1119 | -0.0031 |
| [0.10, 0.25) | 11037 | 1.87 | 0.1923 | 0.1028 | -0.0023 |
| [0.25, 0.50) | 4040 | 1.63 | 0.2219 | 0.1358 | -0.0010 |
| [0.50, 1.01) | 78 | **1.33** | 0.3353 | 0.2526 | -0.0017 |

Top decile 1.54, bottom decile 1.86. Not flat: **monotone decreasing in the share**, so the
interval is better calibrated where the coordinates act most. And not a signal effect -- only a
few hundred of these 15625 voxels carry any truth, so the 4040 voxels in the [0.25, 0.50) band
are nearly all null, and the gradient holds among them. The bias column is now flat at zero
across every band, which is the same retraction again.

**But do not read that as "the coordinates fix the interval."** Look at `sd of g`: it *rises*
with the share, 0.103 to 0.253. The coordinate channel adds real sampling variability -- which
voxels a study reports is itself random across replications -- and `se` rises less than `sd`
does. So two errors are cancelling: a baseline `se` that is too conservative, and a censoring
term whose contribution to the sampling variance is under-counted. `se/sd` near 1 at high share
is an accident of their meeting, not a calibrated interval.

That is a better-posed open question than the one I had: not "why is the se 3x too big" (it is
not) but **"why does the observed information under-count the variance the indicator injects,
while the baseline over-counts?"**

**Scope of the contamination.** Bounded, and worth being precise about. The defect only bites
where a spread of `|g|` is called the estimator's sampling spread, or a bias of `|g|` is read
against a near-zero truth. In `hcp_sdm_vs_cbes.py`, `does_prevalence_decide_it.py` and
`validate_redesign.py` the *truth is taken in absolute value too* (`truth = np.abs(hedges(...))`),
so those are magnitude-against-magnitude comparisons with the same transformation on every arm:
the CBES/SDM-PSI/images-only rankings, the HCP held-out result and the pain split-half rmse all
stand. What falls is the interval section and the share paragraph -- both docstring content.

Habit to add to the list that already has "report n and SE beside any number" and "report
relative error per stratum": **a calibration arm is not optional, and it has to be an arm whose
answer is known analytically.** Arm A costs one line in the arm list and would have caught this
the first time `se/sd` was ever printed.

## The inflation was an unidentified prevalence, and it is fixed

The retraction above left a better-posed question: the baseline was near 1 and the excess did not
follow the share the way it should if the censoring term produced it. The configuration that
settles it turns out to be the one nobody was looking at.

**Give every study an image.** Then the reporting indicator is *structurally* empty -- image
studies are omitted from it at every voxel -- so the coordinate channel cannot be responsible for
anything. Direct check:

    2 images: coordinate_share max 1.000000, nonzero voxels 15625; mean prevalence 0.5605
   20 images: coordinate_share max 0.000000, nonzero voxels     0; mean prevalence 0.5773

Yet that arm measured *worst* of any, in one bed on shared seeds:

| arm (20 images, 20 studies) | null | g=0.2 | g=0.4 | g=0.6 | g=0.8 |
|---|---|---|---|---|---|
| A selection off (textbook IVW) | 1.14 | 1.24 | 1.02 | 1.00 | 1.05 |
| D zero-inflated, indicator empty | **2.21** | **2.30** | 1.28 | 0.99 | 1.06 |

What it was paying for is visible in the numbers above: **a prevalence of 0.577 against a true
1.000, fitted from nothing.** Only the indicator separates "no effect in this study" from "a
small effect plus noise"; a two-component mixture fitted to 20 Gaussian values with known
variances will explain noise as a mixture, and does. The uncertainty in that phantom parameter is
then profiled out of the information about `mu` by the Schur complement. And the signature fits:
the inflation is concentrated at null and weak voxels (2.21, 2.30) and absent at strong ones
(0.99, 1.06), because where the effect is weak "a small effect in every study" and "a large
effect in a few" fit equally well.

**The fix.** Hold `pi` at 1 wherever no study contributes an indicator, and skip the profiling
there. Two subtleties, both load-bearing:

  * The limit cannot be left to the algebra. With `pi` clamped at `1 - 1e-6`, `(1 - r)/(1 - pi)`
    tends to the *ratio of the two component densities*, not to zero, so the cross block stays
    order one and the Schur complement would still subtract a term no parameter earned. My first
    patch assumed it vanished; the test caught it at 7e-5 relative.
  * Forcing the responsibility to exactly 1 (not just near it) makes the collapse an identity.
    The all-image mixture fit now returns the non-mixture fit bit for bit, which is asserted as
    `np.allclose` on `se` and `g` rather than as a loose tolerance.

Post-fix, arm D is bit-identical to arm A in every stratum, and the shipped interval table's
all-image row goes from se/sd 2.23 / 2.26 / 1.55 to **1.16 / 1.22 / 1.19**, with the half-width
at the effect from 0.34 to **0.26**. Arms with any indicator at all are unchanged to the last
digit, which is the check that matters: this is a degenerate-case guard, not a new estimator.

### Two other causes ruled out along the way

**The reporting rule is not it.** The model censors on `|g| < c`; a paper reports local maxima
above a threshold. Those differ, and a voxel on a blob's shoulder can sit far above the cutoff
and still be silent. Feeding the estimator its own rule instead -- report every supra-threshold
voxel, 386 foci per collection against 172 -- changes nothing:

| arm | null | g=0.2 | g=0.4 | g=0.6 | g=0.8 |
|---|---|---|---|---|---|
| model-matched (&#124;g&#124; >= c) | 1.80 | 1.61 | 1.18 | 1.80 | 1.90 |
| as published (maxima) | 1.83 | 1.61 | 1.18 | 1.79 | 1.90 |

Bias is likewise identical (-0.048 against -0.048 at g=0.2). Doubling the reported foci and
making the silence mean exactly what the model thinks it means moves neither the estimate nor its
error. That is a strong negative: the over-shrinkage at the window is not the local-maximum
condition.

**The bed contributes 4%.** The simulator's per-study `g` has an empirical sd across studies of
0.1797 at null voxels against a reported `sqrt(g_var)` of 0.1873, so the noise sits 4% below the
`1/n + g^2/(2n)` the model assumes and every `se` is 4% generous by construction. Small, but it
is part of why a textbook arm reads 1.14 rather than 1.00, and it belongs on the record.

### What is left, and the next prediction

A conservative interval where the coordinates *do* act: 1.28 to 2.03 at two images, falling
monotonically as `coordinate_share` rises (1.87 -> 1.33). That is now the same mechanism seen
from the other side -- the share is how much the indicator says about `pi`, so more share means a
better-determined prevalence and a smaller profiling penalty.

Which makes a sharp prediction. The `pi`/`mu` ridge runs along curves of roughly constant
`pi*mu`, so the *product* should be far better determined than either factor. `g_marginal` and
`se_marginal` are already emitted, so if the ridge is the mechanism, `se_marginal/sd(g_marginal)`
should sit near 1 exactly where `se/sd` for `g` is worst. Running.

### Prediction falsified: `g_marginal` is the steadier estimate and the worse interval

The ridge argument said the product should be the determined combination, so `se_marginal` should
be the sound width. Measured, 32 replications, true `pi` = 1:

| 2 images, 18 tables | sd(g) | se | se/sd | sd(g_m) | se_m | se/sd_m | fitted pi |
|---|---|---|---|---|---|---|---|
| null | 0.1116 | 0.2041 | 1.83 | 0.0762 | 0.3491 | **4.58** | 0.551 |
| g=0.2 | 0.1280 | 0.1877 | 1.47 | 0.1084 | 0.2215 | 2.04 | 0.627 |
| g=0.4 | 0.0803 | 0.1323 | 1.65 | 0.0930 | 0.2305 | 2.48 | 0.858 |
| g=0.6 | 0.0785 | 0.1543 | 1.97 | 0.0877 | 0.2134 | 2.43 | 0.903 |
| g=0.8 | 0.0755 | 0.1300 | 1.72 | 0.0830 | 0.1311 | 1.58 | 0.960 |

**Half the prediction held and it was the useless half.** The product *is* the steadier estimate
-- 0.0762 against 0.1116 at quiet voxels -- so the ridge really does run roughly along constant
`pi*mu`. But its reported error is worse, 0.349 against 0.204, and `se_marginal/sd` beats `se/sd`
in only the strongest stratum.

The reason is mechanical and I should have seen it before running: `Var(mu*pi)` needs the *whole*
2x2 inverse, not a Schur complement, so a near-singular information matrix amplifies there rather
than cancelling. The ridge makes the determinant small; the delta method divides by it.

So the conclusion flips the remedy. The width does not need a different estimand -- it needs a
variance that does not invert a near-singular matrix. A **profile likelihood** does exactly that:
no inverse, and it sidesteps the `dof` question (#62) in the same stroke, since a profile interval
needs no degrees of freedom at all. That is now the one concrete route left for the interval, and
the first thing to try rather than the last.

(The 20-image arm also re-confirms the fix from a second angle: with `pi` held at 1, `se_marginal`
equals `se` exactly in every stratum, as it must when `mu*pi = mu`.)

### The scalar bed cannot speak to the width, and said so on its own calibration row

Before building a profile-likelihood interval into the estimator, the scalar version: the same
censored zero-inflated mixture at one voxel, known truth, both intervals from the same fit, 300
replications. Wald is `1.96 / sqrt(Schur complement)` by finite differences -- a second
implementation of the shipped quantity rather than a copy of the algebra. Profile is
`{mu : 2 (l_p(mu_hat) - l_p(mu)) <= 3.841}` with `pi` re-maximised at each `mu`.

| configuration | cov W | cov P | half W | half P | se/sd W | se/sd P |
|---|---|---|---|---|---|---|
| no indicator, pi fixed at 1 | 0.94 | 0.94 | 0.268 | 0.268 | **0.99** | **0.99** |
| no indicator, pi free | 0.92 | 0.96 | 0.269 | 0.268 | 0.78 | 0.78 |
| 2 images, 18 indicators | 0.99 | 0.99 | 0.200 | 0.162 | 0.74 | 0.60 |
| 2 images, 18 indicators, pi 0.5 | 0.92 | 0.98 | 0.257 | 0.377 | 0.19 | 0.28 |

**The calibration row passes and the next one fails, which is the finding.** Reducing to a
textbook Gaussian mean gives 0.99 for both intervals, so the bed's arithmetic is right. But the
`2 images, 18 indicators` row reads `se/sd` of **0.74**, where the full-scale bed reads 1.83 at
quiet voxels and 1.28 to 1.60 at the foci for the nominally identical configuration. Below 1
instead of above: the 1-D bed does not reproduce the phenomenon at all, so **nothing it says
about the remedy is evidence.** Per the standing rule, that is where to stop rather than where to
start believing it.

Worth naming what differs, because whichever it is, is the actual mechanism at full scale and the
scalar bed has just excluded a long list of candidates:

  * no `tau2` at all in the scalar bed, against a per-voxel DerSimonian-Laird estimate from
    *two* studies at full scale -- noisy and upward-biased there, and it enters the `se` directly.
    Against this: the `tau = 0.3` arm moved `se/sd` from 1.80 to 1.80 at quiet voxels.
  * `1.96` against a `t` on the censoring roster's `dof`.
  * indicator signs drawn i.i.d. per study, against signs generated by a *smooth field* whose
    peaks are local maxima -- spatially correlated across studies and far more variable across
    replications than a binomial.
  * the replication spread taken at a fixed voxel while the whole field and all its peaks are
    redrawn, which is a wider marginal than the scalar bed's.

**What the scalar bed does establish, being a statement about the likelihood rather than about
the data:** the profile interval tracks the likelihood's actual shape where the Wald one does
not. It is narrower when the Wald is too wide (0.162 against 0.200) and *wider* when the Wald is
too narrow (0.377 against 0.257), and its coverage is closer to nominal in both directions --
0.96, 0.99, 0.98 against 0.92, 0.99, 0.92. That is the correct behaviour and a reason to prefer
it. The size of the gain at full scale is simply not measurable here.

Note also the second row: with no indicator and `pi` free, `se/sd` falls from 0.99 to 0.78 and
coverage from 0.94 to 0.92 -- the unidentified-prevalence cost showing up in the scalar bed too,
in the same direction, which is a small independent confirmation of the fix now shipped.

**Next, and it has to be at full scale:** implement the profile interval in the estimator and
measure it on the field bed against the arm table above. The scalar bed has done the one thing it
could do honestly, which was to rule itself out.

### The DerSimonian-Laird `tau2` is part of the conservatism and still cannot be removed

The scalar bed's biggest structural omission was `tau2`, so it was worth testing at full scale --
and unlike the `tau = 0.3` arm, which changes the *truth* and so moves spread and error together,
this changes what the estimator does. Four foci, 32 replications, `tau2_method="none"` against the
default:

| arm | null | g=0.2 | g=0.4 | g=0.6 | g=0.8 |
|---|---|---|---|---|---|
| B images only, tau2 dl | 1.22 | 1.30 | 1.05 | 1.22 | 1.23 |
| F images only, tau2 off | **1.07** | 1.09 | 0.92 | 1.12 | 1.09 |
| C shipped, tau2 dl | 2.03 | 1.60 | 1.28 | 1.39 | 1.29 |
| E shipped, tau2 off | 1.81 | **1.02** | 1.14 | 1.28 | 1.16 |

At two images and `tau = 0.1` that looks like most of the answer: the images-only baseline becomes
calibrated (1.07, and 0.92 to 1.12 across the foci), and the shipped fit becomes nearly so
everywhere except the quiet stratum. A DerSimonian-Laird estimate from *two* studies is noisy and
upward-biased, and it enters the `se` directly while barely widening the estimator's own spread.

**And it is not a lever, which is why the follow-up arms mattered.** Raise the true heterogeneity
to `tau = 0.3` and keep everything else, and switching `tau2` off stops being conservative and
starts being wrong:

| arm (2 images, true tau 0.3) | g=0.2 se/sd | g=0.2 cov | g=0.4 cov | g=0.8 cov |
|---|---|---|---|---|
| I tau2 dl | 0.99 | 0.86 | 1.00 | 0.94 |
| J tau2 off | 0.52 | **0.45** | 0.88 | 0.94 |

Coverage at the weakest focus falls to 0.45 -- an interval that misses the truth more often than
it catches it. At five images the picture is mixed rather than favourable (quiet 2.51 with `tau2`
off against 2.28 with it, `g=0.4` 1.40 against 1.69). So the correct reading is not "prefer none"
but:

**The estimated `tau2` contributes to the conservatism at small k and low heterogeneity, and it
is buying real protection at the same time.** Removing it trades a 1.3x-too-wide interval for a
0.45-coverage one. The lever is tested and closed; it belongs on the list of things that look
like the fix and are not, next to the RFT route and the constant `q`.

One thing to follow up separately: arm I under-covers at the weakest focus even *with* `tau2`
(0.86 at `g=0.2`, se/sd 0.99, bias -0.047 on a spread of 0.37). Conservative at low heterogeneity
and anti-conservative at high, at the same voxel -- which is a different failure from the one this
section is about.

## The algebra was available the whole time, and it corrects two results I simulated

jdkent asks whether some of this can be done algebraically, pointing at a proofs repo. The answer
is yes, and this session is the argument for it: two of its results were reached by staring at
numeric tables, and a third was simply wrong in a way one line of algebra would have prevented.

### 1. The profile interval's asymptote -- derivable, and sharper than the simulation

Write each observation's mixture density as `f_i(mu, pi) = pi a_i(mu) + (1 - pi) b_i`, with `b_i`
the density or probability under an effect of exactly zero. As `|mu| -> infinity`:

  * an image value has `a_i(mu) = phi((g_i - mu)/sigma)/sigma -> 0`;
  * a *silence* has `a_i(mu) = P(|g| < c | mu) -> 0`;
  * a *report* has `a_i(mu) = P(|g| >= c | mu) -> 1`.

So the log-likelihood tends to a limit free of `mu`,

    l_inf = max_pi [ sum_{i not reported} log((1 - pi) b_i)
                     + sum_{i reported}   log(pi + (1 - pi) b_i) ],

a horizontal asymptote. The profile interval on `mu` is bounded exactly when
`2 (l_hat - l_inf)` clears the critical value, and at a voxel where nobody reported -- where the
max is at `pi -> 0` -- that is precisely the likelihood-ratio test of `pi = 0`. This is *sharper*
than what I concluded from the grid: I wrote "bounded iff the data reject pi = 0", which is only
exactly true where no study reported at the voxel. The report limb keeps `pi` in the limit. The
docstring now carries the general form.

Note what the derivation also gives for free: **no reparametrisation escapes it.** `pi*mu` has
the same asymptote, approached along `pi -> 0, mu -> infinity`. I could have known that before
building anything.

### 2. The ridge invariant -- I guessed `pi*mu` and the algebra says otherwise

A silent pair's likelihood depends on `(pi, mu)` only through the probability of the event
observed,

    P_silent(pi, mu) = pi S(mu) + (1 - pi) S(0),     S(mu) = P(|g| < c | mu),

one equation in two unknowns. So the silence channel's indifference curve is the level set of
`P_silent` -- not of `pi*mu`. Checked against the actual profile ridge at a quiet voxel (15
silent pairs, 0 reported):

| mu | pi | pi*mu | P_silent |
|---|---|---|---|
| -0.1420 | 0.5754 | -0.0817 | 0.99814 |
| +0.0182 | 0.4730 | +0.0086 | 0.99959 |
| +0.1784 | 0.2362 | +0.0421 | 0.99843 |
| +0.4988 | 0.0208 | +0.0104 | 0.99529 |
| +1.2998 | 0.0001 | +0.0001 | 0.99950 |

`P_silent` is constant to 0.4%; `pi*mu` moves by a factor of several hundred and changes sign.
**So the #70 entry above is right about the measurement and wrong about the mechanism.** The
product is *not* the invariant, and the reason `g_marginal` is the steadier estimate is not that
it lies along the ridge. That entry should be read as: measured, the product is steadier and its
delta-method interval is worse; the ridge explanation offered for it does not hold.

### 3. The one the algebra would have caught before the test did

Holding `pi` at 1, I claimed the cross block "vanishes through the existing algebra" because the
responsibility goes to 1 and `r(1 - r) -> 0`. With `pi` clamped at `1 - eps`:

    1 - r = eps b_i / ((1 - eps) a_i + eps b_i) ~ eps b_i / a_i,
    so (1 - r)/(1 - pi) -> b_i / a_i,    and  r(1-r)/(pi(1-pi)) -> b_i / a_i,

both order one, not zero. The Schur complement would have kept subtracting a term no parameter
earned. A unit test caught it at 7e-5 relative; two lines of limits would have caught it first.

### What is worth proving, in rough order of what it would buy

1. **Identifiability.** The above, stated properly: from silences alone the likelihood depends on
   `(pi, mu)` through one scalar per distinct `(sigma, c)` pair, so `k` distinct study
   configurations give `k` equations. This predicts *when* the ridge collapses -- heterogeneous
   sample sizes and thresholds should identify both parameters -- and that is a testable
   prediction I have never made, only stumbled toward with `coordinate_share`.
2. **The asymptote and the boundedness criterion.** Done above; belongs written up.
3. **The downward pull of Hedges' variance.** `v(g) = 1/n + g^2/(2n)` is increasing in `|g|`, so
   inverse-variance weights are anticorrelated with the observed magnitude and the pooled mean
   is biased toward zero. This explains the -0.02 to -0.09 bias at the foci that shows up in
   *every* arm table including the textbook one, which I have been reporting for weeks without
   deriving.
4. **`Var(pi mu)` under a near-singular information matrix.** The delta method divides by
   `det = I_mm I_pp - I_mp^2`; the ridge drives `det -> 0`. That is #70's real mechanism and it
   is two lines.
5. **Concavity in `pi` at fixed `mu`**, which the profile's inner solve relies on: `log` of an
   affine function of `pi`, summed. One line, and it licenses the fixed-point iteration I am
   already using as if it were a global maximiser.

Items 1 and 3 are the ones that would change what gets built, rather than explain what already
was.

## The algebra pays: heterogeneous sample sizes identify the magnitude, heterogeneous thresholds do not

This is the first result here that came from the algebra *first* and the simulation second, and it
is also the first that bears on whether to run CBES on a given collection at all.

**Derivation.** A silence constrains `(pi, mu)` only through `P = pi S(mu) + (1 - pi) S(0)` with
`S(mu) = P(|g| < c | mu)`. Studies with different `(sigma, c)` give different equations -- but the
reporting cutoff measured in sampling standard deviations is just the reported statistic back
again. Checked numerically: `c / sigma` is 2.31 to 2.42 for `z = 2.3` across `n` from 12 to 120,
the drift being the t-versus-normal correction. So

    S(0) = 2 Phi(z) - 1,                                  depends on the THRESHOLD alone
    S(mu) = Phi(z - mu sqrt(n)) - Phi(-z - mu sqrt(n)),    depends on the threshold AND mu sqrt(n)

Varying the threshold moves both components together through the same tail, leaving the equations
nearly collinear. Varying the sample size moves the active component through `mu sqrt(n)` while
leaving `S(0)` **exactly** unchanged. That contrast is what separates the two parameters.

**Test.** 20 studies in every arm, so nothing is about having more data. The bounded fraction of
the profile interval is the probe, since the interval is finite exactly where the likelihood pins
`mu` down -- an identifiability measure rather than a proxy for one. True `pi` is 1.0:

| studies | bounded overall | bounded near signal | fitted pi | mean abs g | truth |
|---|---|---|---|---|---|
| alike, n 28-32, one cut | 0.025 | 0.422 | 0.823 | 0.225 | 0.238 |
| **sample size spread 12-120** | 0.035 | **0.691** | **0.913** | 0.238 | 0.238 |
| threshold spread 2.3-4.5 | 0.026 | 0.414 | 0.815 | 0.235 | 0.238 |
| both spread | 0.035 | 0.676 | 0.901 | 0.238 | 0.238 |

Spreading the sample size moves the bounded fraction near signal from 0.42 to 0.69 and the fitted
prevalence from 0.82 to 0.91. Spreading the threshold moves nothing at all -- 0.414 against 0.422,
0.815 against 0.823 -- exactly as the cancellation says. And "both" is indistinguishable from
"sample size only", which is the sharp version of the claim: the threshold contributes nothing
even in combination.

**So: a magnitude is recoverable from a literature of widely differing sample sizes and not from a
literature of uniformly sized studies, however many of them there are.** That is the opposite of
the usual posture toward heterogeneity, it is now in the class docstring, and it is the piece of
guidance a user actually needs before pointing this estimator at a collection.

It also retrospectively explains the pain-versus-HCP disagreement that four dials failed to
explain (prevalence, statistic, tau2, reference construction). The HCP bed was built from
*equal-sized* synthetic studies; the NIDM pain collection has real, widely varying sample sizes.
That is a testable prediction rather than a story, and it is the next thing to check: re-run the
HCP comparison with sample sizes spread over the same range as pain's.

### Prediction falsified: HCP is not in the unidentified regime, so spread cannot rescue it

The prediction at the end of the last section was that HCP's verdict is a property of its
equal-sized studies. Tested with the subject budget, the study count, the held-out truth and the
*image studies' sizes* all held fixed, so only the tables' heterogeneity changes and the
images-only baseline is identical by construction. MOTOR_LH, 2 image studies of 30 plus 14 tables
sharing 420 subjects, 3 splits:

| arm | table sizes | bounded (top decile) | fitted pi | CBES ratio | images-only ratio |
|---|---|---|---|---|---|
| uniform | 14 x 30 | 0.972 | 0.675 | 0.680 | 0.921 |
| spread | 9 to 70 | 0.970 | 0.660 | 0.671 | 0.921 |

`r` and AUC are identical to three decimals in both arms (0.844/0.848, 0.968/0.969). **Spread
does nothing here, and the reason is visible in the `bounded` column: it is 0.97 in both arms.**
HCP's motor signal is strong enough that `pi = 0` is decisively rejected at the voxels being
scored, so `mu` is already identified and there was no deficit for heterogeneity to repair.

Two things follow, and both matter more than the failed prediction.

**The identifiability result is unharmed but its scope is now known.** It was derived and
confirmed on a bed whose foci run 0.2 to 0.8, where the bounded fraction near signal is 0.42 and
heterogeneity lifts it to 0.69. On a collection with a large effect and ~14 peaks per table the
fraction is already 0.97 and the derivation predicts, correctly, that spread then buys nothing.
Identifiability is necessary, not sufficient.

**HCP's shortfall is a genuine bias, not a ridge artifact.** That is the useful half. With `mu`
identified at 97% of the scored voxels, CBES still recovers 0.68 of the truth where images-only
recovers 0.92, and still fits `pi = 0.66` against a true 1.000. A flat likelihood cannot be
blamed for either. So the over-shrinkage at the window and the under-estimated prevalence are
properties of the censoring term's *specification* -- what it assumes a silence means -- rather
than of what the data can support.

**That is now five dials tested and rejected** for the pain-versus-HCP disagreement: prevalence,
the reported statistic, `tau2`, reference construction (2 of 22 points), and sample-size spread.
The remaining structural difference between the two collections that has never been isolated is
the one in #37: pain's tables were *transcribed from papers* while HCP's are extracted by
`reporting.py` from maps. Different peak-selection processes, not different statistics -- and it
is the last candidate standing.

## Correction: the pain tables were never the published ones either

I wrote, twice, that the last structural difference between the pain and HCP beds is that
"pain's tables were transcribed from papers while HCP's are extracted by `reporting.py` from
maps". That is false, and reading `validate_redesign.py` shows it plainly:

    for i in coord_members:
        foci, height = report_peaks(maps[i], mask_bool, shape, zooms, scheme=SCHEME, focus=FOCUS)

**Both beds extract from images with the same function.** The pain collection's own 267
transcribed coordinate rows -- present, across all 21 studies -- were loaded and then ignored.
So #37's premise was wrong and the "last candidate standing" was not a candidate at all.

What it opens instead is better, and is the assumption under *every* real-data result in this
program. The standing instruction is to extract coordinates "the way papers produce them", and
`reporting.py` is that instruction's implementation -- multiplicity correction, whole clusters,
8 mm separation, no cap. It has never once been checked against coordinates a paper actually
printed. Pain supplies both for the same 21 studies, so the check is available:

  * **as published** -- pain's transcribed tables, with `clamp_threshold` supplying the height
    bound from the smallest reported statistic, since a paper's own cut is not recorded;
  * **re-extracted** -- `report_peaks` on the same studies' z maps, which is what every result
    here used.

Descriptively the two differ before any estimator runs: the transcribed tables have a median of
12 peaks per study (range 4 to 24) over 21 studies. The extracted ones are whatever survives
each scheme, which the runs above put nearer 14 per table on HCP. Worth measuring directly:
peaks per study, and the distance from each published peak to the nearest extracted one.

If extraction is faithful, every real-data result stands as measured. If it is not, the
`reporting.py` proxy is doing work nobody has audited, and the pain-versus-HCP disagreement has
a candidate again -- this time a real one.

Also note what pain's sample sizes are: 9, 9, 9, 12, 12, 12, 12, 12, 12, 13, 14, 14, 16, 16, 16,
20, 20, 24, 25, 25, 32. A 3.5-fold spread, all small. Under the identifiability result that is a
*partially* identified regime, which is consistent with HCP (uniform 30) sitting at 0.97 bounded
while the synthetic 0.2-to-0.8 bed sat at 0.42 -- pain has not been measured on that axis at all.

## The extraction is NOT a faithful proxy for published tables, and `cluster` is the worst offender

First audit of the assumption under every real-data result here. NIDM pain carries 267
transcribed peaks across 21 studies *and* their z maps, so `reporting.py` can be compared
against what the papers actually printed.

| scheme | studies extracting | peaks total (published 267) | median per study (published 12) | median published-to-nearest-extracted | within 8 mm | within 20 mm |
|---|---|---|---|---|---|---|
| **cluster** | 21/21 | **119** | 5 | **23.2 mm** | **0.23** | **0.42** |
| fdr | 16/21 | 1764 | 110 | 7.5 mm | 0.47 | 0.84 |
| fwe | 16/21 | 273 | 6 | 48.1 mm | 0.21 | 0.33 |

**The positive control is the `fdr` row, and it is what licenses reading the rest.** If study
ids had been mismatched -- each paper's peaks compared against another study's map -- no scheme
could recover 84% of published peaks within 20 mm. It does, so the coordinates genuinely
correspond to these maps and the failures below are properties of the scheme rather than of the
harness. (The id matching was checked directly as well: the images frame has an integer index, so
the study id comes from the `study_id` column.)

Three things, in order of how much they hurt.

**`cluster` -- the scheme behind essentially every real-data number in this program -- recovers
23% of published peaks within 8 mm and 42% within 20 mm.** So a *majority of published peaks lie
outside even the coverage radius a silence asserts over*. Those are fabricated silences: the
model is told "no study reported near here" at voxels where a paper did report. It also finds 119
peaks where the papers printed 267, and 5 per study where they printed 12.

**`fwe` matches the count almost exactly and the locations not at all.** 273 extracted against
267 published, which looks like a bullseye, while the median published peak sits 48 mm from
anything it found and only 21% land within 8 mm. Count agreement is no evidence of fidelity, and
I would have accepted it as such.

**`fdr` is spatially closest and quantitatively absurd**, at 110 peaks per study against a
published 12 -- and it silences 5 of 21 studies entirely, which `cluster` does not.

So no scheme reproduces published tables. What this does and does not invalidate:

  * It does **not** invalidate comparisons *between arms* -- CBES against images-only against
    SDM-PSI all saw the same extracted tables, so the rankings stand as measured.
  * It does mean **none of those results describes what CBES would do on a transcribed
    literature**, which is what a user would actually hand it. The headline pain result -- rmse
    0.210 against images' 0.269 -- was measured on tables bearing this much resemblance to the
    papers it came from.
  * It does not explain the pain-versus-HCP disagreement, since both used `cluster`.

**The next run is therefore the important one**, and it is now cheap: re-run the pain split-half
with pain's *published* tables, `clamp_threshold` supplying the height bound from the smallest
reported statistic, against the same reference. If the -22% rmse survives on real tables, the
claim is about the estimator. If it does not, it was about `report_peaks`.

## The headline result survives on the published tables, and is better centred there

The audit above made this the most consequential outstanding run: CBES's central real-data claim
had only ever been measured on `report_peaks` output that recovers 23% of pain's published peaks
within 8 mm. Pain carries its own 267 transcribed peaks, so the claim can be tested on them.

Published tables get no recorded height -- a paper does not print its cut -- so the assumed cut
is the conventional p < 0.001 (z = 3.09) and `clamp_threshold` lowers it per study to that
study's smallest reported statistic, which is a hard upper bound on whatever it really used.
Each transcribed point is given its study's own z-map value at that location, the closest honest
stand-in for the statistic the table would have carried. Eight splits, same reference, same arms:

| estimate | bias | at top | rmse | rank r | AUC |
|---|---|---|---|---|---|
| images only, pooled | +0.136 | +0.137 | 0.269 | 0.484 | 0.893 |
| CBES silence off | +0.142 | +0.155 | 0.272 | 0.499 | 0.899 |
| **CBES `g`, published** | **+0.069** | **+0.005** | **0.230** | 0.480 | 0.886 |
| CBES `g`, extracted | +0.075 | -0.065 | 0.210 | 0.486 | 0.889 |
| CBES `g_marginal`, published | -0.064 | -0.291 | 0.184 | 0.483 | 0.887 |

Paired against images-only over the 8 splits, published arm: rmse -0.040 (p = 0.0024), mean
error -0.067 (p = 0.0001), error at top -0.132 (p = 0.0001), rank r -0.004 (p = 0.84), AUC
-0.006 (p = 0.30), Pearson r -0.030 (**p = 0.075**, no longer significant).

**So the claim is about the estimator and not about `report_peaks`.** Three things worth keeping:

  * The bias at the strongest voxels goes to **+0.005** against images-only's +0.137. That is
    the best-centred top stratum of any arm in the whole program, and it is on real coordinates.
  * The Pearson cost that was significant on extracted tables (0.042, p = 0.017) is smaller and
    not significant on published ones (0.030, p = 0.075). The extraction was costing pattern.
  * The rmse *gain* is smaller on published tables (14% against 22%), and that is the expected
    direction given the audit: extraction finds 119 peaks where papers printed 267, so it reads
    more voxels as silent, shrinks harder, and scores better on an rmse dominated by the 9204
    near-zero voxels while over-shrinking the top (-0.065 against +0.005).

The docstring now carries both rows and says the published one is the one to quote. This is the
first time the central real-data claim has rested on coordinates a human transcribed from a
paper rather than on my own proxy for them.

## The implemented derivatives check out: 23 symbolic claims, none failing

jdkent's second example verifies the *implementation's* closed forms against sympy rather than
properties of the model, which is the higher-value kind. CBES has three such hand-derived forms
and none had ever been checked. `proofs/censored_information.py` does it:

  * `_censoring_terms`' `first` and `second` are `dS/dmu` and `d2S/dmu2`. Worth checking rather
    than trusting: `second = -(u phi(u) - l phi(l))/sigma^2` is not the chain rule applied to
    `first` on its face, and works only through `d(phi)/dx = -x phi`.
  * the `sign` generalisation really is the complement's derivative. One sign factor has to
    serve *both* derivatives of `1 - S`, which is true because `1 - S` differentiates to `-S'`
    and `-S''` alike -- asserted in the code, proved here, for both limbs.
  * `d2_over_prob - score**2` is `d2/dmu2 log p`, i.e. the code's curvature term really is
    `p''/p - (p'/p)^2` with the sign carried correctly through the reported limb.
  * an image pair's `(g - mu) * precision` and `-precision`.
  * **all three observed-information blocks.** `i_mu = -w(r h + r(1-r) s^2)`,
    `cross = -w s r(1-r)/(pi(1-pi))` and `i_pi = w(r/pi - (1-r)/(1-pi))^2` each equal minus the
    corresponding second derivative of `w log(pi a + (1-pi) b)`, with `a` left as an abstract
    function of `mu` so the check is not about a particular density. Anti-vacuity asserts that
    none of the three collapses to zero for a concrete `a`.

All 11 pass, alongside the 12 in `censoring.py`. So the estimator's algebra is right; what was
wrong this session was a *limit* taken outside these formulas (the `1 - r` over `1 - pi`
business) and the harness around them, not the derivatives themselves. That is worth knowing
precisely, because it says the remaining problems are specification and measurement rather than
calculus.

### What is left that algebra could still settle

Ordered by whether it would change what gets built:

1. **Why the local-maximum correction comes out with the sign backwards.** The model uses
   `P(|g| >= c)` for a report; a paper reports a voxel only if it cleared `c` *and* was a local
   maximum. That is the one specification error now localised -- over-statement 1.00, 1.78,
   1.08, 1.13 at true g of 0.2 to 0.8, non-monotone so no constant reweighting absorbs it -- and
   #74 says it, not identifiability, is what limits the magnitude. **But the obvious fix is
   already rejected empirically:** the RFT expected-maxima route was tried and its density has
   the sign backwards at a signal peak (#63), and no constant `q` helps (8.9 to 12.2% error as
   `q` falls from 1.00 to 0.56). So the algebraic question is not "write down the correction" --
   it is *why* the smooth-field density points the wrong way here, which would either rescue the
   route or close it for good. Doing it symbolically is the only way left to tell those apart,
   since the numerics have already said "no" without saying why.
2. **The identifiability condition properly stated.** I have the two-configuration case. What
   would be useful is the rank of the Jacobian of `k` distinct `(sigma, c)` pairs in
   `(pi, mu)` -- i.e. how many genuinely distinct study configurations are needed, and how the
   condition number degrades as the configurations cluster. That turns "spread helps" into a
   number a user could compute from their own collection before running anything.
3. **The EM's fixed point versus the profile maximum.** The profile interval revealed the
   estimate is stabilised by the start point and `max_iter` rather than by the data. Whether the
   EM's fixed point coincides with the profile maximiser where the likelihood is flat is a
   question about the map's contraction, and it decides whether `g` is an estimator or a
   regularised summary.
4. **`Var(pi mu)` and the determinant.** Two lines, explains #70 mechanically, worth having
   written down next to the delta-method code.

Items 1 and 3 are the ones that would change the estimator.

## Why the maxima density has the sign backwards: a term it does not contain

The last specification error, and the numerics had said "no" without saying why. `#63` recorded
that the RFT expected-maxima route made every focus worse and that no constant `q` helped. The
derivation (`proofs/reporting_probability.py`, 4 claims) locates the reason.

Write the observed field as `Z = m + e`, `e` smooth stationary Gaussian of mean zero, `m` the
signal. A local maximum at the origin needs `Z'(0) = 0` and `Z''(0) < 0`. At a *signal peak*:

  * `Z'(0) = m'(0) + e'(0) = e'(0)`, because `m'(0) = 0` by definition of a peak. The
    first-order condition says exactly what it says under the null.
  * `Z''(0) = m''(0) + e''(0) = -kappa + e''(0)` with `kappa = -m''(0) > 0`. So

        P(Z''(0) < 0) = P(e''(0) < kappa) = Phi(kappa / sigma_2),

    where `sigma_2 = sd(e''(0))`.

**The zero-mean field is the case `kappa = 0`, where that probability is exactly 1/2.** And
`Phi(kappa/sigma_2)` is strictly increasing in `kappa`, so a null-field density is not merely
*different* at a signal peak -- it is a **lower bound**. The understatement factor is
`2 Phi(kappa/sigma_2)`, equal to 1 at `kappa = 0` and tending to **2** as the peak sharpens.

So: exactly one of the two conditions defining a local maximum picks up the signal, and it is the
one the standard density pins at one half. Applying that density at a signal peak understates the
chance of a maximum by up to a factor of two, understates it most where the effect is sharpest,
and therefore supplies a correction that is too large precisely where the effect is -- which is
the direction the measurements found, and why the route looked like it had the sign reversed.

**What this does and does not settle.** It gives the direction of the error and a factor-of-two
bound, and it says the fix is not a different density but the *same* density with the signal's
own curvature in it -- a quantity the model does not carry and which is not recoverable from a
coordinate table. It does **not** derive the non-monotone pattern (1.00, 1.78, 1.08, 1.13), which
would need the interaction with the height threshold; that is stated in the proof rather than
glossed.

So the honest status of the reported limb: the error is understood, bounded, and **not fixable
from tables alone**, because the correction needs a per-study peak sharpness that papers do not
report. That is a better place to leave it than "the RFT route has the sign backwards", which was
true and unexplained.

One process note. A fifth claim in the first draft was `simplify(m'(0)) == 0` at a peak, which is
a tautology dressed as a proof -- it would have printed `[ok]` and established nothing. Removed.
jdkent's examples guard against exactly this with anti-vacuity asserts, and it is worth copying
the habit rather than only the format.

## Profiling: the win was in the image loader, and the algebra had nothing to give

jdkent asks for profiling and for faster algebraic strategies with equivalent output. Profiled a
whole-brain fit at 4 mm, 20 studies, 2 images.

**Without the null, 2.0 s total:** `_censoring_terms` 0.62 s (32%), **`gc.collect` 0.43 s (22%)**,
`_working_sets` 0.19 s, `_update_prevalence` 0.13 s.

**With a permutation null at only 40 iterations, 58.8 s total:** `_censoring_terms` **34.5 s
(59%)** over 1647 calls, `_update_prevalence` 8.5 s, `_mu_derivatives` 4.4 s. At the default
`n_iters = 1000` that is most of the wall clock, since each permutation runs a full EM.

### The loader: 61.7x, and it was free

The `gc.collect` is nilearn's. `safe_get_data` runs a full collection on every call, and
`masker.transform` reaches it -- NiMARE's own conftest says so, and `meta/cbma/utils.py` already
avoids it elsewhere with direct boolean indexing. `_load_image_studies` was not:

    masker.transform (aligned image, 29398 voxels) : 136.72 ms
    boolean index instead                          :   2.22 ms     61.7x

Four reads at two image studies is 0.55 s of a 2.0 s fit, and it grows linearly with the image
count. Shipped, with the fast path declined when the grids differ or the masker carries
standardisation, detrending, smoothing or a target grid -- those change the values, and skipping
them for speed would be a different estimator. Values agree to 1.8e-7: nilearn resamples even
onto an identical grid, so the fast path skips an interpolation rather than reproducing it. Below
the float32 the maps are stored in, but **not bit-identical**, and said so in the commit.

### The algebra had nothing to give, and the docstring was wrong about why

`_censoring_terms`' docstring said the cost is "spread over four memory-bound kernels, so what
pays is removing a pass rather than speeding one up". Measured on 600,000 pairs:

| kernel | ms |
|---|---|
| `ndtr(upper)` | 10.66 |
| `ndtr(lower)` | 8.13 |
| `exp` density, upper | 2.98 |
| `exp` density, lower | 3.01 |
| one multiply pass | 0.40 |
| one divide pass | 0.44 |
| `np.where(sign > 0, ...)` | 2.60 |

**45% of it is two normal CDFs**, which no rearrangement of the surrounding arithmetic reaches.
Three were tried and benchmarked against the shipped version for speed *and* agreement:

| variant | ms | speedup | max relative drift |
|---|---|---|---|
| shipped | 41.82 | 1.00 | -- |
| one reciprocal for two divisions | 41.68 | 1.00 | 2e-16 |
| `second` from `first` by the identity | 38.32 | **1.09** | 5e-14 |
| both | 39.68 | 1.05 | 5e-14 |

The identity is `u phi(u) - l phi(l) = u(phi(u) - phi(l)) + k phi(l)` with `l = u - k`, so
`second = (u/sigma) * first - k phi(l)/sigma^2` reuses `first` instead of forming the two
products. Correct, and worth 9%. **Not taken:** 9% for a rewrite of the estimator's most
delicate function, with drift, is a bad trade.

### Two further ideas, checked and rejected

  * **Drop the lower tail's CDF.** Its median contribution is 2e-5 of the silent probability,
    which sounds negligible -- but its maximum is 0.30, and `phi(lower)` exceeds 1% of the
    score's numerator for **38%** of pairs. It is not droppable.
  * **Cache the silent working set across permutations.** The null holds the silence pattern
    fixed, so the indicator pairs look permutation-invariant and reusable. They are not:
    `inv_sigma_ind = 1/sqrt(null_var + tau2[voxel])` and `tau2` is re-estimated from each
    permutation's image values. Checked before implementing, which is the only reason it did
    not become a wrong-answer bug.

So the honest summary: one large real win in the loader, nothing available in the censoring
algebra, and the lever for the permutation null is `n_cores` -- it is a thousand independent
fits and that is inherent to the design, not an inefficiency.

## Every error rate on record was measured under a null that no longer exists

Found while rewriting the PR description, which still advertised the old design in full.

The familywise figures in that description -- the 0.070 and 0.060 cells, the **0.180** at two foci
per study, and the whole max-statistic guard calibration above -- were measured under the
*arrangement* null, which shuffled reported effect sizes among the foci of an analysis. The
redesign reads no magnitude from a table, so that null was deleted. The shipped one scrambles each
image study's values among that study's own voxels and holds the silence pattern fixed.

**Different mechanism, so the numbers do not transfer, and the defect they exposed was a property
of the deleted null.** The two-foci pathology came from swapping two values inside a study being a
tiny perturbation of the max statistic; there is nothing to swap now, because the tables carry no
values. Whether the image permutation has its own pathology is simply unknown.

`experiments/fpr_after_redesign.py` measures it fresh: global null, nothing with an effect
anywhere, thresholds varying study to study as in a literature search, three arms (20 studies with
2 images, 20 with 1, 12 with 2), reporting the two quantities separately because they are easy to
confuse -- the mean *share of voxels* under an uncorrected p < 0.05, and the share of
*simulations* with any voxel surviving FWE.

A two-simulation pilot put the uncorrected rate at 0.0501, 0.0379 and 0.0464, which is the right
neighbourhood; the familywise rate needs the full run and is not quotable from two draws. The 40
simulation run is in flight.

**The PR description now says this is in progress rather than carrying the stale figures**, which
is the part that mattered: a reviewer reading them would have been assessing an estimator that no
longer exists.

## SDM has never run on NeuroVault, and on the cached collection it cannot

jdkent asked whether a NeuroVault-based studyset had been included in the SDM-PSI comparison. It
had not: SDM had been run on HCP (2 images) and on NIDM pain (1 image, published tables) only.
Writing `animal_sdm_vs_cbes.py` to close that produced a more useful answer than the comparison
would have.

**All six splits were skipped for want of studies.** Diagnosing it rather than tuning around it:

```
11 animal studies; peaks per study by scheme
   study  max|z|  cluster    fdr    fwe
    8836    3.89        0      0      0
    8838    0.70        0      0      0
    ...
   TOTAL                0      0      0
```

**Zero peaks, every scheme, every study.** Seven of eleven have a whole-map maximum below 2.

Suspecting the harness first, the derived z maps were checked against NeuroVault's own raw t maps
for the same collections: 5.23 -> 4.23, 4.63 -> 3.80, 1.73 -> 1.61, the shrinkage being the t-to-z
map at df ~ 19. The conversion is sound. **The collection is genuinely weak** -- no study in it
would have printed a single focus in a paper -- so it cannot support *any* coordinate-based
method, and the absence of an SDM arm there is a property of the data.

### Which reopens #31, the NeuroVault validation recorded as complete

`convergence_second_collection.py` extracts with the same `report_peaks` call, so it should
produce nothing now. Run: **zero splits for every arm** -- `CBES g`, `g_marginal`, `prevalence`,
ALE, MKDA density and KDA all blank. It then prints its concluding paragraph about the pain
ordering underneath, which is exactly how an empty result gets read as a finding.

The reason the validation ever produced numbers is in the older script it grew from.
`second_collection.py` thresholds with `ImagesToCoordinates(z_threshold=U)` at `U = 3.2905` --
**a bare uncorrected height, no multiplicity correction and no extent test**, which the standing
protocol forbids in as many words. Under that rule exactly 3 of 11 studies report at all, one of
them a single voxel:

```
study 8836: 74 suprathreshold voxels (max |z| 4.23)
study 8893:  1 suprathreshold voxel  (max |z| 3.44)
study 8962: 31 suprathreshold voxels (max |z| 3.80)
```

So the one NeuroVault real-data validation in this program **rests on an extraction the protocol
forbids, on three studies, and cannot be redone compliantly on this collection.** #31 should not
be counted as a NeuroVault validation. What remains is pain -- the only bed with transcribed
tables -- and HCP, whose tables are extracted but whose signal is strong enough to survive
correction.

To get a real NeuroVault arm would need a different collection: studies whose maps actually clear
a corrected threshold. That is a fetch over the network and a new selection step, not a rerun.
