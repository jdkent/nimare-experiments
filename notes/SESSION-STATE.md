# Where CBES stands — estimands and uncertainty

Written for a check-in. Everything is measured in this repo; the PR carries only the algorithmic
conclusions. Detail in `notes/cbes-open-program.md`, `notes/cbes-first-principles.md`,
`notes/cbes-success-criteria.md`.

## Read this first: I had a bug that inflated every magnitude number for most of the session

My simulators built a study's statistic as `(truth + noise/sqrt(n)) * sqrt(n)` — a **normal
statistic with known variance**. CBES, correctly for real data, reads a reported statistic as a
**t on n − 1 degrees of freedom** and maps it back before converting to an effect size. In the far
tail, where every reported peak lives, that map is steep: at z = 5.33 with n = 30 it turns
d = 0.98 into g = 1.25. So I was handing the estimator a number that did not mean what it thought,
and about a third of the bias I attributed to it was mine.

```
                                  bias at the same foci
known-variance z declared as Z            +0.630
proper t declared as T                    +0.296
```

I then spent a long time explaining the inflated number — inventing a "pooling step" contribution
(it is +0.001) and then a "conversion convexity" one (which *is* this bug; the genuine Jensen term
is +0.043). Beds now generate a real t and assert their own convention before measuring anything.
The lesson is in `PROTOCOL.md`: name the one-minute check that the measurement is wrong, and run
it before proposing a mechanism.

## The definitive coverage table

Truth 0.800, prevalence 1 at every site, calibrated reporting regime (every study reports, 3-4
clusters each), **100 replications per arm**, studies reporting a genuine *t*, every coordinate
table produced by a cluster-forming threshold with no cap on the number of foci. `half/truth` is
the interval's half-width over the effect, because coverage without width is not a measurement.

Two earlier versions of this table are superseded. The first had a mis-calibrated reporting
regime (18% of studies reported anything). The second ran every arm at `peak_bias=None` rather
than the configuration the docstring recommends, which is how a whole session's magnitude numbers
came to describe a variant. **This one runs the documented configuration
(`peak_bias="per-study"`, `peak_bias_scale="images"`) as the primary arm set, with
`peak_bias=None` kept alongside as a labelled variant so the difference is visible.**

| studies | images | configuration | bias | se/sd | coverage | half/truth |
| --- | --- | --- | --- | --- | --- | --- |
| 12 | 0 | — | +0.255 | 2.14 | 0.75 | 0.42 |
| 12 | 0 | fwhm 16 | +0.248 | 1.70 | **0.28** | 0.27 |
| 12 | 0 | fwhm 24 | +0.249 | 1.50 | **0.05** | 0.21 |
| 12 | 2 | calibrated | −0.038 | 1.32 | 0.99 | 0.26 |
| 12 | 2 | `peak_bias=None` | +0.127 | 1.36 | 0.87 | 0.33 |
| 12 | 6 | calibrated | −0.026 | 1.26 | 0.98 | 0.20 |
| 12 | 6 | `peak_bias=None` | +0.027 | 1.27 | 0.99 | 0.22 |
| 12 | 12 | calibrated | −0.018 | 1.10 | 0.94 | 0.16 |
| 24 | 0 | — | +0.246 | 2.01 | **0.35** | 0.29 |
| 24 | 0 | fwhm 16 | +0.243 | 1.79 | **0.03** | 0.19 |
| 24 | 0 | fwhm 24 | +0.244 | 1.59 | **0.00** | 0.15 |
| 24 | 2 | calibrated | −0.064 | 1.35 | 0.91 | 0.19 |
| 24 | 2 | `peak_bias=None` | +0.163 | 1.48 | **0.58** | 0.25 |
| 24 | 6 | calibrated | −0.049 | 1.41 | 0.92 | 0.17 |
| 24 | 6 | `peak_bias=None` | +0.073 | 1.38 | 0.89 | 0.20 |
| 24 | 24 | calibrated | −0.022 | 1.11 | 0.97 | 0.11 |
| 12 | 0 | τ 0.3 | +0.268 | 1.70 | 0.84 | 0.62 |
| 12 | 6 | τ 0.3 | −0.048 | 1.06 | 0.92 | 0.35 |
| 12 | 12 | τ 0.3 | −0.015 | 1.12 | 0.93 | 0.27 |

Six things this settles.

**1. The interval works with donors and does not without them.** Every calibrated arm with two or
more image donors covers 0.91 to 0.99 against a nominal 0.95. Every coordinates-only arm at the
default kernel covers 0.35 to 0.84. The bias is what separates them: about +0.25 (31% of the
effect) coordinates-only against −0.02 to −0.06 calibrated.

**2. `se/sd` never drops below 1.06.** No interval in the table is too narrow for the estimator's
own variability, so every coverage shortfall is bias. That kills "the standard error is wrong" as
an explanation for anything here.

**3. All of it is one ratio.** A shifted normal on (bias, se, sd) predicts coverage across the
sixteen τ=0 arms to a mean absolute error of 0.034, spanning measured coverage from 0.00 to 0.99.
The interval fails exactly when and as much as the bias-to-width ratio says. Earlier descriptions
of it as "too narrow" or "conservative" were restating `b/se` in words. On the τ=0.3 arm the
prediction is off by +0.10, at the edge of that band — the model degrades under heterogeneity.

**4. Coverage degrades as a collection grows, in every configuration.** `se` shrinks roughly as
`1/sqrt(studies)` and the bias does not shrink at all: 0.75 → 0.35 coordinates-only, 0.99 → 0.91
with two donors, 0.98 → 0.92 with six. A large coordinate-only collection gives a tighter
interval around the wrong value. This inverts the usual reassurance and is now in the docstring.

**5. Widening the kernel destroys the interval — which settles #56 as a genuine trade.** fwhm 16
takes coverage to 0.28 and fwhm 24 to 0.05 at twelve studies (0.03 and 0.00 at twenty-four),
because widening shrinks `se` by 40-50% and leaves the bias untouched. Widening was measured on
real data to improve both the accuracy of the map and its spatial extent. So the two goals point
opposite ways and a single default cannot serve both: **a wide kernel for the map, a narrow one
for the interval.**

**6. Neither `peak_bias` setting dominates, and they fail differently.** The calibrated one is
*stable* — −0.018 to −0.064 across image shares from 0.08 to 1.00. The variant *swings* with the
share — +0.163 at 2 of 24 down to −0.018 at 12 of 12 — so it wins at 6 of 12 (0.99 against 0.98,
where its +0.027 happens to be smaller than the calibrated −0.026) and loses badly at 2 of 24
(0.58 against 0.91). Stability across the share is the property worth having, because a real
collection's share is not something the analyst chooses.

And one oddity worth keeping: **heterogeneity improves coverage while making the point estimate
worse.** Coordinates-only at τ=0.3 has a larger bias (+0.268 against +0.255) and better coverage
(0.84 against 0.75), because it widens the interval more than it moves the estimate. The clearest
case in the whole table that coverage alone is not a metric.

Two arms are a behavioural check rather than a measurement: 12-of-12 and 24-of-24 give *identical*
numbers under `calibrated` and `peak_bias=None`, to every digit. That is what the all-donor
calibration fix should do — with every study imaged there is no coordinate value for a scale to
act on, so it returns 1.0 and the two paths become arithmetically the same fit.

## Ledger: what stands, what was retracted

I corrected myself a lot today. Two root causes account for most of it — a simulator that reported
the wrong kind of statistic, and characterising the estimator in a configuration its own
documentation warns against — and both were in the input path rather than in the model.

### Stands

| finding | evidence |
| --- | --- |
| The pooling and the observed-information SE are correct | all-donor arms: `se/sd` 1.00–1.11, coverage 0.94–0.97 against nominal 0.95 |
| `se/RMSE` predicts coverage; `se/sd` does not | 0.61 / 0.45 / 1.18 / 0.97 against coverage 0.72 / 0.28 / 0.99 / 0.97, while `se/sd` ranks them backwards |
| Coverage degrades as studies accumulate under a fixed bias | 0.72 at 12 studies → 0.28 at 24, coordinates only |
| `prevalence` = reporting fraction inflated by explicable silence | fitted π runs 1.000 → 0.500 as the reported value moves from just above the cut to far above it, with the naive fraction fixed at 0.500 |
| Threshold inference is badly wrong under cluster-extent reporting | `study-min` infers z = 4.0 against a true forming cut of 3.1; prevalence error 0.201 against 0.008 for a fixed constant |
| Ordering within one map is unreliable | exactly right in 19% of maps at a strong effect, 6% at a weak one; threshold-independent |
| A naive count beats fitted `prevalence` on ordering, loses on level | 82% vs 40% exact ordering; bias 0.148 vs 0.075 |
| `g_marginal` beats convergence maps on CBES's support; `g` does not | +0.111 AUC, p 0.001 for the product; +0.001, p 0.980 for `g` |
| The default kernel is too narrow against an IBMA-like reference | fwhm 16 covers 64% of the truth's top decile against 34%, and improves AUC on already-covered voxels too |
| Kernel width trades the map against the interval | bias flat, `se` halves, coverage 1.00 → 0.00 across 10 → 24 mm |
| The documented configuration pins the scale to ~5% regardless of donor count | −0.047 at two donors, −0.041 at six |
| `peak_information` is independently corroborated | +0.17 to +0.23 z excess in a favourable bed, agreeing with a separate measurement by another route |
| Prevalence/magnitude are separable only in a window of detectability | both recovered where `dD/dμ` is 1.9–2.2, neither where it is 0.0–0.5 |
| Power spread identifies the split | fitted `g` swings 0.387–0.603 under a fixed roster where the truth is constant 0.6; flat at 0.55–0.60 once `n` and `u` vary |

### Retracted

| I said | actually |
| --- | --- |
| Coordinates-only bias is +0.51 (63%) | +0.26 (32%); the rest was my simulator reporting a known-variance z where the estimator assumes a t |
| The pooling step adds as much bias as the winner's curse | it adds +0.001 |
| A convex conversion amplifies the curse by +0.32 | that step *was* the convention bug; the genuine Jensen term is +0.043 |
| Prevalence ordering is right in 50% / 12% of maps | 19% / 6% on the corrected convention — I was too generous |
| Deleting reported heights improves the map by +0.284 | it costs 0.268 with a correct threshold, and drops AUC to chance on real data. Feature request withdrawn |
| A couple of images does not distort the relative map | it does — U-shaped ratio spread 1.17× / 1.36× / 1.03×. My original prediction was right and I had retracted it on mis-instrumented evidence |
| Over the whole brain ALE is best and every CBES map worst | true only at the default kernel; at 16 mm CBES beats MKDA, at 24 mm it beats ALE |
| The coordinate channel is diluted, never corrected | false of the estimator; true only of `peak_bias=None`, which is what all my arms used |
| A shared scratch directory dropped arms from a run | it cannot have — the dropped set includes an arm that writes no images. Cause unresolved |

### The two habits these argue for

1. **Name the one-minute check that the measurement is wrong, and run it before proposing a
   mechanism.** Three separate errors were a surprising number, a mechanism sought in the
   estimator because that is what I was studying, and a one-minute check skipped.
2. **Make the documented default the first arm of every comparison.** Not "read the docs" — I had
   read the warning and written it into these notes before ignoring it.

## The biggest single finding: the default threshold inference is what miscalibrates `prevalence`

Under cluster-extent reporting — the commonest scheme in the literature — the smallest value a
study reports is its smallest cluster **maximum**, not its threshold. So `threshold="study-min"`
(the default) infers a cut ~0.9 z too high, and since `prevalence` is governed by where the
fitted magnitude sits relative to the assumed cut, it inflates badly. Mean absolute error over
four sites at true prevalence 0.25/0.50/0.75/1.00:

| threshold setting | median cut | prevalence | g | g_marginal |
| --- | --- | --- | --- | --- |
| `study-min` (default) | 4.015 | 0.201 | **0.173** | 0.270 |
| `pooled-min` | 3.532 | 0.056 | 0.211 | 0.166 |
| the true forming cut | 3.090 | 0.037 | 0.251 | **0.116** |
| the library constant 3.2905 | 3.291 | **0.008** | 0.235 | 0.132 |

With the library's own constant, prevalence comes back **0.248 / 0.501 / 0.752 / 0.972** against
a truth of 0.25 / 0.50 / 0.75 / 1.00. With the default: 0.482 / 0.860 / 0.961 / 0.998. So the
documented compression — "a true 0.25 comes back as 0.49 to 0.60" — is largely this default, not
the censored likelihood.

`g` prefers the default, for a bad reason: it's biased high by peak selection, and a cut assumed
too high makes the censoring term pull it down. That's one error cancelling another, and `g` stays
contaminated by prevalence at every setting. Recommendation with evidence attached in task #53:
supply a plausible fixed threshold; keep `study-min` only for voxelwise-height tables.

## The three findings I'd stand behind

**1. What `prevalence` actually computes.** Measured on the estimator's own EM with 6 reporting and
6 silent studies, so the naive reporting fraction is exactly 0.500, varying only the reported
height against a fixed cutoff: fitted prevalence runs 1.000, 0.964, 0.518, 0.501, 0.500 as the
reported value goes from just above the cut to far above it. So

> `prevalence` is the fraction of covering studies that reported, inflated by however much of the
> silence the censoring term can explain — and that is governed entirely by where the fitted
> magnitude sits relative to the reporting threshold.

This unifies four things previously recorded as separate puzzles: prevalence near 1 at strongly
reported sites; the floor and non-monotonicity at weak effects; why the likelihood beats a naive
count on bias but loses on ordering (it *is* the count times a noisy multiplicative correction,
which is the exact recipe for preserving a mean and destroying a ranking); and why `g_marginal`
works — it is magnitude times an empirical reporting probability, not two biases cancelling.

**2. A naive count beats the fitted prevalence at the reading the docs endorse.** Same collections,
both estimators. The naive count is `(studies with a focus within 15 mm) / (studies)`:

| | mean abs bias | rank corr | exact 4-site ordering |
| --- | --- | --- | --- |
| CBES prevalence, effect 0.8 | **0.061** | +0.820 | 42% |
| naive count, effect 0.8 | 0.102 | **+0.951** | **80%** |
| CBES prevalence, effect 0.5 | **0.129** | +0.570 | 22% |
| naive count, effect 0.5 | 0.414 | **+0.877** | **60%** |

(measured at the default threshold and on the old statistic convention; being re-run. The bias
column should improve for CBES with a supplied threshold — the ordering column will not, since
ordering is threshold-independent.)

The censored likelihood does what it was built for — it corrects the count's downward bias, by a
factor of three — and loses decisively on ordering, which is the only reading the documentation
endorses. That combination isn't tenable: either stand behind the level and retract "read it
ordinally", or emit the count too and say which to use for what. (Being re-run on the t
convention; the comparison is like-for-like either way since both read the same data.)

**3. ~~Reported heights are worse than useless.~~ RETRACTED — see the third correction below.**
The reported magnitudes carry little but they are not harmful, and they are doing most of the
localising on real data.

## A second retraction: the ordering numbers I shipped were too generous

I documented the prevalence ordering caveat with numbers from the buggy statistic convention. On
the corrected convention it is worse: the four-site ranking is exactly right in **19%** of maps at
a strong effect (not 50%) and **6%** at a weak one (not 12%), with rank correlations +0.76 and
+0.33 (not +0.86 and +0.59). The docstring has been updated. The ordering is *threshold*-
independent — identical to three decimals under both settings — because a change of cutoff
applies a roughly common inflation across voxels, so this stands on its own and is not fixed by
task #53.

## A third retraction: I recommended deleting the reported heights, and that was wrong

I found that flattening every reported height improved the map's correlation with the truth by
+0.284 (paired p < 0.001), called it the program's strongest result, and filed a feature request.
It was an artefact of a threshold interaction. Once every focus carries the same number,
`threshold="study-min"` infers a cutoff *equal to that number*, so every observation sits exactly
at its own censoring bound and the map's variation comes from censoring geometry rather than data.

```
                    threshold: study-min (inferred)   threshold: supplied 3.2905
as reported                    0.429                             0.297
study-flattened                0.638  (+0.209, p 0.000)          0.133  (-0.164, p 0.000)
all-flattened                  0.713  (+0.284, p 0.000)          0.029  (-0.268, p 0.000)
```

The sign flips. And the real-data check — NIDM pain, 10 splits, held-out reference, real threshold
handed in — agrees with the supplied-threshold column: flattening takes r from 0.230 to 0.04 and
the AUC for the truth's top decile from 0.653 to 0.48, i.e. chance, both at paired p < 0.001.

So the magnitudes carry essentially all of this fit's localising ability on real data. What
survives from the other measurements is only that they carry *little in absolute terms* — 0.055 z
per unit of true g at a 3.29 cut, 5–9% of the variance at a focus's own location. "Little" is not
"negative", and I conflated those. The feature request is withdrawn (task #54).

The protocol caught this: it mandates confirming on real data, and that is what reversed it.

## The kernel width, which turned out to matter more than anything else measured

Task #38 asked whether CBES localises better than ALE/MKDA against a held-out image reference.
The answer went through three versions, and the third is the one to keep.

| | brain covered | truth's top decile covered | AUC on covered | AUC whole mask |
| --- | --- | --- | --- | --- |
| fwhm 6 mm | 0.013 | 0.060 | 0.677 | 0.527 |
| **fwhm 10 mm (default)** | 0.090 | 0.336 | 0.757 | 0.642 |
| fwhm 16 mm | 0.263 | 0.636 | **0.798** | 0.745 |
| fwhm 24 mm | 0.608 | 0.879 | 0.780 | **0.782** |
| MKDA density | 1.000 | 1.000 | — | 0.667 |
| ALE | 1.000 | 1.000 | — | 0.753 |

At the default kernel CBES estimates 9% of the brain, captures a third of the truth's strongest
voxels, and **loses** to MKDA over the whole brain. At 16 mm it beats MKDA; at 24 mm it beats ALE.
And it is not a trade — the AUC on the *already-covered* voxels improves too (0.757 → 0.798), so a
wider kernel makes the estimate better where it existed and also extends it.

So the "CBES is the worst arm over the whole brain" finding I reported was a default-parameter
artefact. The qualification that survives: the optimum depends on how extensive the pooled effects
are, and this reference is a smooth field, so a wide kernel is rewarded for matching its extent.
10 mm is too narrow against an IBMA-like reference; whether 16 mm is right in general isn't
settled by one collection. Task #56.

Along the way I swept the wrong knob first — `coverage_radius`, which governs which studies count
as *silent*, not where estimates exist. The tell was a covered-share of 0.089 identical at every
radius from 8 to 45 mm while the AUC column moved. Rule added to `PROTOCOL.md`: confirm a
parameter moved the thing you're attributing to it before reading the sweep.

## Theory that held up

A coordinate observation is an **occupancy record with imperfect detection**: study *k* reports
near voxel *v* with probability `pi_v * D(mu_v, n_k, u_k)`. Fitting that likelihood on exact
records with the true detection function — a ceiling on what CBES could do — recovers *both*
factors well where the detection gradient is large and fails at both ends:

```
        true mu     0.00    0.20    0.40    0.60    0.80
      median mu   -0.325   0.208   0.425   0.611   0.769
      median pi    0.000   0.133   0.687   0.767   0.542
   dD/dmu at truth  0.000   0.341   1.876   2.195   0.546
```

So there is a **window of detectability** — roughly where `mu*sqrt(n)` is within ~1.5 of the
threshold. Below it nothing is detected and prevalence is unidentifiable; above it detection
saturates and prevalence absorbs the level. A map spans the window, which is why comparing two
voxels compares quantities identified to different degrees.

**Confirmed prediction:** a spread of study power widens the window. Under a fixed roster, fitted
`g` swings 0.481→0.734 across a prevalence sweep where the truth is a constant 0.6; once sample
sizes vary it is flat at 0.62–0.67, and the prevalence slope on the truth goes +0.534 → +1.019.
Real collections do have heterogeneous sample sizes, so the regime that breaks the separation is
not the common one — and it is now a checkable precondition.

**A constructive proposal (validated once):** report the **reporting probability** standardised to
a stated reference design — "in a study of 30 at p<0.001 corrected, the chance of a focus within
10 mm of this voxel is 0.42". It is the expectation of an observed Bernoulli, so it needs no scale
constant and no factoring; it sits in the strong (count) channel; and it answers the power
question directly. Estimated/true agreement is within 0.04 inside the collection's range of n and
degrades smoothly to 0.11 outside it.

## Shipped to the PR this session

- `se_marginal`: `g_marginal` — the one magnitude an IBMA also estimates — now has a standard
  error, by the delta method on the full 2×2 observed information, verified against a numerical
  Hessian of the estimator's own likelihood.
- `scale_interval_` is now a log-scale confidence interval for the scale, not the donors' sample
  range (which was 0.51× an honest interval at 2 donors and 4.32× at 20 — the error changed sign).
- The prevalence ordering claim corrected: reliable on average over many maps, not within one
  (exactly right in 50% of maps at a strong effect, 12% at a weak one), with the
  window-of-detectability reason stated.
- The degrees-of-freedom caveat: reading a reported z as a t on n−1 df is load-bearing, and the
  sensitivity grows with the height (1.08× at z=3.3, 1.38× at z=6). Published z-maps often have a
  higher effective df, so the likely direction is a further over-statement.

- The voxel-level family-wise correction is now withheld when the permutation distribution of
  the maximum barely moves (both a low distinct-value share and a low coefficient of variation).
  The arrangement-count guard passed a configuration that rejected at 0.150 against a nominal
  0.050, because two-focus studies admit 2^k arrangements while swapping two similar magnitudes
  moves the map's maximum hardly at all.

- **Three silent losses on the legacy `Dataset` conversion**, found by asking why a round-trip
  changed the calibrated scale. `nimare.io`'s supported set omitted `g`/`g_var`, so a collection
  carrying effect-size maps lost them at `to_dataset()`. Point values were dropped entirely --
  `Point.values` is a dict keyed by column name and the converter only handled the raw-JSON
  list-of-dicts shape, so `isinstance(..., dict)` rejected every one and every reported peak
  height vanished. And `_statistic_column` quietly returned "images carry the fit; the
  coordinates are unused" when no coordinate carried a statistic, which is the
  images-without-coordinates case the estimator refuses by name elsewhere -- one decision made
  two ways, with the quiet one reachable by accident. Any of the three turns a CBES fit into a
  different estimator with no error and no warning. A fourth difference is recorded rather than
  fixed: `to_dataset()` does not carry the mask, so the fit silently falls back to the default
  MNI template, which moved the calibrated scale 0.933 to 1.138.
- **Where the image channel's 2-3% downward bias comes from**, documented rather than corrected.
  Inverse-variance pooling with Hedges' variance `1/n + g^2/(2(n-1))` computes the weight from the
  observed effect, so a study that drew high gets less weight and the pooled estimate is pulled
  toward zero. A `g_var` map converted from a test statistic carries the term; one from a
  per-voxel mixed model does not; the estimator cannot tell which it was handed, so correcting it
  would be wrong as often as right.
- **The end-to-end coverage of the interval on `g`.** The only coverage number the docstring
  carried was measured where the model is exactly true (94.5-98.4%). What a user gets is in the
  table below, and the counterintuitive half is now stated: coverage *degrades* as a collection
  grows, because `se` shrinks and the bias does not.

CI green on every push. 104 tests pass.

## What needs your decision

- **#53** Should `threshold="study-min"` stop being the default? It costs 0.20 of prevalence
  accuracy on cluster-extent tables where a fixed constant costs 0.008.
- **#56** Should the default `fwhm` move off 10 mm? It costs both accuracy and the comparison
  against convergence estimators on the one collection that can be measured.
- **#51** The naive count vs fitted prevalence: stand behind the level, or emit both?
- **#47** The images-only refusal is inconsistent with the all-donor path, and its stated reason
  ("no foci to permute") is now false since `_permute_image_values` exists. Allow both, or refuse
  both? Evidence favours allowing.
- **#49** Emit a per-voxel window-of-detectability diagnostic?
- **#48** Should a simulation-derived bias number go in a warning, or only the mechanism?
- **#57** The image-*fraction* framing looks wrong now. With `r` between 3 and 4.6 (one image
  study carries three to five coordinate studies' worth of pooling weight, measured out of
  sample) what matters is the weight share `r*n_img / (r*n_img + n_coord)`. That falls when
  coordinate-only studies are added, so **adding coordinate studies to a mixed collection makes
  the magnitude worse while narrowing the interval** -- both pushing coverage down. Should the
  guidance be recast in those terms?
- **#43** How to present scale uncertainty for `g_absolute` — second interval, combined, or
  documented multiplication?
- **#36** Where `g` comes from when a collection has images.

## What the two root errors did and did not touch

Worth separating, because the retractions cluster in one part of the work and leave another
untouched.

**The statistic-convention error and the configuration error both hit the *magnitude accuracy*
numbers.** Every bias figure, the weight-share model, the scale-error ratios, the comparative
localisation margins — those are the ones that moved, some of them twice.

**They did not touch the estimand analysis**, which is what the mandate asked about, because that
work did not depend on either:

- `g` estimates a conditional magnitude and no image-based meta-analysis estimates it, so
  `g_absolute` cannot be validated against one even in principle. An argument about what the
  estimators target, not a measurement.
- `g_marginal` shares an IBMA's estimand, `π·μ`. Same.
- `prevalence` is the fraction of covering studies that reported, inflated by however much of the
  silence the censoring term can explain. Measured on the estimator's own EM at a single voxel
  with the reporting pattern held fixed — no simulated fields, no conversion, no `peak_bias`.
- The window of detectability, and that prevalence and magnitude are separable only inside it.
  Measured with an exact occupancy likelihood on exact detection records, which never passes
  through CBES at all.
- Reporting probability as an identified, standardisable estimand, and achievability as its value
  at a planned design. Derivation plus an exact-model check.
- `peak_information` is independently corroborated. It works on the z scale, where the convention
  error never applied.

So the picture of *what CBES can and cannot estimate* has held up through both corrections, while
the picture of *how accurately it estimates it* has been rebuilt twice. That is a reassuring
division in one way — the conceptual work was not resting on the broken harness — and a warning in
another: the numbers are the part that will end up in a paper, and they are the part that kept
being wrong.
