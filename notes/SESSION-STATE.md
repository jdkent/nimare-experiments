# Where CBES stands — estimands and uncertainty

## Current state, latest session

Everything below this section predates it and several of its numbers have since been retracted;
the retractions are listed here and the detail is at the end of `notes/cbes-evidence.md`.

### Shipped to the PR

  aaf1557  two bugs found by a brute-force oracle: the EM retiring at its start value where the
           reported limb makes the curvature positive, and selection_model="none" profiling out a
           prevalence it never fitted
  6d4252f  cite MetaNSUE, SDM-PSI and ES-SDM as antecedents
  fca72c9  document max_iter as shrinkage; drop a null_method explanation whose sign was backwards
  4bd49c9  se is not reliably conservative; which way it errs depends on the prevalence
  79d914f  retract the family-wise inflation figures
  f55b919  call the small-collection family-wise rate unverified
  f33b630  the profile bound's finiteness is the diagnostic; the regime usually cannot be tested
  fd6d1fb  correct the prior-art claim: SDM-PSI maximises the same censored likelihood
  275b5fc  cover the generator's preconditions and its per-study variation

### Proofs, in `proofs/`, each with a calibration check

  mixture_identification            the reporting indicator's information matrix is rank 1, so
                                    tables cannot separate mu from pi however many there are;
                                    for unlike studies det I = pi^2 sum (u_j v_k - u_k v_j)^2,
                                    which reproduces #73's measured 35x sample-size advantage
  roughness_cannot_be_liberal       roughening a null raises its maximum, so destroying spatial
                                    autocorrelation makes a max-statistic test conservative
  count_versus_fitted_prevalence    P(report) = 1 - s_0 - pi(s_a - s_0), so the naive count ranks
                                    perfectly wherever the magnitude is flat
  boundary_test_for_full_prevalence testing pi = 1 is a boundary problem; the level-0.05 cut is
                                    2.7055, not 3.8415
  early_stopping_is_shrinkage       t EM iterations give (1 - rho^t) theta* + rho^t theta_0, a
                                    ridge whose weight grows with the missing information

### Retracted this session, all of them mine

  * the family-wise inflation. Every 0.150 came from 40 simulations, where the binomial se is
    0.034. Three estimates (0.150 at 40, 0.055 at 200 pre-fix, 0.100 at 100 post-fix) do not
    separate. **Unverified below about 20 studies**, neither sound nor broken.
  * "the censored likelihood beats imputation". The point estimates tie at 25,600 fits; only the
    interval differs, and partly through finite M.
  * "converging fixes the prevalence". It overshoots: 0.541 to 0.723 against a true 0.6. The
    ordinality caveat stands.
  * "se is conservative". 1.25 to 1.52 at prevalence 1, 0.80 to 0.83 at 0.6 — anti-conservative
    in the regime the estimator is for.
  * the identified_share diagnostic, withdrawn before shipping when simulation contradicted it.
  * "the mu bias is O(1/k)". It is a fixed 22% of the sd, a boundary effect, worth 5% of MSE.
  * "SDM-PSI imputes rather than using a censored likelihood". It uses one, citing Tobin.
  * "a defaulted sample size misreports g". Flattening a real 3.56x spread left g unharmed.

Each failed the same way: a ratio looked constant or a difference looked large, and I had not
checked the denominator. **No rate goes into documentation without its standard error, and 40
replicates sizes a run rather than producing a result.**

### Failure modes as they now stand

  1. The estimand collapses where every study has the effect, and the collection cannot tell you
     which regime it is in — the boundary test needs about 20 images for 0.62 power at a true 0.6.
  2. The interval is wrong in both directions and no setting fixes it; read it only where
     g_lower/g_upper are finite, where coverage is 0.97 against 0.87 elsewhere.
  3. The default does not converge, and that early stopping is buying real accuracy on mu.
  4. Coordinate tables cannot separate mu from pi. Proven, not measured.
  5. Reported peak heights are discarded, so the report limb is over-stated by up to a factor of
     two, and that is not correctable from tables.

### Open

  #31 NeuroVault validation not protocol-compliant. #40/#41 point-process estimator.
  #58 residual study-count drift. #62 dof for a censored-likelihood se. #66 q(mu) routes blocked.
  #86 match SDM-PSI's post-fit imputation. #88 re-measure every quoted rate at a stated se.

Wanted, and not obtainable here: Costafreda 2012 (doi:10.1016/j.jneumeth.2012.07.016), which may
already contain the rank-1 result; the MetaNSUE methods paper (doi:10.1177/0962280218811349);
Schnedler 2005 on censored random vectors.

---

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

## The coverage table, scored against the interval the documentation recommends

Truth 0.800, prevalence 1 at every site, calibrated reporting regime, **100 replications per
arm**, studies reporting a genuine *t*, every coordinate table produced by a cluster-forming
threshold with no cap on the number of foci. `cov(z)` is `g +/- 1.96 se`; **`cov(t)` is the
interval the class docstring recommends** -- a *t* on each fit's own `dof`, taken per replication
because the critical value is nonlinear in `dof`. `width` is that interval's half-width as a
fraction of the effect, because coverage without width is not a measurement.

Three earlier versions are superseded: one with a mis-calibrated reporting regime, one that ran
every arm at `peak_bias=None` rather than the documented configuration, and one that scored only
the normal interval. Raw log in `results/coverage_both_intervals_21arm.log`.

| studies | images | configuration | bias | se/sd | dof | cov(z) | cov(t) | width |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 12 | 0 | — | +0.255 | 2.14 | 4.5 | 0.75 | **0.98** | 0.57 |
| 12 | 0 | fwhm 16 | +0.248 | 1.70 | 7.8 | 0.28 | 0.52 | 0.32 |
| 12 | 0 | fwhm 24 | +0.249 | 1.50 | 9.7 | 0.05 | **0.13** | 0.24 |
| 12 | 2 | calibrated | −0.038 | 1.32 | 4.4 | 0.99 | 1.00 | 0.36 |
| 12 | 2 | `peak_bias=None` | +0.127 | 1.36 | 4.4 | 0.87 | 0.99 | 0.45 |
| 12 | 6 | calibrated | −0.026 | 1.26 | 7.0 | 0.98 | 1.00 | 0.25 |
| 12 | 6 | `peak_bias=None` | +0.027 | 1.27 | 7.0 | 0.99 | 1.00 | 0.27 |
| 12 | 12 | calibrated | −0.018 | 1.10 | 11.0 | 0.94 | 0.97 | 0.18 |
| 24 | 0 | — | +0.246 | 2.01 | 9.1 | 0.35 | **0.56** | 0.33 |
| 24 | 0 | fwhm 16 | +0.243 | 1.79 | 16.2 | 0.03 | 0.04 | 0.20 |
| 24 | 0 | fwhm 24 | +0.244 | 1.59 | 20.2 | 0.00 | **0.00** | 0.16 |
| 24 | 2 | calibrated | −0.064 | 1.35 | 8.6 | 0.91 | 0.96 | 0.22 |
| 24 | 2 | `peak_bias=None` | +0.163 | 1.48 | 8.6 | 0.58 | 0.74 | 0.29 |
| 24 | 6 | calibrated | −0.049 | 1.41 | 10.8 | 0.92 | 0.96 | 0.19 |
| 24 | 6 | `peak_bias=None` | +0.073 | 1.38 | 10.8 | 0.89 | 0.93 | 0.22 |
| 24 | 24 | calibrated | −0.022 | 1.11 | 23.0 | 0.97 | 0.97 | 0.12 |
| 12 | 0 | τ 0.3 | +0.268 | 1.70 | 3.6 | 0.84 | 0.94 | **0.91** |
| 12 | 6 | τ 0.3 | −0.048 | 1.06 | 6.7 | 0.92 | 0.94 | 0.43 |
| 12 | 12 | τ 0.3 | −0.015 | 1.12 | 11.0 | 0.93 | 0.96 | 0.30 |

(The `peak_bias=None` rows at 12-of-12 and 24-of-24 are identical to the calibrated ones to every
digit, which is the behavioural check on the all-donor calibration fix rather than a measurement.)

**The headline, stated for the interval a user following the documentation would build.** With
image donors the interval is usable: coverage 0.96 to 1.00 at half-widths of 0.12 to 0.36 of the
effect, across every image share and study count measured. Coordinates-only it is not, but not for
the reason I said earlier today. It *covers* at twelve studies (0.98) -- at a half-width of 0.57
of the effect, an interval running [0.60, 1.51] for a truth of 0.80. At twenty-four studies it
drops to 0.56. And under heterogeneity it covers 0.94 at a half-width of **0.91 of the effect**,
which is the clearest case in the table of coverage meaning nothing on its own.

**What actually discriminates is bias and width, and that is all.** The bias-to-width model
predicts `cov(t)` across all twenty-one arms to a mean absolute error of **0.024** (0.036 for
`cov(z)`), so there is no residual behaviour to explain. `se/sd` never falls below 1.06, so no
interval is too narrow for the estimator's own variability and every shortfall is bias.

**Widening the kernel is the one thing that breaks it outright**, and it breaks it twice over:
`se` falls *and* `n_eff` rises toward the study count so the critical value shrinks too. Coverage
runs 0.98 / 0.52 / 0.13 across fwhm 10 / 16 / 24 at twelve studies and 0.56 / 0.04 / 0.00 at
twenty-four, on a bias that does not move. A wide kernel was measured on real data to improve both
the accuracy of the map and its extent. One default cannot serve both.

**Neither `peak_bias` setting dominates and they fail differently.** The calibrated one is stable
across the image share (−0.018 to −0.064 from f = 0.08 to 1.00); the variant swings with it
(+0.163 at 2 of 24 down to −0.018 at 12 of 12). So the variant wins at 6 of 12 and loses badly at
2 of 24, and stability is the property worth having because the share is not the analyst's choice.

## Ledger: what stands, what was retracted

I corrected myself a lot today. Two root causes account for most of it — a simulator that reported
the wrong kind of statistic, and characterising the estimator in a configuration its own
documentation warns against — and both were in the input path rather than in the model.

### Stands

| finding | evidence |
| --- | --- |
| The pooling and the observed-information SE are correct | all-donor arms: `se/sd` 1.10–1.11, coverage 0.94–0.97 against nominal 0.95. Their −0.018 to −0.022 residual is now traced to inverse-variance weighting with a draw-dependent Hedges variance, not to the estimator: a draw-independent variance gives −0.002 and +0.001 |
| Coverage is a deterministic function of the bias-to-width ratio | a shifted normal on (bias, `se`, `sd`) predicts all sixteen τ=0 arms to a mean absolute error of 0.034, over measured coverage from 0.00 to 0.99. Supersedes the earlier `se/RMSE` vs `se/sd` framing: `se/sd` ranks arms backwards because it omits the bias, which is the whole story |
| Coverage degrades as studies accumulate under a fixed bias | 0.75 → 0.35 coordinates-only on doubling the studies (100 reps). At a *fixed* donor weight share it is slower: 0.97 / 0.97 / 0.72 at 12 / 24 / 48 studies, because `se/sd` rises (1.27 / 1.44 / 1.53) and offsets the shrinking interval |
| `prevalence` = reporting fraction inflated by explicable silence | fitted π runs 1.000 → 0.500 as the reported value moves from just above the cut to far above it, with the naive fraction fixed at 0.500 |
| Threshold inference is badly wrong under cluster-extent reporting | `study-min` infers z = 4.0 against a true forming cut of 3.1; prevalence error 0.201 against 0.008 for a fixed constant |
| Ordering within one map is unreliable | exactly right in 19% of maps at a strong effect, 6% at a weak one; threshold-independent |
| A naive count beats fitted `prevalence` on ordering, loses on level | 82% vs 40% exact ordering; bias 0.148 vs 0.075 |
| `g_marginal` beats convergence maps on CBES's support; `g` does not | +0.111 AUC, p 0.001 for the product; +0.001, p 0.980 for `g` |
| The default kernel is too narrow against an IBMA-like reference | fwhm 16 covers 64% of the truth's top decile against 34%, and improves AUC on already-covered voxels too |
| Kernel width trades the map against the interval | bias flat (+0.255 / +0.248 / +0.249), `se` halves, coverage 0.75 → 0.28 → 0.05 across fwhm 10 → 16 → 24 at twelve studies and 0.35 → 0.03 → 0.00 at twenty-four |
| The documented configuration's residual is small and *stable* across the image share | −0.038 / −0.026 / −0.018 at 2 / 6 / 12 of twelve studies and −0.064 / −0.049 / −0.022 at 2 / 6 / 24 of twenty-four, against +0.163 to −0.018 for `peak_bias=None` over the same range. Stability is the property worth having, since the share is not the analyst's choice |
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
| The interval has no validated operating point | with the documented configuration and two donors it covers 0.91–0.99 across every image share measured |
| The calibrated configuration overshoots twice what dilution predicts | about −0.008 of it is specific to the calibration; −0.005 is a drift present with no images at all |
| `g_relative`'s interval covers 0.00 | my test, not the estimator: I normalised the truth over 27,000 mostly-zero voxels and the fit over its smoothed support, so the two quantiles differ sevenfold. `g_relative` has no external truth to be scored against |
| Coverage degrades as a collection grows (citing the two-donor rows) | true, but those rows also halved the donors' weight share. At fixed share it is 0.97 / 0.97 / 0.72 at 12 / 24 / 48 studies — real, and slower than I implied |

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

Ordered by how much turns on the answer. Task numbers refer to the tracked list.

**1. #62 — the degrees of freedom for a censored-likelihood standard error.** `n_eff` counts only
studies whose foci reach a voxel; the likelihood also uses the studies that stayed silent, whose
silence is what bounds the effect. So the same studies are credited with information in the
numerator and denied it in the denominator: at a well-covered voxel the roster is 10 studies,
`n_eff` is 4.76, `dof` is 3.76, and the critical value is 2.87 where the roster would give 2.26.
Worth ~15-20% of interval width. **It is not a fix for the calibration** -- projected across
sixteen arms, coverage under a roster `dof` still spans 0.00 to 0.99. Options: the roster
(simple, overstates the other way), a Kish count weighting silences by their censoring
information (matches the likelihood, but the weight is parameter-dependent), a profile-likelihood
interval (needs no `dof`, fixes the `dof = 0` case where the recipe currently yields `nan`, most
work), or leave and document, which is the current state.

**2. #53 — should `threshold="study-min"` stop being the default?** It costs 0.20 of prevalence
accuracy on cluster-extent tables where a fixed constant costs 0.008. Newly measured: it costs
only about **0.044 of bias**, so the decision turns almost entirely on `prevalence`, not on `g`.

**3. #56 — should the default `fwhm` move off 10 mm?** Now a sharp trade rather than a preference.
Widening improves the map on real data and **destroys the interval**: coverage 0.98 / 0.52 / 0.13
across 10 / 16 / 24 mm at twelve studies, on a bias that does not move, because `se` falls *and*
`n_eff` rises toward the study count so the critical value shrinks too. A single default cannot
serve both goals.

**4. #61 — the max-statistic guard's threshold.** The guard fires on 78% of fits in the bad cell
and takes the overall familywise rate from 0.150-0.180 to 0.050 -- but that is
`(1 - 0.783) x 0.231` and not error control. Among fits it *passes*, rejection is 0.231 (exact
binomial p = 0.0245). The six-foci control is clean, so tightening is cheap. Which of its two
conditions to move is being measured.

**5. #57 — recast the mixed-collection guidance as weight share.** `r` is between 3.0 and 4.6
(measured out of sample), so one image donor carries three to five coordinate studies' worth of
pooling weight. Adding coordinate-only studies to a mixed collection therefore moves the
magnitude *away* from the truth while narrowing the interval. Shipped to the docstring; flagged
here in case you want it said differently.

**Smaller, and unchanged:**

- **#51** The naive count vs fitted `prevalence`: stand behind the level, or emit both?
- **#47** The images-only refusal is inconsistent with the all-donor path, and its stated reason
  ("no foci to permute") is false since `_permute_image_values` exists. Evidence favours allowing.
- **#49** Emit a per-voxel window-of-detectability diagnostic?
- **#48** Should a simulation-derived bias number go in a warning, or only the mechanism?
- **#43** How to present scale uncertainty for `g_absolute` — second interval, combined, or
  documented multiplication?
- **#36** Where `g` comes from when a collection has images.

**And one thing I cannot do for you:** `nimare-experiments` has **no git remote**. Every note and
script in it exists only in this container, which is reclaimed after inactivity. I am scoped to
`neurostuff/nimare` and cannot push elsewhere, and per your instruction this material does not
belong in the PR. I have sent you this file and a tarball of `notes/`, `PROTOCOL.md`, `results/`
and all the experiment scripts, but the repository needs a remote.

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
