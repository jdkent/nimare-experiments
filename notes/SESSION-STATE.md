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

## The corrected coverage picture

Truth 0.800, prevalence 1 at every site, calibrated reporting regime (every study reports, 3–4
clusters each), 100 replications, studies reporting a genuine t. `half/truth` is the interval's
half-width over the effect.

| studies | images | bias | se/sd | coverage | half/truth |
| --- | --- | --- | --- | --- | --- |
| 12 | 0 | +0.259 | 2.19 | 0.72 | 0.40 |
| 12 | 0, `peak_bias='per-study'` | +0.255 | 2.14 | 0.75 | 0.42 |
| 12 | 2 | +0.127 | 1.36 | 0.87 | 0.33 |
| 12 | 6 | +0.027 | 1.27 | 0.99 | 0.22 |

(the 12-of-12, heterogeneity and 24-study arms were still running at the time of writing)

So: the coordinates-only magnitude runs about **32% high**, not 63%, and the interval covers 0.72
against a nominal 0.95 — a real failure, not a collapse. Two images halves the bias; six
essentially removes it. `se/sd` stays at or above 1 everywhere, so no interval is too narrow for
the estimator's own variability: every coverage shortfall is bias.

`peak_bias='per-study'` moves the bias from +0.259 to +0.255. That is consistent with what the
docstring already claims (it corrects the between-study part, not the common scale), but it puts a
number on it: the between-study part is negligible. Advice that stops at "use `peak_bias`" is
advice to do nothing.

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

CI green on every push. 100 tests pass.

## What needs your decision

- **#53** Should `threshold="study-min"` stop being the default? It costs 0.20 of prevalence
  accuracy on cluster-extent tables where a fixed constant costs 0.008.
- **#51** The naive count vs fitted prevalence: stand behind the level, or emit both?
- **#47** The images-only refusal is inconsistent with the all-donor path, and its stated reason
  ("no foci to permute") is now false since `_permute_image_values` exists. Allow both, or refuse
  both? Evidence favours allowing.
- **#49** Emit a per-voxel window-of-detectability diagnostic?
- **#48** Should a simulation-derived bias number go in a warning, or only the mechanism?
- **#43** How to present scale uncertainty for `g_absolute` — second interval, combined, or
  documented multiplication?
- **#36** Where `g` comes from when a collection has images.
