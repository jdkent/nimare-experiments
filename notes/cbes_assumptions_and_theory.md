# CBES: the assumptions, the derivation, and why it works

A companion to `effect_size_cbma.md`. That note records what was tried; this one states what is
being *assumed*, derives what follows, and places those assumptions next to the ones the rest of
neuroimaging meta-analysis already makes. The last section explains the whole thing without
equations.

---

## Part I. Notation and assumptions

### Notation

| symbol | meaning |
|---|---|
| $v$ | a voxel |
| $k = 1 \dots K$ | studies |
| $\theta_k(v)$ | study $k$'s true (population) standardized effect at $v$ |
| $\mu(v)$ | the meta-analytic mean effect at $v$ |
| $\tau^2(v)$ | between-study variance at $v$ |
| $\pi(v)$ | fraction of studies for which the effect at $v$ is genuinely non-zero |
| $N_k$ | study $k$'s sample size |
| $\hat g_{ki}, s^2_{ki}$ | Hedges' $g$ and its sampling variance at reported peak $i$ |
| $c_k$ | study $k$'s reporting threshold, on the $g$ scale |
| $w_k(v)$ | kernel weight, $1$ at the reported peak, decaying with distance |
| $R(v)$ | the set of studies that reported *something* near $v$ |

### The assumptions, stated sharply

**A1 — Statistic-to-effect-size is exact.**
A reported $t$ (or $z$) together with $N_k$ determines $\hat g$ and $s^2$:
$$\hat g = J(\nu)\,\frac{t}{\sqrt{N}}, \qquad s^2 = \frac{1}{N} + \frac{\hat g^2}{2N}$$
with $J$ the Hedges correction. This assumes the reported number is the *test statistic* of a
one- or two-sample comparison on $N_k$ subjects — not a contrast value, a percent signal change,
or a statistic from a model with nuisance structure the formula doesn't know about.

**A2 — Local homogeneity of the effect field.**
$\theta_k$ varies slowly relative to the kernel bandwidth $h$: for $\|v - x\| \lesssim h$,
$\theta_k(v) \approx \theta_k(x)$. This is what licenses treating a peak's value as an
observation about *neighbouring* voxels, and it is the assumption behind **local likelihood**
estimation (Tibshirani & Hastie 1987) rather than mere local averaging.

**A3 — Conditional exchangeability (random effects).**
Given that study $k$ has an effect at $v$, $\theta_k(v) \sim N(\mu(v), \tau^2(v))$, independently
across studies. Studies are interchangeable draws; none is privileged.

**A4 — Selection is by height, on the observed statistic.**
Study $k$ reports a peak near $v$ if and only if $|\hat g_k| > c_k$. Three things are being
assumed: that selection depends on the *observed* statistic (not the true effect), that the rule
is a *height* threshold, and that $c_k$ is the same across the brain within a study.

**A5 — Zero inflation.**
With probability $\pi(v)$ study $k$ has a real effect drawn per A3; with probability
$1 - \pi(v)$ it has *exactly none*. This is the substantive modelling choice, and §II.2 shows
what breaks without it.

**A6 — Full coverage.**
A study that reported nothing near $v$ nonetheless *examined* $v$. Silence is evidence, not
absence of data. Violated by ROI studies.

**A7 — Independence across studies.**
No shared subjects, no correlation induced by shared pipelines or labs.

**A8 — Spatial exchangeability under the null.**
Under $H_0$, reported coordinates are exchangeable with positions drawn uniformly in the mask.
This is precisely the null ALE and MKDA use, adopted here unchanged.

**A9 — The reported value is a truncated draw.**
The model treats a reported $\hat g$ as a draw from the sampling distribution *conditioned on
exceeding $c_k$*. **This one is known to be false**, and knowing exactly how it fails is what
makes the residual bias predictable — see §II.5.

---

## Part II. What follows

### II.1 The estimator is a local random-effects meta-analysis

Under A1–A3, at each voxel the inverse-variance weighted mean with kernel weights,

$$\hat\mu(v) = \frac{\sum_k W_k(v)\, \hat g_k}{\sum_k W_k(v)}, \qquad
W_k(v) = \frac{w_k(v)}{s_k^2 + \tau^2(v)}$$

is the weighted-likelihood estimator for $\mu(v)$. Setting every $w_k = 1$ recovers textbook
DerSimonian–Laird exactly — which is the precise sense in which **the spatial model is a
weighting scheme, not a separate algorithm**. Everything meta-analytic about CBES is standard;
everything spatial is in $w$.

### II.2 Why zero inflation is not optional

Consider the model *without* A5 — a plain Tobit, one $\mu$ per voxel, silence read as censoring.
A silent study contributes $\log P(|\hat g| < c_k \mid \mu)$, which is maximised at $\mu = 0$ and
decreasing in $|\mu|$. So at a voxel where most studies are silent, the censoring term pushes
$\hat\mu$ toward zero *without limit*, and it cannot distinguish two situations that are
scientifically opposite:

- a study whose effect is real but sub-threshold, and
- a study whose effect is genuinely zero.

Both produce silence; the Tobit must explain both as "small $\mu$". Empirically this drove
estimates to $-0.16$ and $-0.34$ where the truth was positive. The mixture in A5 gives silence a
second explanation, and the *relative* number of reporting and silent studies identifies which
one the data support. This is the same structure as an occupancy model in ecology: a species not
detected at a site may be absent, or present and missed, and only a mixture can tell them apart.

### II.3 The scale is provably non-identified

Take any $\alpha > 0$ and rescale the *entire problem*:
$$\hat g \to \alpha \hat g, \quad s^2 \to \alpha^2 s^2, \quad c_k \to \alpha c_k,
\quad \mu \to \alpha\mu, \quad \tau^2 \to \alpha^2\tau^2 .$$

- The censoring terms $P(|\hat g| < c_k \mid \mu)$ are **exactly invariant**: every quantity in
  the standardised limit $(c_k - \mu)/\sigma$ scales by $\alpha$ top and bottom.
- The density terms pick up a Jacobian $\alpha^{-1}$ per reporting observation — a constant that
  does not depend on $\mu$, so it cannot shift the maximiser.

Therefore $\hat\mu \to \alpha\hat\mu$ and the likelihood is otherwise unchanged: **no amount of
coordinate data distinguishes $\alpha$ from $1$.** The absolute effect-size scale is not weakly
identified, not hard to estimate — it is *structurally unavailable*. Seven attempts to recover it
failed for this reason, not for want of a better optimiser.

What survives is everything that doesn't depend on $\alpha$:

| quantity | scales with $\alpha$? | usable without images? |
|---|---|---|
| spatial pattern of $\hat\mu$ | no (common factor) | **yes** |
| $z = \hat g / \mathrm{se}$ | no (cancels) | **yes** |
| p-values, all corrected maps | no | **yes** |
| prevalence $\pi$ | no (a probability) | **yes** |
| absolute magnitude of $\hat g$ | **yes** | no |

This is the formal content of "treat $g$ as a relative map." It is a statement about
identification, in Manski's sense, and no confidence interval closes the gap: a CI shrinks with
$K$, identification bounds do not.

### II.4 Why inference must come from a spatial null

$\hat\mu / \mathrm{se}$ is **not** null-referenced. The peaks being pooled were selected for
being large, so under a global null the pooled effect *at a reported focus* is large by
construction. Referring it to a normal distribution rejected 41% of voxels at a nominal 5%.

The fix is A8: refer the statistic to the distribution obtained when coordinates fall at random.
Measured over 30 global-null simulations, the uncorrected rate is 0.039 against a nominal 0.050
and 0.0011 against 0.0010, with corrected rejection in 0.03–0.07 of whole simulations.

### II.5 The residual bias, and why it is predictable

A9 is false in a specific, knowable way: a reported peak is not a draw conditioned on exceeding
$c_k$, it is a **local maximum** of a smooth field that exceeded $c_k$. Maxima are biased upward
beyond simple truncation — the classic winner's curse, and in a spatial setting a *selection over
a neighbourhood* rather than over a point.

Measured on the 21 NIDM pain studies: the pooled truth at a reported peak is $0.22$–$0.27$ where
the reporting study says $1.22$, i.e. a peak overstates the local effect about **fivefold**, and
the pooled map ends up about **twice** the image-based reference. The `peak_bias="per-study"`
correction removes the component that varies between studies — because $\rho_k$ is a function of
$(c_k, N_k)$, both observable — and leaves the common factor, which is exactly the $\alpha$ of
§II.3.

### II.6 Where the usefulness comes from

1. **It uses information the convergence methods discard.** ALE and MKDA use only *where* foci
   are. A study reporting $t = 8$ at $N = 60$ and one reporting $t = 3.2$ at $N = 12$ contribute
   identically to ALE. CBES reads the magnitudes, and on the pain data reaches
   $\rho = 0.84$ against the image-based truth where ALE gets $0.19$ — not a criticism of ALE,
   which is measuring convergence and is being scored on a quantity it never claimed to estimate.
2. **It models the process that generated the data.** Thresholded reporting is the defining
   feature of the coordinate literature, and a method that ignores it is fitting the wrong
   likelihood.
3. **It separates two questions that convergence conflates.** $\pi(v)$ answers *how many studies
   show an effect here*; $\mu(v)$ answers *how big is it where it exists*. A focus reported by
   30 of 30 studies at modest magnitude and one reported by 8 of 30 at large magnitude are
   scientifically different and give the same ALE value.
4. **It degrades gracefully.** With no images: a relative magnitude map with valid inference.
   With five images out of twenty-one: magnitude ratio $0.98$ against the full-image answer.

---

## Part III. How these sit against the field's existing assumptions

### Assumptions CBES shares with standard neuroimaging

**Smoothness.** Every RFT-based analysis assumes the statistic field is smooth Gaussian. A2 is
the same assumption doing a different job: here it justifies *weighting*, not thresholding.
Notably, an attempt to use RFT properly — modelling silence as "the maximum over a neighbourhood
failed to clear $u$" — was correct in theory and worse in practice, because it needs a smoothness
that coordinates do not carry and that varies *more within a brain* (1.99×) than *between
studies* (1.6×). Pointwise censoring is biased in a way that can be characterised; the regional
version is biased in a way that cannot.

**Spatial normalisation.** All CBMA assumes coordinates are comparable across studies after
normalisation to a common template. The kernel absorbs part of this error; none of the methods
model it explicitly.

**Independence of studies (A7).** Shared by every meta-analysis, routinely violated by
overlapping author groups and shared datasets, and no CBMA method addresses it.

### Assumptions CBES adds

**A4 (height thresholding)** is the substantive addition, and it is the most fragile. Reporting
one local maximum per cluster satisfies it — on the pain studies the imposed threshold is
recovered to within $0.02\,z$. A **cluster-extent threshold does not**: requiring $k \ge 10, 20,
50$ voxels lifts the smallest reported statistic by $0.19$, $0.32$ and $0.57\,z$, and the
order-statistic correction recovers only about half of that. Since extent thresholding is at
least as common as height thresholding in the literature, this is the assumption most likely to
be violated in a real corpus.

**A6 (full coverage)** is assumed silently by every method that treats non-reporting as
informative. An ROI study that never examined a voxel is not evidence of a null effect there.

### Against the alternatives

| method | what it estimates | selection modelled? | magnitudes used? | key assumption |
|---|---|---|---|---|
| **ALE** | convergence of foci | no | no | foci are a point process; null is spatial randomness |
| **MKDA** | convergence, study-weighted | no | no | as ALE, with study as the unit |
| **CBMR** | foci intensity | no | no | foci counts are Poisson/NB with a spline intensity |
| **SDM-PSI** | effect size | via **imputation** within bounds | yes | bounds are correct; imputation is MAR-like |
| **IBMA** | effect size | not needed | yes | full images available |
| **CBES** | effect size + prevalence | via **censored likelihood** | yes | A1–A9 above |

The sharpest contrast is with **SDM-PSI**, the only other coordinate method estimating effect
size. SDM computes per-study *identification bounds* from the peaks and then multiply-imputes
within them; CBES never imputes an image, and lets non-reporting enter only through the
probability of silence. Head to head on the same 21 studies and the same derived coordinates,
scored against the images both were trying to recover:

| | SDM-PSI | CBES |
|---|---|---|
| localization where signal is ($r$, top quartile) | 0.68–0.73 | **0.80–0.81** |
| magnitude ratio (1.0 is perfect) | 0.42–0.47 | 2.04–2.06 |
| with 5 of 21 studies as images | 0.76 / 0.53 | **0.91 / 0.98** |

Neither is calibrated on coordinates alone, and they fail in **opposite directions** —
imputation toward the interior of bounds that straddle zero pulls SDM down; uncorrected
peak-height inflation pushes CBES up. Worth knowing: SDM's own published validation reports
*relative* mean square error against its own earlier kernel, never an absolute calibration slope,
so a systematic factor-of-two shared by both arms would not show up in it.

---

## Part IV. The same thing, without equations

### Mountains in the fog

Imagine a mountain range you can never see directly. Hundreds of hikers have walked it, and each
one sent back a postcard.

Here's the catch: **every hiker walked in fog**, and each one's fog sat at a different height.
A hiker whose fog was low could see a lot of the range. A hiker whose fog was high saw almost
nothing. And the postcards only ever mention the summits that poked *above* that hiker's fog —
"there's a peak at these coordinates, and it looked about this tall." Everything below the fog
line went unmentioned, whether it was a respectable hill or dead-flat prairie.

Your job: reconstruct the mountain range. That is coordinate-based meta-analysis.

**The older methods count postcards.** If forty hikers all mention a summit near the same spot,
there's probably really a mountain there. That's ALE and MKDA, and it works well — but notice
what it throws away. A postcard saying "enormous peak, unmistakable" counts exactly the same as
one saying "small bump, barely visible." Pins on a map, no heights.

**CBES reads the heights too.** That's the whole idea. But reading them honestly means dealing
with four problems.

**Problem one: every hiker used a different fog.** A hiker whose fog was very high only mentions
enormous summits, so their postcards look impressive — not because their mountains were bigger,
but because their fog hid everything else. Before comparing postcards you have to work out how
high each hiker's fog was. Usually they don't say, so CBES infers it: *the smallest summit a
hiker bothered to mention is just above their fog line*, and if they mentioned only a few, the
smallest of those is probably somewhat above it — an amount you can calculate.

**Problem two: a summit is not the mountain.** When a hiker reports "this peak was 1,200 metres,"
that is the single highest point they could see. The ridge around it is lower. So if you average
the reported summit heights, you get a range that's far too tall — roughly **five times** too
tall, measured. This is the same reason the winner of a race is usually a bit lucky as well as a
bit fast: you're looking at the *maximum*, and maxima are flattering. CBES discounts each
hiker's summits by an amount that depends on their fog height and how good their eyesight was.

**Problem three — the interesting one: silence is ambiguous.** A hiker says nothing about the
eastern valley. Two completely different things could be true. There might be a real hill there,
just under their fog. Or the east might be **genuinely flat**.

Older thinking forced one answer: silence means "small." But then if thirty hikers say nothing
about a place where five hikers saw a clear peak, the maths concludes the mountain must be tiny —
even *negative*, a hole in the ground, which is nonsense.

CBES allows both explanations at once and lets the postcards decide the mix. It estimates two
separate things for every spot:

- **how often** there's a mountain there at all (call it prevalence), and
- **how tall** it is where it does exist.

This matters more than it sounds. "Thirty out of thirty hikers found a modest hill" and "eight
out of thirty found a huge peak" are completely different findings about the world — and the
old pin-counting methods give them the same score.

**Problem four: the fundamental one.** Here's the part no cleverness fixes.

Suppose every hiker's altimeter was miscalibrated by the same factor — all reading twice the true
height. Could you tell from the postcards? **No.** The tall mountains would still be tall
relative to the short ones, the pattern of where mountains are would be perfect, but every
absolute number would be doubled and *nothing in the data would reveal it*. That's not a
limitation of the method — it's a property of the information. The postcards pin down the
*shape* of the range and not its *scale*.

So CBES is honest about it: the map of where the mountains are, which are bigger than which, how
confident you can be — all trustworthy. The claim "this peak is exactly 1,200 metres" is not,
unless someone hands you one properly calibrated altimeter.

And here's the good news, measured: **one** hiker with a real altimeter — one study that shares
its full data instead of just a postcard — fixes most of it. Five fixes it almost completely.
Ten is barely better than five.

**Last piece: how do you know a cluster of postcards isn't a coincidence?** Take all the reported
summits, scatter them to random spots on the map, and rebuild the range from the scattered
version. Do that a thousand times. If the real pile-up is bigger than almost every random one,
it's a real mountain. That's the permutation test, and CBES borrows it unchanged from the methods
that came before it.

### In one paragraph

CBES treats each reported coordinate as a measurement of effect *size*, not just a pin in a map;
corrects each study for how strictly it was thresholded and how many subjects it had; models
silence as genuinely ambiguous between "no effect" and "effect too small to report," which lets
it estimate how *often* an effect occurs separately from how *big* it is; and gets its p-values
by shuffling coordinates rather than trusting a formula that the selection has already broken.
It recovers *where* effects are and *how they rank* far better than counting foci does. It cannot
recover their absolute scale from coordinates alone — that is a theorem, not a bug — and one or
two studies contributing real images is enough to supply it.
