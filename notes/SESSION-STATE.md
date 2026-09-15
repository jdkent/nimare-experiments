# Where CBES stands — estimands and uncertainty

Written for a check-in. Everything here is measured in this repo; the PR carries only the
algorithmic conclusions. Detail and method are in `notes/cbes-open-program.md`,
`notes/cbes-first-principles.md` and `notes/cbes-success-criteria.md`.

## The three things worth knowing

**1. The interval on `g` has no working operating point, and the failure is bias, not width.**
100 replications against a known truth of 0.800 with prevalence 1 (so no estimand mismatch), in
a reporting regime calibrated to real practice (every study reports, 3-4 clusters each):

| studies | images | bias | coverage | interval half-width / truth |
| --- | --- | --- | --- | --- |
| 12 | 0 | +0.506 | 0.10 | 0.44 |
| 12 | 0, `peak_bias='per-study'` | +0.514 | 0.16 | 0.47 |
| 12 | 2 | +0.248 | 0.67 | 0.40 |
| 12 | 6 | +0.079 | 0.85 | 0.23 |
| 12 | 12 | +0.010 | 0.96 | 0.14 |
| 24 | 0 | +0.489 | **0.00** | 0.28 |
| 24 | 6 | +0.176 | 0.56 | 0.23 |

The `se` is honest about the estimator's own variability everywhere (reported width / replication
sd runs 1.0 to 1.7, i.e. slightly conservative). Coverage nonetheless collapses because the point
estimate is biased, and **more studies make coverage worse**, not better — a fixed bias with a
shrinking interval. The 12-of-12 row is the reference, not a result: donor coordinates are
dropped, so with every study imaged the coordinate table is empty and that fit is an IBMA. It
does confirm that the pooling and the standard error code are correct.

**2. The bias is one formula, and it depends on the image *share*, not the number of images.**

    bias(f) = 0.506 * (1 - f) / ((1 - f) + 5.4 f)        f = fraction of studies with images

Fitted with one free parameter; max residual 0.010 across biases spanning 0.496. The same six
donors give +0.079 at 6-of-12 and +0.176 at 6-of-24. So the coordinate channel is never being
*corrected*, only *diluted* — which is exactly why `peak_bias='per-study'` moves nothing
(+0.506 → +0.514). Practical reading: "at least 2 images" is necessary and nowhere near
sufficient. Two images among twenty studies is f = 0.10 and leaves about +0.37. Getting under 10%
of the effect needs half the collection imaged.

**3. Why, in one sentence: the information in a coordinate table is in the counts and locations,
not the heights.** Reported peak height moves 0.055 z per unit of true g at a 3.29 cut while the
cluster count moves about ten times as much; one coordinate explains 5–9% of the variance in the
truth at its own location. Every approach here that used counts and locations worked (the
point-process fit recovered the scale to 1–5%); every approach that used reported magnitudes
failed in the direction selection predicts. CBES puts its whole estimating equation in the weak
channel and uses the strong one only to decide which voxels to touch.

## The theoretical result that ties it together

A coordinate observation is an **occupancy record with imperfect detection**: study *k* reports
near voxel *v* with probability `pi_v * D(mu_v, n_k, u_k)`. In that model, occupancy and detection
are separately identified only through repeat visits or through covariates that shift detection
without shifting occupancy. Studies are the repeat visits; sample size and reporting threshold are
the detection covariates.

So `prevalence` and `g` are separately identified **only by the spread of sample sizes and
reporting thresholds across the studies covering a voxel**. With a homogeneous roster only the
product is identified. That single mechanism explains four things previously recorded as separate
puzzles: why `g_marginal` is the best-behaved magnitude map (it is close to the identified
combination, not a lucky cancellation), why `prevalence` has a floor near 0.2 and goes
non-monotone at weak effects, why the split is unstable while the product is stable, and why no
amount of work on the heights fixes it.

It also suggests a diagnostic the estimator can compute from the collection in hand — the spread
of sample sizes among studies contributing at each voxel — and report the product while refusing
the split where that spread is small. The test that would confirm or kill it is running.

## Corrections to things I previously told you

- **`prevalence` ordering "survives".** Tested as stated — ordering across voxels *within one
  map*, which is the claim the docstring makes. Strong effect: mean Spearman +0.862 but the
  four-site ranking is exactly right in only **50%** of maps. Weak effect: +0.586 with sd 0.341
  and exactly right in **12%**. And the compression changes sign across the range (a true 0.25
  reads 0.36, a true 1.00 reads 0.93), so spacing is not interpretable. "Reliable on average over
  many maps, for strong effects" is the defensible claim.
- **"Coordinates never help alongside images."** Measured at one effect size. At a weak effect
  with a single image, 40 coordinate studies take the correlation from 0.804 to 0.893 — a real
  gain. At two or more images, at either effect size, coordinates add nothing and eventually
  subtract. The ≥2-image scoping lands precisely on the regime where the joint model has no
  advantage, which is worth saying out loud.
- **I claimed a bug where there is a design inconsistency.** I thought the all-donor path
  smuggled a refused configuration past the guard, and I patched it; that was wrong and is
  reverted. `_permute_image_values` randomises images over their own voxels, which tests the same
  hypothesis, so that path is legitimate — it is the *images-only refusal's stated reason* ("no
  foci to permute") that is now false. Needs your decision, not a patch (task #47).
- **I ran a 100-replication coverage study in the wrong regime** before calibrating it: 18% of
  studies reported anything, so a 12-study analysis was really a 2-study one and a third of
  replications were unfittable, with survivors selected for high signal. Discarded and rerun
  after calibrating. The rule is tightened in the notes: a calibration check means a printed
  measurement of the reporting regime, refused unless it matches published practice.

## What needs your decision

- **#47** Images-only refusal vs the all-donor path: allow both (evidence favours this) or refuse
  both. Either way the current error message's reasoning is wrong.
- **#48** Whether a simulation-derived bias number belongs in a warning, or only the mechanism
  and the share.
- **#36** Where `g` comes from when a collection has images.
- Whether `g`'s target should be restated as relative-only, with `se` withdrawn or documented as
  sampling-only, and `g_marginal` and `prevalence` presented as the deliverables.

## Still running or queued

Bias decomposition by stage; whether the coordinates-only bias is one constant or distorts the
pattern; `g_marginal`'s cancellation against reporting density; fitted `prevalence` against a
naive count of studies with a nearby focus; and the identification-by-power-spread test above.
