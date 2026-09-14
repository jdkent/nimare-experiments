"""Item 5: does the censored likelihood's curvature SE cover, and can it be corrected?

Under the default selection model, `se` is 1/sqrt(-curvature) of the censored mixture likelihood
at the fitted mu -- the observed-information standard error. HKSJ fixed the *other* SE, the
pooled inverse-variance one, lifting coverage from 86.4% to 95.9% at three studies; but the
zero-inflated path discards that SE entirely, so the default output has no small-sample
correction and has never been checked.

Observed information is asymptotic in the number of contributing studies, and here that number
is a handful. Three things are expected to hurt, and are separated below:

  * the asymptotic normal reference, where a t on few effective studies is the honest one;
  * mu and pi being estimated jointly, so the curvature in mu alone ignores the cost of not
    knowing pi -- the profile curvature overstates the information;
  * the plateau at one reporting study, where curvature is near zero and the SE explodes.

Simulated at the level the estimator actually solves: one voxel, known mu and pi, studies that
report or stay silent according to the model, fitted through the estimator's own _fit_chunk.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from scipy import stats
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import null_effect_variance

CUTOFF = 0.55
TRIALS = 3000


def one_trial(rng, n_studies, true_mu, true_pi):
    """Draw one voxel's worth of data from the model, then fit it."""
    sizes = rng.integers(20, 41, n_studies).astype(float)
    null_var = null_effect_variance(sizes, design="one-sample")[:, None]
    has_effect = rng.random(n_studies) < true_pi
    weights = np.zeros((n_studies, 1))
    g_obs = np.zeros((n_studies, 1))
    var_obs = np.ones((n_studies, 1))
    covered = np.zeros((n_studies, 1), dtype=bool)

    n_reporting = 0
    for study in range(n_studies):
        sd = np.sqrt(float(null_var[study, 0]))
        drawn = rng.normal(true_mu if has_effect[study] else 0.0, sd)
        if abs(drawn) >= CUTOFF:  # cleared the reporting threshold
            weights[study, 0] = 1.0
            g_obs[study, 0] = drawn
            var_obs[study, 0] = float(null_var[study, 0])
            covered[study, 0] = True
            n_reporting += 1
    if n_reporting == 0:
        return None

    estimator = CBES(fwhm=8.0, null_method="none", max_iter=200)
    mu, pi, se = estimator._fit_chunk(
        weights=weights,
        g_obs=g_obs,
        var_obs=var_obs,
        covered=covered,
        tau2=np.zeros(1),
        null_var=null_var,
        cutoffs=np.full((n_studies, 1), CUTOFF),
        start=np.array([float(np.mean(g_obs[covered]))]),
    )
    return float(mu[0]), float(se[0]), n_reporting


print(f"{TRIALS} trials per cell; cutoff |g| = {CUTOFF}; nominal 95% intervals\n")
print(f"{'studies':>8s} {'true mu':>8s} {'true pi':>8s} {'reporting':>10s} {'bias':>7s} "
      f"{'normal':>7s} {'t(k-1)':>8s} {'t(k-2)':>8s} {'k>1':>9s}")
for n_studies in (10, 20, 40):
    for true_mu, true_pi in ((0.8, 1.0), (0.8, 0.5), (0.4, 0.5)):
        rng = np.random.default_rng(7)
        # Each reference gets its own denominator. A t on k - 1 degrees of freedom is undefined
        # when one study reported, so those trials are excluded from its rate rather than
        # counted as misses -- doing the latter read as 41.9% coverage in the sparse cells and
        # was an artefact of the accounting, not of the interval.
        hits = {"z": 0, "t1": 0, "t2": 0}
        counted = {"z": 0, "t1": 0, "t2": 0}
        biases, reporting = [], []
        for _ in range(TRIALS):
            got = one_trial(rng, n_studies, true_mu, true_pi)
            if got is None:
                continue
            mu, se, k = got
            if not np.isfinite(se) or se <= 0:
                continue
            biases.append(mu - true_mu)
            reporting.append(k)
            error = abs(mu - true_mu)
            counted["z"] += 1
            hits["z"] += error <= 1.959964 * se
            if k > 1:
                counted["t1"] += 1
                hits["t1"] += error <= stats.t.ppf(0.975, k - 1) * se
            if k > 2:
                counted["t2"] += 1
                hits["t2"] += error <= stats.t.ppf(0.975, k - 2) * se
        if not counted["z"]:
            continue

        def rate(key):
            return f"{100 * hits[key] / counted[key]:.1f}%" if counted[key] else "n/a"

        print(f"{n_studies:8d} {true_mu:8.2f} {true_pi:8.2f} {np.mean(reporting):10.1f} "
              f"{np.mean(biases):+7.3f} {rate('z'):>7s} {rate('t1'):>8s} {rate('t2'):>8s} "
              f"{100 * counted['t1'] / counted['z']:8.0f}%", flush=True)
