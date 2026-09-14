"""Does the observed-information SE cover, where the EM curvature SE did not?

The review's finding 2: `se` was the curvature of the EM's Q function with responsibilities
held fixed, so it kept r*h per observation and dropped r(1-r)s^2 -- the part attributable to not
knowing which mixture component an observation came from -- and ignored the jointly estimated
prevalence entirely. Dropping a positive term from a negative curvature overstates information,
so the reported error was too small. The review measured 62.5-89.8% coverage of nominal-95%
intervals, not improving with more studies.

This runs the review's own design against the fixed code: the estimator's idealized censored
mixture, latent active state Bernoulli(pi), known sampling SD, tau2 = 0, reporting iff the
absolute observed effect clears a cutoff, unit kernel weights. Fixing heterogeneity and
variances at their true values isolates mixture uncertainty, which is the point.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from scipy import stats
from nimare.meta.cbma import CBES

TRIALS = 2000
SD = 0.2


def one_trial(rng, k, true_mu, true_pi, cutoff):
    """One voxel's worth of data from the model the estimator assumes, then fit it."""
    var = SD**2
    null_var = np.full((k, 1), var)
    active = rng.random(k) < true_pi
    weights = np.zeros((k, 1))
    g_obs = np.zeros((k, 1))
    var_obs = np.ones((k, 1))
    covered = np.zeros((k, 1), dtype=bool)
    n_reporting = 0
    for study in range(k):
        drawn = rng.normal(true_mu if active[study] else 0.0, SD)
        if abs(drawn) >= cutoff:
            weights[study, 0] = 1.0
            g_obs[study, 0] = drawn
            var_obs[study, 0] = var
            covered[study, 0] = True
            n_reporting += 1
    if n_reporting < 2:
        return None

    estimator = CBES(fwhm=8.0, null_method="none", max_iter=400)
    mu, pi, se = estimator._fit_chunk(
        weights=weights,
        g_obs=g_obs,
        var_obs=var_obs,
        covered=covered,
        tau2=np.zeros(1),
        null_var=null_var,
        cutoffs=np.full((k, 1), cutoff),
        start=np.array([float(np.mean(g_obs[covered]))]),
    )
    if not np.isfinite(se[0]) or se[0] <= 0:
        return None
    interior = _PREVALENCE_MARGIN < pi[0] < 1.0 - _PREVALENCE_MARGIN
    if not interior:
        return None
    return float(mu[0]), float(se[0]), n_reporting


_PREVALENCE_MARGIN = 1e-3

print(f"{TRIALS} trials per cell; sampling SD {SD}; nominal 95% intervals")
print("regular fits only: more than one report, interior prevalence, positive information\n")
print(f"{'K':>4s} {'mu':>5s} {'pi':>5s} {'c':>6s} {'fits':>6s} {'bias':>7s} "
      f"{'coverage (t on k-1)':>21s}")
for k, true_mu, true_pi, cutoff in (
    (40, 0.4, 0.5, 0.350),
    (40, 0.6, 0.5, 0.550),
    (100, 0.4, 0.5, 0.350),
    (100, 0.3, 0.7, 0.350),
    (40, 0.8, 0.5, 0.550),
    (40, 0.8, 0.5, 0.658),
    (100, 0.6, 0.5, 0.658),
):
    rng = np.random.default_rng(11)
    hits = used = 0
    biases = []
    for _ in range(TRIALS):
        got = one_trial(rng, k, true_mu, true_pi, cutoff)
        if got is None:
            continue
        mu, se, reporting = got
        used += 1
        biases.append(mu - true_mu)
        critical = stats.t.ppf(0.975, max(reporting - 1, 1))
        hits += abs(mu - true_mu) <= critical * se
    if not used:
        continue
    print(f"{k:4d} {true_mu:5.2f} {true_pi:5.2f} {cutoff:6.3f} {used:6d} "
          f"{np.mean(biases):+7.3f} {100 * hits / used:20.1f}%", flush=True)
