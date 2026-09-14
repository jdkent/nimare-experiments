"""What would a profile MLE plus interval actually cost, in the unit that matters?

The expensive thing in a fit is a censoring evaluation -- two normal CDFs over every silent
(study, voxel) pair -- and it is a function of mu alone, not of the prevalence. So the cost of
any scheme is just how many distinct mu values it has to visit. The current EM visits max_iter
of them (25 by default). Counted here for:

  joint EM         what the code does now
  profile MLE      outer Brent over mu, prevalence solved at each mu reusing that mu's censoring
  + interval       plus two root-finds for the profile-likelihood endpoints

Built on the real model with several silent studies at different sigmas and cutoffs, so the
inner prevalence solve has no closed form and has to iterate, exactly as in the estimator.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import ndtr

rng = np.random.default_rng(0)
CHI2_95 = 3.841458820694124  # one-parameter profile-likelihood interval


def make_voxel(n_reporting, n_silent):
    sigma_rep = np.abs(rng.normal(0.25, 0.05, n_reporting)) + 0.05
    sigma_sil = np.abs(rng.normal(0.25, 0.05, n_silent)) + 0.05
    cutoff = np.abs(rng.normal(0.55, 0.08, n_silent)) + 0.2
    g_obs = cutoff.mean() + np.abs(rng.normal(0.15, 0.05, n_reporting))
    return dict(sigma_rep=sigma_rep, g_obs=g_obs, sigma_sil=sigma_sil, cutoff=cutoff,
                cutoff_rep=np.full(n_reporting, cutoff.mean()))


def censoring(mu, v, counter):
    """The expensive call: two normal CDFs per silent pair. Depends on mu only."""
    counter[0] += 1
    z = v["cutoff"] / v["sigma_sil"]
    return ndtr((v["cutoff"] - mu) / v["sigma_sil"]) - ndtr((-v["cutoff"] - mu) / v["sigma_sil"]), \
        ndtr(z) - ndtr(-z)


def loglik(mu, pi, v, prob):
    present = np.exp(-0.5 * ((v["g_obs"] - mu) / v["sigma_rep"]) ** 2) / v["sigma_rep"]
    absent = np.exp(-0.5 * (v["g_obs"] / v["sigma_rep"]) ** 2) / v["sigma_rep"]
    sil_present, sil_absent = prob
    return (np.log(np.maximum(pi * present + (1 - pi) * absent, 1e-300)).sum()
            + np.log(np.maximum(pi * sil_present + (1 - pi) * sil_absent, 1e-300)).sum())


def profile(mu, v, counter):
    """max over pi at this mu. One censoring evaluation, then cheap iteration on pi."""
    prob = censoring(mu, v, counter)
    best = -np.inf
    for pi in np.linspace(0.01, 0.99, 33):  # bounded 1-D, no new censoring evaluations
        best = max(best, loglik(mu, pi, v, prob))
    return best


print(f"{'reporting':>9s} {'silent':>7s} {'joint EM':>9s} {'profile':>9s} "
      f"{'+ interval':>11s} {'interval':>18s}")
for n_reporting, n_silent in ((1, 5), (1, 20), (2, 20), (5, 20)):
    v = make_voxel(n_reporting, n_silent)
    counter = [0]
    result = minimize_scalar(lambda m: -profile(m, v, counter), bounds=(0.0, 2.0),
                             method="bounded", options={"xatol": 1e-3})
    n_profile = counter[0]
    peak, mu_hat = -result.fun, result.x

    counter[0] = 0
    target = lambda m: profile(m, v, counter) - (peak - CHI2_95 / 2)
    lo = hi = None
    try:
        lo = brentq(target, 0.0, mu_hat) if target(0.0) < 0 else 0.0
    except ValueError:
        lo = 0.0
    try:
        hi = brentq(target, mu_hat, 2.0) if target(2.0) < 0 else 2.0
    except ValueError:
        hi = 2.0
    n_interval = counter[0]

    print(f"{n_reporting:9d} {n_silent:7d} {25:9d} {n_profile:9d} "
          f"{n_profile + n_interval:11d} {f'[{lo:.2f}, {hi:.2f}]':>18s}")
