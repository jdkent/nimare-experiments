"""Is the absolute effect size identified by reporting *rates*, without images?

The peak heights carry no magnitude information, but whether a study reports at all does:
P(report) = Phi(g*sqrt(N) - u) is a probit in g whose slope is pinned by the sample size. So
across studies of differing N, the pattern of who reported should identify g on an absolute
scale -- the sample size acting as an exclusion restriction in the Heckman sense.

This checks the claim directly, with no estimator involved: simulate studies with a known g,
throw the reported heights away entirely, and fit g from the reporting indicators alone.
If the profile likelihood has a sharp peak at the truth, the information is there.
"""
import warnings; warnings.simplefilter("ignore")
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr

rng = np.random.default_rng(0)

def simulate(true_g, n_studies, n_range, thresholds, tau=0.1, seed=0):
    """Who reports? A study reports if its observed peak z clears its own threshold."""
    r = np.random.default_rng(seed)
    sizes = r.integers(*n_range, size=n_studies).astype(float)
    u = r.choice(thresholds, size=n_studies)
    study_g = r.normal(true_g, tau, size=n_studies)
    observed_g = r.normal(study_g, 1.0 / np.sqrt(sizes))
    reported = np.abs(observed_g) * np.sqrt(sizes) > u
    return sizes, u, reported

def loglik(g, sizes, u, reported, tau=0.1):
    """Log-likelihood of the reporting pattern alone, at effect size g."""
    sd = np.sqrt(1.0 / sizes + tau**2)
    cut = u / np.sqrt(sizes)                       # threshold on the g scale
    p_silent = np.clip(ndtr((cut - g) / sd) - ndtr((-cut - g) / sd), 1e-12, 1 - 1e-12)
    return np.sum(np.where(reported, np.log1p(-p_silent), np.log(p_silent)))

def mle(sizes, u, reported, lo=0.0, hi=2.5):
    grid = np.linspace(lo, hi, 501)
    ll = np.array([loglik(g, sizes, u, reported) for g in grid])
    return grid[int(np.argmax(ll))], grid, ll

print("Recovering g from reporting indicators only -- the heights are never used.\n")
print(f"{'true g':>7s} {'studies':>8s} {'N range':>10s} {'reported':>9s} "
      f"{'g-hat':>7s} {'95% CI':>16s}")
for true_g in (0.2, 0.35, 0.5, 0.8):
    for n_studies, n_range in ((21, (9, 33)), (40, (10, 60)), (100, (10, 60))):
        estimates = []
        for seed in range(40):
            sizes, u, rep = simulate(
                true_g, n_studies, n_range, [2.3263, 3.0902, 3.2905, 4.2649], seed=seed
            )
            if rep.sum() == 0 or rep.sum() == len(rep):
                continue
            g_hat, _, _ = mle(sizes, u, rep)
            estimates.append(g_hat)
        estimates = np.array(estimates)
        if not estimates.size:
            continue
        lo, hi = np.percentile(estimates, [2.5, 97.5])
        frac = rep.mean()
        print(f"{true_g:7.2f} {n_studies:8d} {str(n_range):>10s} {frac:9.2f} "
              f"{estimates.mean():7.3f} {f'[{lo:.2f}, {hi:.2f}]':>16s}")
    print()
