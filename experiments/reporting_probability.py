"""Can the reporting probability be estimated and standardised to a reference design?

Section 19 proposes an estimand that avoids everything that has gone wrong with `g`:

    P_k(v) = probability that study k reports a focus within r of voxel v
           = pi_v * D(mu_v, n_k, u_k)

It is the expectation of an observed Bernoulli, so it needs no scale constant, no image donor and
no factoring of prevalence from magnitude. The question is whether it can be *standardised*: fit
the detection dependence on sample size and threshold, then report the probability a *reference*
study would find a focus there. Without that it is only a description of the collection in hand.

The truth is computable exactly here. For a study of size n reporting at threshold u, the chance
of a supra-threshold observation at a site of true magnitude mu is

    D = P(|Z| >= u)  with  Z ~ Normal(mu*sqrt(n), 1)      (the bed's generating process)

so `P = pi * D` needs no simulation to know.

Scored two ways, because a model can fit the collection and still not transport:
  calibration   over the collection's own designs -- does a predicted 0.4 happen 40% of the time
  transport     at reference designs inside and outside the collection's range of n, which is
                where the proposal's honesty claim lives: it should be right in range and it
                should be visibly wrong out of it, not quietly wrong.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats, optimize

SITES_MU = np.array([0.0, 0.2, 0.4, 0.6, 0.8])
SITES_PI = np.array([1.0, 1.0, 0.75, 0.75, 0.50])
N_SIMS = int(os.environ.get("NSIMS", 200))
N_STUDIES = 24
THRESHOLDS = (2.5758, 3.0902, 3.2905, 3.7190)
#: The collection's own sample sizes. Transport is then tested inside and outside this.
N_RANGE = (20, 60)
REFERENCE_N = (15, 30, 45, 60, 100, 200)
REFERENCE_U = 3.2905


def true_P(mu, pi, n, u):
    """Exact reporting probability under the bed's generating process."""
    lam = mu * np.sqrt(n)
    return pi * (stats.norm.sf(u - lam) + stats.norm.cdf(-u - lam))


def simulate(rng):
    """One collection: per study a design, and per site whether it reported."""
    designs = [(int(rng.integers(*N_RANGE)), float(rng.choice(THRESHOLDS)))
               for _ in range(N_STUDIES)]
    reports = np.zeros((N_STUDIES, SITES_MU.size), dtype=int)
    for k, (n, u) in enumerate(designs):
        has = rng.random(SITES_MU.size) < SITES_PI
        z = rng.normal(SITES_MU * np.sqrt(n), 1.0)
        reports[k] = has & (np.abs(z) >= u)
    return designs, reports


def fit_site(designs, reported):
    """Fit (pi, mu) at one site from the detection records alone -- no heights used.

    The likelihood is the occupancy one: each study is a Bernoulli with probability
    pi * D(mu, n_k, u_k). Both parameters enter through the *same* records, so this is exactly
    the separation section 18 says needs a spread of designs.
    """
    n = np.array([d[0] for d in designs], dtype=float)
    u = np.array([d[1] for d in designs], dtype=float)

    def negll(theta):
        pi = 1.0 / (1.0 + np.exp(-theta[0]))
        mu = theta[1]
        p = np.clip(true_P(mu, pi, n, u), 1e-9, 1 - 1e-9)
        return -float(np.sum(reported * np.log(p) + (1 - reported) * np.log1p(-p)))

    best = None
    for pi0 in (-1.0, 0.0, 1.5):
        for mu0 in (0.2, 0.5, 0.9):
            out = optimize.minimize(negll, [pi0, mu0], method="Nelder-Mead",
                                    options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 2000})
            if best is None or out.fun < best.fun:
                best = out
    pi = 1.0 / (1.0 + np.exp(-best.x[0]))
    return pi, float(best.x[1])


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    print(f"{N_STUDIES} studies per collection, n ~ U{N_RANGE}, {len(THRESHOLDS)} thresholds, "
          f"{N_SIMS} collections")
    print(f"sites: mu = {SITES_MU},  pi = {SITES_PI}\n")

    fitted_pi, fitted_mu = [], []
    for _ in range(N_SIMS):
        designs, reports = simulate(rng)
        row_pi, row_mu = [], []
        for s in range(SITES_MU.size):
            p, m = fit_site(designs, reports[:, s])
            row_pi.append(p); row_mu.append(m)
        fitted_pi.append(row_pi); fitted_mu.append(row_mu)
    fitted_pi = np.array(fitted_pi); fitted_mu = np.array(fitted_mu)

    print("--- the factors: means are useless here, so the distribution is shown ---")
    print("        true mu  " + " ".join(f"{v:7.2f}" for v in SITES_MU))
    print("      median mu  " + " ".join(f"{v:7.3f}" for v in np.median(fitted_mu, axis=0)))
    print("        mean mu  " + " ".join(f"{v:7.3f}" for v in fitted_mu.mean(0)))
    print("     90th pct mu " + " ".join(f"{v:7.3f}"
                                         for v in np.percentile(fitted_mu, 90, axis=0)))
    # A fit with mu*sqrt(n) far above the cut has detection pinned at 1, so mu is unidentified
    # from above: everything beyond that point predicts the same thing. That is the runaway.
    runaway = np.mean(fitted_mu * np.sqrt(np.mean(N_RANGE)) > THRESHOLDS[-1] + 3.0, axis=0)
    print("   fits with detection saturated  " + " ".join(f"{v:7.2f}" for v in runaway))
    print("        true pi  " + " ".join(f"{v:7.2f}" for v in SITES_PI))
    print("      median pi  " + " ".join(f"{v:7.3f}" for v in np.median(fitted_pi, axis=0)))
    print("        mean pi  " + " ".join(f"{v:7.3f}" for v in fitted_pi.mean(0)))
    # How flat is the detection curve at the truth? d D / d mu says how much the data can
    # distinguish a change in magnitude from a change in prevalence.
    nbar = float(np.mean(N_RANGE))
    grad = np.array([(true_P(m + 0.05, 1.0, nbar, REFERENCE_U)
                      - true_P(m - 0.05, 1.0, nbar, REFERENCE_U)) / 0.1 for m in SITES_MU])
    print("   dD/dmu at the truth (N=%d)   " % nbar
          + " ".join(f"{v:7.3f}" for v in grad))

    print("\n--- the estimand: reporting probability at a reference design ---")
    print(f"reference threshold u = {REFERENCE_U}; the collection spans n in {N_RANGE}")
    header = "  ".join(f"mu={m:.1f}" for m in SITES_MU)
    print(f"{'ref N':>7s}  {'in range':>8s}   {header}")
    for n0 in REFERENCE_N:
        truth = true_P(SITES_MU, SITES_PI, n0, REFERENCE_U)
        est = np.array([true_P(fitted_mu[:, s], fitted_pi[:, s], n0, REFERENCE_U).mean()
                        for s in range(SITES_MU.size)])
        inside = "yes" if N_RANGE[0] <= n0 <= N_RANGE[1] else "NO"
        cells = "  ".join(f"{e:.2f}/{t:.2f}" for e, t in zip(est, truth))
        print(f"{n0:7d}  {inside:>8s}   {cells}")
    print("\nCells are estimated/true. The proposal needs the in-range rows to match and the")
    print("out-of-range rows to degrade in a way that can be disclosed, not hidden.")
