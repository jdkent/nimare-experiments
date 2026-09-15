"""Would a profile-likelihood interval be calibrated where the Wald one is not? Scalar first.

The shipped width comes from the observed information with the prevalence profiled out by a
Schur complement, and it is conservative wherever the indicator is thin: `se/sd` runs 1.28 to
2.03 on the field bed, and the delta-method width on `pi*mu` is worse still because it inverts
the whole near-singular 2x2. A profile likelihood inverts nothing -- the interval is the set of
`mu` whose profile log-likelihood sits within 1.921 of the maximum -- and needs no degrees of
freedom, so it would settle the `dof` question at the same time.

Whether it is actually better is a property of this likelihood, not of profiling in general, and
it answers in one dimension in milliseconds. So it is tested here before anything is built into
the estimator: the same censored zero-inflated mixture, one voxel, known truth, both intervals
from the same fit.

  * **Wald** -- 1.96 * sqrt(1 / I_profiled), the shipped recipe, where I_profiled is the Schur
    complement of the 2x2 observed information.
  * **Profile** -- {mu : 2 (l_p(mu_hat) - l_p(mu)) <= 3.841}, with `pi` re-maximised at each mu.

Reported per configuration: coverage of the true `mu`, mean half-width, and the ratio of the
reported half-width to 1.96 times the actual spread of `mu_hat` across replications, which is
the same `se/sd` diagnostic the full-scale beds use, so the two are comparable.

Calibration check: the `no indicator, pi fixed at 1` row reduces to a textbook Gaussian mean,
where both intervals must read close to 1.00. If they do not, the harness is wrong.
"""
import numpy as np
from scipy import optimize, stats

RNG = np.random.default_rng(20250915)
CRIT = stats.chi2.isf(0.05, 1)  # 3.841


def silent_prob(mu, sigma, cut):
    """P(|g| < cut) for g ~ N(mu, sigma^2)."""
    return stats.norm.cdf((cut - mu) / sigma) - stats.norm.cdf((-cut - mu) / sigma)


def log_likelihood(mu, pi, values, value_sd, null_sd, signs, ind_sd, cut):
    """Log-likelihood of the zero-inflated censored mixture at one voxel.

    ``values`` are the image studies' observed g. ``signs`` are the indicator studies:
    +1 silent, -1 reported. The null component is an effect of exactly zero, so its densities
    and probabilities are the same expressions at mu = 0.
    """
    pi = np.clip(pi, 1e-9, 1.0 - 1e-9)
    total = 0.0
    if len(values):
        active = stats.norm.pdf(values, mu, value_sd)
        null = stats.norm.pdf(values, 0.0, null_sd)
        total += np.sum(np.log(pi * active + (1.0 - pi) * null + 1e-300))
    if len(signs):
        p_active = silent_prob(mu, ind_sd, cut)
        p_null = silent_prob(0.0, ind_sd, cut)
        active = np.where(signs > 0, p_active, 1.0 - p_active)
        null = np.where(signs > 0, p_null, 1.0 - p_null)
        total += np.sum(np.log(pi * active + (1.0 - pi) * null + 1e-300))
    return total


def profile(mu, fixed_pi, *args):
    """max over pi of the log-likelihood at this mu, or the likelihood at a fixed pi."""
    if fixed_pi is not None:
        return log_likelihood(mu, fixed_pi, *args)
    grid = np.linspace(0.01, 0.99, 25)
    best = max(grid, key=lambda p: log_likelihood(mu, p, *args))
    out = optimize.minimize_scalar(
        lambda p: -log_likelihood(mu, p, *args),
        bounds=(max(1e-6, best - 0.08), min(1.0 - 1e-6, best + 0.08)),
        method="bounded",
    )
    return -out.fun


def fit(args, fixed_pi):
    """Maximise the profile in mu by a coarse grid then a bounded refine."""
    grid = np.linspace(-1.5, 1.5, 61)
    scores = [profile(m, fixed_pi, *args) for m in grid]
    seed = grid[int(np.argmax(scores))]
    out = optimize.minimize_scalar(
        lambda m: -profile(m, fixed_pi, *args),
        bounds=(seed - 0.06, seed + 0.06),
        method="bounded",
    )
    return float(out.x), -float(out.fun)


def wald_half_width(mu_hat, pi_hat, args, fixed_pi):
    """1.96 / sqrt(Schur complement of the observed information), by finite differences.

    Finite differences rather than the closed form deliberately: this bed is asking whether the
    *quantity* the estimator reports is calibrated, so a second implementation of it that agrees
    with the algebra is worth more here than a copy of the algebra.
    """
    h = 1e-4
    def ll(m, p):
        return log_likelihood(m, p, *args)
    i_mm = -(ll(mu_hat + h, pi_hat) - 2 * ll(mu_hat, pi_hat) + ll(mu_hat - h, pi_hat)) / h**2
    if fixed_pi is not None:
        return 1.96 / np.sqrt(i_mm) if i_mm > 0 else np.inf
    i_pp = -(ll(mu_hat, pi_hat + h) - 2 * ll(mu_hat, pi_hat) + ll(mu_hat, pi_hat - h)) / h**2
    i_mp = -(
        ll(mu_hat + h, pi_hat + h) - ll(mu_hat + h, pi_hat - h)
        - ll(mu_hat - h, pi_hat + h) + ll(mu_hat - h, pi_hat - h)
    ) / (4 * h**2)
    profiled = i_mm - (i_mp**2 / i_pp if i_pp > 0 else 0.0)
    return 1.96 / np.sqrt(profiled) if profiled > 0 else np.inf


def profile_half_width(mu_hat, peak, args, fixed_pi):
    """Half the width of {mu : 2 (l_p(mu_hat) - l_p(mu)) <= 3.841}, averaged over the two sides."""
    def deficit(m):
        return 2.0 * (peak - profile(m, fixed_pi, *args)) - CRIT

    edges = []
    for direction in (-1.0, 1.0):
        step, far = 0.05, None
        for _ in range(80):
            candidate = mu_hat + direction * step
            if deficit(candidate) > 0:
                far = candidate
                break
            step *= 1.3
        if far is None:
            return np.inf
        edges.append(optimize.brentq(deficit, mu_hat, far, xtol=1e-4))
    return float((edges[1] - edges[0]) / 2.0)


CONFIGS = [
    ("no indicator, pi fixed at 1", 2, 0, 1.0, True),
    ("no indicator, pi free", 2, 0, 1.0, False),
    ("2 images, 18 indicators", 2, 18, 1.0, False),
    ("2 images, 18 indicators, pi 0.5", 2, 18, 0.5, False),
]
N_REPS = 300
TRUE_MU, N_SUB, CUT = 0.5, 30.0, 0.65

if __name__ == "__main__":
    sd_value = np.sqrt(1.0 / N_SUB + TRUE_MU**2 / (2 * N_SUB))
    sd_null = np.sqrt(1.0 / N_SUB)
    print(f"true mu {TRUE_MU}, n {int(N_SUB)}, cutoff {CUT} g, {N_REPS} replications")
    print("  se/sd is the reported half-width over 1.96 times the actual spread of mu_hat\n")
    head = (f"{'configuration':>32} {'cov W':>6} {'cov P':>6} {'half W':>7} {'half P':>7} "
            f"{'se/sd W':>8} {'se/sd P':>8}")
    print(head)
    for label, n_images, n_ind, true_pi, fix in CONFIGS:
        hats, wald, prof = [], [], []
        cov_w, cov_p = [], []
        for _ in range(N_REPS):
            live = RNG.random(max(n_images, 1)) < true_pi
            values = np.where(
                live[:n_images],
                RNG.normal(TRUE_MU, sd_value, n_images),
                RNG.normal(0.0, sd_null, n_images),
            )
            if n_ind:
                live_ind = RNG.random(n_ind) < true_pi
                p = np.where(live_ind, silent_prob(TRUE_MU, sd_null, CUT),
                             silent_prob(0.0, sd_null, CUT))
                signs = np.where(RNG.random(n_ind) < p, 1.0, -1.0)
            else:
                signs = np.array([])
            args = (values, sd_value, sd_null, signs, sd_null, CUT)
            fixed = 1.0 if fix else None
            mu_hat, peak = fit(args, fixed)
            pi_hat = 1.0
            if not fix:
                grid = np.linspace(0.01, 0.99, 99)
                pi_hat = float(grid[int(np.argmax(
                    [log_likelihood(mu_hat, p, *args) for p in grid]))])
            hw = wald_half_width(mu_hat, pi_hat, args, fixed)
            hp = profile_half_width(mu_hat, peak, args, fixed)
            hats.append(mu_hat)
            wald.append(hw)
            prof.append(hp)
            cov_w.append(abs(mu_hat - TRUE_MU) <= hw)
            cov_p.append(abs(mu_hat - TRUE_MU) <= hp)
        spread = 1.96 * float(np.std(hats, ddof=1))
        mw, mp = float(np.median(wald)), float(np.median(prof))
        print(f"{label:>32} {np.mean(cov_w):6.2f} {np.mean(cov_p):6.2f} {mw:7.3f} {mp:7.3f} "
              f"{mw / spread:8.2f} {mp / spread:8.2f}")
