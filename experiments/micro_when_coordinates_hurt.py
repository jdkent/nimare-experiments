"""The crossing where a biased channel stops helping, in closed form and then in simulation.

Two channels estimate the same mu. The images are unbiased with variance `1 / (k n)`, falling as
image studies accumulate. The coordinates carry a bias `b` from the assumptions they are read
through -- an assumed reporting threshold, an assumed coverage radius -- and that bias does not
fall as coordinate studies accumulate, only their variance `1 / (m I)` does.

Inverse-variance weighting puts `w = R / (1 + R)` on the coordinates, with
`R = m I / (k n)` the ratio of the two channels' information. Writing out the mean squared error
of the combination and asking when it beats `1 / (k n)`:

    (w b)^2 + (1 - w)^2 / (k n) + w^2 / (m I)  <  1 / (k n)

and every term in `R` cancels, leaving

    b^2  <  1 / (m I) + 1 / (k n).

So the combination helps exactly while the squared channel bias sits below the *sum* of the two
channels' variances -- which is the variance of their difference. That makes the crossing a
Hausman statistic:

    Q = (g_img - g_coord)^2 / (var_img + var_coord),      helps while Q < 1.

Both quantities in Q are things an estimator already reports, and neither depends on knowing the
truth, so this is testable on real data. What it predicts: the crossing image count is
`k* = 1 / (n (b^2 - 1 / (m I)))`, and it is infinite -- coordinates never hurt -- whenever the
coordinate channel's own variance already exceeds its squared bias.

This bed checks the algebra against simulation, with the bias arising the way it really does:
the studies report at a true threshold and the estimator reads them at a wrong one.
"""
import numpy as np
from scipy import optimize, stats

import os

N_PER_STUDY = 30
M_TABLES = 18
REPS = int(os.environ.get("REPS", 4000))

#: Multiplier on the reported standard error, to model the estimator's known SE inflation
#: (se/sd measured at 1.2 to 2.2). Q divides by a variance, so an se inflated by `f` deflates Q
#: by `f**2` and the `Q > 1` crossing lands late. The weights are a ratio of informations and so
#: are untouched by a uniform inflation, which is why only Q's denominator carries it here.
SE_INFLATION = float(os.environ.get("INFLATE", 1.0))


def report_rate(mu, cut, n):
    """Two-sided probability that a study of size n reports, at threshold `cut` on the g scale."""
    root = np.sqrt(n)
    return stats.norm.sf((cut - mu) * root) + stats.norm.cdf((-cut - mu) * root)


def coordinate_channel(mu, cut_true, cut_assumed, n, m):
    """Asymptotic bias and variance of the estimate read off m reporting indicators.

    The indicators are generated at `cut_true`; the estimator inverts the rate at `cut_assumed`,
    so it converges to the mu that would have produced the observed rate under the wrong cut.
    """
    p_true = report_rate(mu, cut_true, n)
    # The two-sided rate is symmetric in mu and so not monotone on a bracket spanning zero;
    # it increases in |mu|, so the inversion has to run on one side of it.
    lo, hi = 0.0, 6.0
    if report_rate(lo, cut_assumed, n) > p_true or report_rate(hi, cut_assumed, n) < p_true:
        return np.nan, np.nan
    mu_limit = optimize.brentq(
        lambda x: report_rate(x, cut_assumed, n) - p_true, lo, hi, xtol=1e-10)
    step = 1e-5
    slope = ((report_rate(mu_limit + step, cut_assumed, n)
              - report_rate(mu_limit - step, cut_assumed, n)) / (2 * step))
    if abs(slope) < 1e-9:
        return np.nan, np.nan
    variance = p_true * (1.0 - p_true) / (m * slope**2)
    return mu_limit - mu, variance


def simulate(mu, cut_true, cut_assumed, n, m, k, rng):
    """One realisation: k images pooled, m indicators inverted, then combined."""
    g_img = mu + rng.normal(0.0, 1.0 / np.sqrt(k * n))
    var_img = 1.0 / (k * n)
    reports = rng.random(m) < report_rate(mu, cut_true, n)
    rate = reports.mean()
    floor, ceiling = report_rate(0.0, cut_assumed, n), report_rate(6.0, cut_assumed, n)
    if not floor < rate < ceiling:
        return np.nan, np.nan, np.nan
    g_coord = optimize.brentq(
        lambda x: report_rate(x, cut_assumed, n) - rate, 0.0, 6.0, xtol=1e-8)
    step = 1e-5
    slope = ((report_rate(g_coord + step, cut_assumed, n)
              - report_rate(g_coord - step, cut_assumed, n)) / (2 * step))
    var_coord = max(rate * (1.0 - rate) / (m * slope**2), 1e-9)
    weight = (1.0 / var_coord) / (1.0 / var_coord + 1.0 / var_img)
    combined = weight * g_coord + (1.0 - weight) * g_img
    # The Hausman statistic the algebra says marks the crossing, in the form an estimator can
    # compute: it needs only the two channels' disagreement and the variance of that difference.
    q = ((g_coord - g_img) ** 2
         / ((var_img + var_coord) * SE_INFLATION**2))
    return combined, g_img, q


if __name__ == "__main__":
    mu, cut_true = 0.5, 0.55
    print(f"mu {mu}, n {N_PER_STUDY}, m {M_TABLES} tables, true cut {cut_true} g, "
          f"{REPS} replications, se inflation {SE_INFLATION:.2f}\n")
    print(f"{'assumed cut':>12} {'bias b':>8} {'var_c':>8} {'b^2':>8} "
          f"{'predicted k*':>13}")
    scenarios = []
    for cut_assumed in (0.55, 0.60, 0.70, 0.85, 1.10):
        b, var_c = coordinate_channel(mu, cut_true, cut_assumed, N_PER_STUDY, M_TABLES)
        excess = b**2 - var_c
        kstar = (1.0 / (N_PER_STUDY * excess)) if excess > 0 else np.inf
        scenarios.append((cut_assumed, b, var_c, kstar))
        print(f"{cut_assumed:12.2f} {b:+8.3f} {var_c:8.4f} {b**2:8.4f} "
              f"{kstar:13.1f}")

    print(f"\nmeasured mse of the combination minus mse of the image pool, by image count")
    print(f"{'assumed cut':>12} " + " ".join(f"{f'k={k}':>9}" for k in (1, 2, 4, 8, 16, 32))
          + f" {'crossing':>9} {'predicted':>10}")
    rng = np.random.default_rng(0)
    for cut_assumed, b, var_c, kstar in scenarios:
        diffs, qs, crossing, q_crossing = [], [], None, None
        for k in (1, 2, 4, 8, 16, 32):
            pair = np.array([simulate(mu, cut_true, cut_assumed, N_PER_STUDY, M_TABLES, k, rng)
                             for _ in range(REPS)])
            ok = np.isfinite(pair).all(axis=1)
            mse_combined = np.mean((pair[ok, 0] - mu) ** 2)
            mse_images = np.mean((pair[ok, 1] - mu) ** 2)
            diffs.append(mse_combined - mse_images)
            qs.append(float(np.median(pair[ok, 2])))
            if crossing is None and diffs[-1] > 0:
                crossing = k
            if q_crossing is None and qs[-1] > 1.0:
                q_crossing = k
        print(f"{cut_assumed:12.2f} " + " ".join(f"{d:+9.5f}" for d in diffs)
              + f" {str(crossing):>9} {kstar:10.1f}")
        print(f"{'median Q':>12} " + " ".join(f"{v:9.3f}" for v in qs)
              + f" {str(q_crossing):>9} {'(Q > 1)':>10}")

    print("\nA positive difference means the coordinates are costing accuracy. The crossing "
          "column is the\nfirst image count at which that happens; predicted is "
          "1 / (n (b^2 - var_c)) from the algebra.")
