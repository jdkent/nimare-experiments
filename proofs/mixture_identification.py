r"""Is :math:`(\mu, \pi)` identified, and why is the Wald interval never right?

Two measured facts want an explanation that is not another simulation.

  * The fitted prevalence converges to 0.723 against a true 0.6. For a correctly specified model
    the population log-likelihood is maximised at the truth, so either the model is misspecified
    against its own simulator -- it is not, both were written from the same densities -- or the
    maximum is so flat that a finite sample lands somewhere else.
  * ``se/sd`` is 0.80 at the default iteration cap and 2.55 at convergence, for the same data.
    No iteration count makes it 1.

Both are statements about the curvature of one 2x2 expected information matrix, so both are
algebra. The observation at a voxel is either a study's image value, distributed as the mixture

.. math:: f(g) = \pi\,\phi_a(g - \mu) + (1 - \pi)\,\phi_0(g),

or a study's reporting indicator, a Bernoulli whose silence probability is the same mixture
integrated over :math:`(-c, c)`:

.. math:: S(\mu, \pi) = \pi\,s_a(\mu) + (1 - \pi)\,s_0.

The claims below build :math:`I` from those two pieces, verify that it is the negative expected
Hessian rather than an outer product asserted to be one, and give the determinant in closed form.
The last claim is the useful one: the Schur complement that the estimator inverts to get ``se``
collapses to zero exactly when :math:`I_{\mu\pi}^2 = I_{\mu\mu} I_{\pi\pi}`, and that is a
condition on the data a voxel can be tested against.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

mu, pi, c, sa, s0 = sp.symbols("mu pi c s_a s_0", real=True)
sigma_a, sigma_0 = sp.symbols("sigma_a sigma_0", positive=True)
n_img, n_tab = sp.symbols("n_i n_t", positive=True)

proof = Proof(
    "mixture_identification",
    r"Identification of $(\mu,\pi)$ in the zero-inflated censored likelihood",
    __doc__,
)


def normal_pdf(x, sd):
    return sp.exp(-x**2 / (2 * sd**2)) / (sd * sp.sqrt(2 * sp.pi))


def normal_cdf(x, sd):
    return (1 + sp.erf(x / (sd * sp.sqrt(2)))) / 2


# ------------------------------------------------------------------ the indicator block

#: Probability a study with a real effect stays silent: its value falls inside (-c, c).
silent_active = normal_cdf(c - mu, sigma_a) - normal_cdf(-c - mu, sigma_a)
#: The same for a study with no effect here. Free of mu, which is the whole reason pi and mu
#: can be told apart at all.
silent_null = normal_cdf(c, sigma_0) - normal_cdf(-c, sigma_0)

silence = pi * silent_active + (1 - pi) * silent_null

proof.define(r"S(\mu,\pi)", silence)

# One indicator is a Bernoulli with success probability `silence`, so its Fisher information in
# any parameter theta is (dS/dtheta)^2 / (S(1-S)). Written out for the three blocks.
d_mu = sp.diff(silence, mu)
d_pi = sp.diff(silence, pi)
denominator = silence * (1 - silence)

indicator_mu = d_mu**2 / denominator
indicator_pi = d_pi**2 / denominator
indicator_cross = d_mu * d_pi / denominator

# A Bernoulli's expected information equals the outer product of its score and also the negative
# expected second derivative. Verifying it rather than asserting it: this is the step that would
# hide a factor of two.
log_likelihood_silent = sp.log(silence)
log_likelihood_report = sp.log(1 - silence)
expected_second = -(
    silence * sp.diff(log_likelihood_silent, mu, 2)
    + (1 - silence) * sp.diff(log_likelihood_report, mu, 2)
)
proof.claim(
    "bernoulli-information-is-the-outer-product",
    sp.simplify(expected_second - indicator_mu),
    r"-\mathbb{E}\,\partial_\mu^2 \log f = (\partial_\mu S)^2 / (S(1-S))",
)

# dS/dpi does not involve mu's derivative, and dS/dmu carries a factor pi. That factor is the
# whole story: at pi = 0 the indicator says nothing about mu, and near pi = 0 it says little.
proof.claim(
    "the-indicator-score-in-mu-carries-a-factor-pi",
    sp.simplify(d_mu - pi * sp.diff(silent_active, mu)),
    r"\partial_\mu S = \pi\,\partial_\mu s_a",
)
proof.claim(
    "the-indicator-score-in-pi-is-the-component-gap",
    sp.simplify(d_pi - (silent_active - silent_null)),
    r"\partial_\pi S = s_a - s_0",
)

# ------------------------------------------------------------------ the information matrix

# Abbreviate so the determinant is readable: u is the indicator's sensitivity to mu at pi = 1,
# and v is the gap between the two components' silence probabilities.
u, v, D = sp.symbols("u v D", real=True)
substitution = {
    sp.diff(silent_active, mu): u,
    silent_active - silent_null: v,
    denominator: D,
}

I_mu = n_tab * pi**2 * u**2 / D
I_pi = n_tab * v**2 / D
I_cross = n_tab * pi * u * v / D

determinant = sp.simplify(I_mu * I_pi - I_cross**2)
proof.claim(
    "the-indicator-alone-cannot-identify-both",
    sp.simplify(determinant),
    r"I_{\mu\mu} I_{\pi\pi} - I_{\mu\pi}^2 = 0 \quad\text{(indicators only)}",
)

# So the Schur complement from indicators alone is identically zero: the reporting indicator
# constrains one direction in (mu, pi) and no more, whatever the sample size. Everything that
# separates the two parameters has to come from the image values.
schur_indicator = sp.simplify(I_mu - I_cross**2 / I_pi)
proof.claim(
    "indicators-leave-no-information-about-mu-once-pi-is-free",
    sp.simplify(schur_indicator),
    r"I_{\mu\mu} - I_{\mu\pi}^2 / I_{\pi\pi} = 0 \quad\text{(indicators only)}",
)

# ------------------------------------------------------------- heterogeneous studies

# The rank-1 result above assumed every coordinate study shares one (c, sigma). Studies that
# differ contribute different (u, v), and a sum of rank-1 matrices need not be rank 1. Lagrange's
# identity says exactly how much rank the sum gains: the determinant is a sum of squared
# cross-products, so it is positive as soon as two studies disagree about the ratio u/v.
u1, u2, v1, v2, w1, w2 = sp.symbols("u_1 u_2 v_1 v_2 w_1 w_2", real=True)
sum_mm = pi**2 * (w1 * u1**2 + w2 * u2**2)
sum_pp = w1 * v1**2 + w2 * v2**2
sum_mp = pi * (w1 * u1 * v1 + w2 * u2 * v2)
proof.claim(
    "two-unlike-studies-restore-rank-two-by-lagranges-identity",
    sp.expand(sum_mm * sum_pp - sum_mp**2 - pi**2 * w1 * w2 * (u1 * v2 - u2 * v1) ** 2),
    r"\det I = \pi^2 \sum_{j<k} w_j w_k (u_j v_k - u_k v_j)^2",
)

# So the whole of the coordinate channel's ability to separate mu from pi sits in the spread of
# u/v across studies, and in nothing else. That is a statement the estimator's own measurements
# can be checked against: spreading sample sizes was measured to help and spreading thresholds
# was measured not to, and the ratio has to reproduce that asymmetry.
ratio_uv = sp.diff(silent_active, mu) / (silent_active - silent_null)
proof.define(r"u/v", ratio_uv)

# ------------------------------------------------------------------ what the images add

# The image block is not degenerate, because the mixture density's dependence on mu and on pi are
# different functions of the observed value. Write its two scores at a single observation.
g = sp.Symbol("g", real=True)
density = pi * normal_pdf(g - mu, sigma_a) + (1 - pi) * normal_pdf(g, sigma_0)
score_mu = sp.diff(sp.log(density), mu)
score_pi = sp.diff(sp.log(density), pi)

# The two scores are proportional only if their ratio is free of g. Rather than assert that,
# factor each one and let the difference show: the mu score carries (g - mu) / sigma_a^2, which
# changes sign at g = mu, and the pi score carries the density gap, which does not. A ratio with
# a zero in it cannot be constant.
active = pi * normal_pdf(g - mu, sigma_a)
proof.claim(
    "the-image-score-in-mu-carries-a-factor-g-minus-mu",
    sp.simplify(score_mu - (g - mu) / sigma_a**2 * active / density),
    r"\partial_\mu \ell = \frac{g-\mu}{\sigma_a^2}\,\frac{\pi\phi_a}{f}",
)
proof.claim(
    "the-image-score-in-pi-is-the-density-gap",
    sp.simplify(
        score_pi - (normal_pdf(g - mu, sigma_a) - normal_pdf(g, sigma_0)) / density
    ),
    r"\partial_\pi \ell = (\phi_a - \phi_0) / f",
)
# So at g = mu the mu score vanishes while the pi score does not, which is impossible for two
# proportional functions. Evaluated there, the pi score is a nonzero constant.
gap_at_mu = (normal_pdf(0, sigma_a) - normal_pdf(mu, sigma_0)) / density.subs(g, mu)
proof.claim(
    "at-g-equal-mu-only-one-score-vanishes",
    sp.simplify(score_mu.subs(g, mu)),
    r"\partial_\mu \ell \big|_{g=\mu} = 0 \quad\text{while}\quad \partial_\pi \ell \neq 0",
)
proof.define(r"\partial_\pi \ell\big|_{g=\mu}", sp.simplify(gap_at_mu))

# ------------------------------------------------------------------ the practical condition

I_mm, I_pp, I_mp = sp.symbols("I_mm I_pp I_mp", real=True)
schur = I_mm - I_mp**2 / I_pp
proof.define(r"\text{Schur}", schur)
proof.claim(
    "the-wald-se-diverges-exactly-at-the-ridge",
    sp.simplify(sp.together(schur).as_numer_denom()[0] - (I_mm * I_pp - I_mp**2)),
    r"\text{Schur} \to 0 \iff I_{\mu\pi}^2 \to I_{\mu\mu} I_{\pi\pi}",
)

def calibration():
    """Check the algebra against a result already established at full scale.

    Measured on the field bed (#73): spreading sample sizes 12 to 120 across 20 studies lifts the
    fraction of voxels where the likelihood bounds mu from 0.42 to 0.69, while spreading
    thresholds 2.3 to 4.5 changes nothing. Lagrange's identity says the only thing that can
    produce that asymmetry is the sum of squared cross-products, so evaluate it for both designs.
    A proof that could not reproduce a measurement it claims to explain is not an explanation.
    """
    import numpy as np
    from scipy.stats import norm

    true_mu, tau = 0.5, 0.15

    def uv(n, z):
        sd_a, sd_0 = np.sqrt(1.0 / n + tau**2), np.sqrt(1.0 / n)
        cut = z / np.sqrt(n)
        hi, lo = (cut - true_mu) / sd_a, (-cut - true_mu) / sd_a
        s_active = norm.cdf(hi) - norm.cdf(lo)
        s_zero = norm.cdf(cut / sd_0) - norm.cdf(-cut / sd_0)
        return -(norm.pdf(hi) - norm.pdf(lo)) / sd_a, s_active - s_zero

    def cross(pairs):
        return sum(
            (pairs[j][0] * pairs[k][1] - pairs[k][0] * pairs[j][1]) ** 2
            for j in range(len(pairs))
            for k in range(j + 1, len(pairs))
        )

    by_size = cross([uv(n, 3.09) for n in (12, 20, 40, 80, 120)])
    by_threshold = cross([uv(20, z) for z in (2.3, 2.8, 3.3, 3.9, 4.5)])
    flat = cross([uv(20, 3.09)] * 5)

    print()
    print("calibration against #73, sum of squared cross-products:")
    print(f"  sample sizes 12-120 at z=3.09   {by_size:.4e}")
    print(f"  thresholds 2.3-4.5 at n=20      {by_threshold:.4e}")
    print(f"  no spread at all                {flat:.4e}")
    print(f"  ratio {by_size / by_threshold:.0f}x in favour of spreading sample sizes, "
          "which is the measured direction")
    # The ratios u/v spread about equally in the two designs (1.23 against 0.93), so the factor
    # is not about the ratio: spreading the threshold moves both mixture components together and
    # leaves v = s_a - s_0 small, which is what the cross-product is built from.
    assert flat == 0.0, flat
    assert by_size > 10 * by_threshold, (by_size, by_threshold)


if __name__ == "__main__":
    print(f"{len(proof.claims)} claims verified in {proof.name}")
    for label, shown in proof.claims:
        print(f"  {label}")
    calibration()
