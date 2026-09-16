r"""The scalar interval-censored random-effects reference, before any of it is implemented.

This is step 1 of the design document's build order (its section 6.3): "a Gaussian or
noncentral-t scalar censored model with real bounds, jointly treating heterogeneity. This
establishes transformations, derivatives, inference, and reference agreement." Everything the
implementation will compute is derived here first, because a derivative coded from a guess and
checked against a simulation is a derivative fitted to a bed.

Model. A study's latent effect is :math:`\theta_i \sim N(m, \tau^2)` and its estimate is
:math:`Y_i \mid \theta_i \sim N(\theta_i, s_i^2)`, so marginally
:math:`Y_i \sim N(m, s_i^2 + \tau^2)`. An observation is an *event* about :math:`Y_i`, not a
value: an interval :math:`(L_i, U_i)`, possibly one-sided, possibly degenerate.

What is proved, and why each one is load-bearing for the code:

  * the two-level model collapses to one Gaussian, and the interval probability is that
    Gaussian's increment (claims 1-2) -- this is the whole likelihood;
  * a degenerate interval reproduces the ordinary random-effects log-density up to a constant
    free of the parameters, so the same optimiser recovers standard random-effects
    meta-analysis when every study is observed (claims 3-4). Section 11 of the document lists
    that agreement as a validation requirement; here it is an identity, not a test;
  * the one-sided limits (claims 5-6);
  * the exact score in :math:`(m, \tau^2)` (claims 7-8), which is what the implementation codes
    and what section 11 asks to be checked against numerical derivatives;
  * **the reporting partition** (claims 9-10). The document's section 4.1 says the shipped
    estimator drops an annulus and then uses complementary probabilities for the two events it
    keeps. Claim 9 gives the error exactly: using :math:`1 - p_r` in place of :math:`p_s`
    overstates the silence probability by the dropped event's own probability, and claim 10
    gives the conditional probabilities that would be correct among retained events. So the
    defect is quantified rather than asserted, and the size of it is an expression the
    implementation can evaluate;
  * **the mixture-variance identity** (claims 11-12). Section 4.3 states
    :math:`\operatorname{Var}(\theta) = \pi\tau_a^2 + \pi(1-\pi)\mu^2`; claim 12 turns that into
    the exact amount by which an image-only between-study variance overstates the active
    component's variance, which is the inconsistency the document objects to;
  * **what identifies** :math:`(m, \tau^2)` **from intervals alone** (claim 13). Two studies
    with the *same standardised bounds* still identify both parameters as long as their sampling
    variances differ, and the cross-determinant is exactly proportional to
    :math:`\sigma_1 - \sigma_2`. This is the :math:`\tau^2` analogue of the established result
    that sample-size spread identifies the magnitude (#73), reached independently, and it is not
    a restatement of the rank-one theorem in ``mixture_identification.py`` -- that file proves
    the determinant formula, this one evaluates a particular gradient in it.

**Stated calibration target, before any simulation is written.** The document reports its own
4,000-replication scalar experiment (its section 8.1, true mean .4, within-study SD .2,
between-study SD .15, directional threshold .6, half of exceedances retained, Monte Carlo SE
0.34 percentage points near 95%):

    1 + 20   image-only RMSE .256   correct model .129   ignoring retention .157
                                    coverage 93.95%      coverage 77.60%
    1 + 500  image-only RMSE .251   correct model .023   ignoring retention .114
                                    coverage 94.93%      coverage 0.00%
    8 + 100  image-only RMSE .090   correct model .047   ignoring retention .103
                                    coverage 94.75%      coverage 23.45%

The implementation is calibrated against that table, not against a bed of my own construction.
A disagreement is a defect in one of the two and must be resolved before anything else here is
believed.

**Regimes this cannot speak to.** Every claim is scalar and Gaussian with a *known* bound. It
says nothing about where the bounds come from, which is the document's central objection and
remains unaddressed by any of it; nothing about the noncentral-t observation model of its
section 2; and nothing about spatial dependence.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

m, y, theta, h = sp.symbols("m y theta h", real=True)
L, U = sp.symbols("L U", real=True)
s2, tau2 = sp.symbols("s^2 tau^2", positive=True)
sigma = sp.Symbol("sigma", positive=True)

proof = Proof(
    "scalar_censored_reference",
    "The scalar interval-censored random-effects model, derived before it is implemented",
    __doc__,
)


def normal_pdf(value, mean, sd):
    return sp.exp(-((value - mean) ** 2) / (2 * sd**2)) / (sd * sp.sqrt(2 * sp.pi))


def normal_cdf(value, mean, sd):
    return (1 + sp.erf((value - mean) / (sd * sp.sqrt(2)))) / 2


# ------------------------------------------------- the two levels collapse to one Gaussian

total_sd = sp.sqrt(s2 + tau2)
convolution = sp.integrate(
    normal_pdf(y, theta, sp.sqrt(s2)) * normal_pdf(theta, m, sp.sqrt(tau2)),
    (theta, -sp.oo, sp.oo),
)
proof.claim(
    "the-two-level-model-marginalises-to-one-gaussian",
    sp.simplify(convolution - normal_pdf(y, m, total_sd)),
    r"\int \phi(y;\theta,s)\,\phi(\theta;m,\tau)\,d\theta = \phi(y;m,\sqrt{s^2+\tau^2})",
)

interval_probability = normal_cdf(U, m, total_sd) - normal_cdf(L, m, total_sd)
proof.define(r"P(L < Y < U)", interval_probability)
proof.claim(
    "and-the-interval-probability-differentiates-back-to-that-density",
    sp.simplify(sp.diff(interval_probability, U) - normal_pdf(U, m, total_sd)),
    r"\partial_U\left[\Phi\!\left(\tfrac{U-m}{\sigma}\right)"
    r"-\Phi\!\left(\tfrac{L-m}{\sigma}\right)\right] = \phi(U;m,\sigma)",
)

# ---------------------------------- a degenerate interval is the ordinary random-effects term

# Centre a width-h interval on the observed value and expand. The h^1 coefficient is the
# density and the h^2 coefficient vanishes, so log P = log h + log phi + O(h^2).
narrow = normal_cdf(y + h / 2, m, sigma) - normal_cdf(y - h / 2, m, sigma)
series = sp.series(narrow, h, 0, 3).removeO()
proof.claim(
    "a-narrow-interval-carries-the-density-in-its-leading-term",
    sp.simplify(series.coeff(h, 1) - normal_pdf(y, m, sigma)),
    r"P(|Y-y|<h/2) = h\,\phi(y;m,\sigma) + O(h^3)",
)
proof.claim(
    "with-no-second-order-term-so-the-width-factors-out-of-the-log-likelihood",
    sp.simplify(series.coeff(h, 2)),
    r"\log P = \log h + \log\phi(y;m,\sigma) + O(h^2),"
    r"\quad \partial_{(m,\tau^2)}\log h = 0",
)

# --------------------------------------------------------------- the one-sided limits

proof.claim(
    "sending-the-lower-bound-to-minus-infinity-leaves-a-left-censored-term",
    sp.simplify(
        sp.limit(interval_probability, L, -sp.oo) - normal_cdf(U, m, total_sd)
    ),
    r"P(Y<U) = \Phi\!\left(\tfrac{U-m}{\sigma}\right)",
)
proof.claim(
    "and-the-upper-bound-to-plus-infinity-a-right-censored-one",
    sp.simplify(
        sp.limit(interval_probability, U, sp.oo) - (1 - normal_cdf(L, m, total_sd))
    ),
    r"P(Y>L) = 1 - \Phi\!\left(\tfrac{L-m}{\sigma}\right)",
)

# --------------------------------------------------------------------- the exact score

a = (L - m) / total_sd
b = (U - m) / total_sd
phi = sp.Function("phi")
increment = normal_cdf(U, m, total_sd) - normal_cdf(L, m, total_sd)
loglik = sp.log(increment)

standard_pdf = lambda t: sp.exp(-t**2 / 2) / sp.sqrt(2 * sp.pi)  # noqa: E731
score_m = -(standard_pdf(b) - standard_pdf(a)) / (total_sd * increment)
proof.claim(
    "the-score-in-the-mean-is-the-density-difference-over-the-interval-mass",
    sp.simplify(sp.diff(loglik, m) - score_m),
    r"\partial_m \ell = -\frac{\varphi(b)-\varphi(a)}{\sigma\,[\Phi(b)-\Phi(a)]}",
)
score_tau2 = -(b * standard_pdf(b) - a * standard_pdf(a)) / (2 * (s2 + tau2) * increment)
proof.claim(
    "and-the-score-in-the-heterogeneity-weights-each-bound-by-its-own-standardised-value",
    sp.simplify(sp.diff(loglik, tau2) - score_tau2),
    r"\partial_{\tau^2} \ell = -\frac{b\varphi(b)-a\varphi(a)}"
    r"{2\sigma^2\,[\Phi(b)-\Phi(a)]}",
)

# ------------------------------------------------------- the reporting partition (section 4.1)

p_r, p_a, p_s = sp.symbols("p_r p_a p_s", positive=True)
partition = sp.Eq(p_r + p_a + p_s, 1)
proof.define(r"p_r + p_a + p_s = 1", partition)
silence = 1 - p_r - p_a
proof.define(r"p_s", silence)

# Substituting the complement of the report probability for the silence probability overstates
# it by exactly the probability of the event that was dropped.
proof.claim(
    "using-one-minus-the-report-probability-for-silence-overstates-it-by-the-dropped-event",
    sp.simplify((1 - p_r) - silence - p_a),
    r"(1 - p_r) - p_s = p_a",
)
# The probabilities that would be correct among retained events are the conditional ones, and
# they differ from the marginals by a factor the dropped event controls.
retained = p_r + silence
proof.claim(
    "and-the-correct-conditional-report-probability-exceeds-the-marginal-by-a-known-factor",
    sp.simplify(p_r / retained - p_r - p_r * p_a / retained),
    r"\frac{p_r}{p_r+p_s} - p_r = \frac{p_r\,p_a}{p_r+p_s}",
)

# ------------------------------------------------ the mixture-variance identity (section 4.3)

pi, mu, tau_a2 = sp.symbols("pi mu tau_a^2", positive=True)
first_moment = pi * mu
second_moment = pi * (mu**2 + tau_a2)
mixture_variance = second_moment - first_moment**2
proof.claim(
    "a-zero-inflated-effect-has-variance-pi-tau-a-squared-plus-pi-one-minus-pi-mu-squared",
    sp.simplify(mixture_variance - (pi * tau_a2 + pi * (1 - pi) * mu**2)),
    r"\operatorname{Var}(\theta) = \pi\tau_a^2 + \pi(1-\pi)\mu^2",
)
proof.claim(
    "so-dividing-an-image-only-between-study-variance-by-the-prevalence-overstates-the-active-one",
    sp.simplify(mixture_variance / pi - tau_a2 - (1 - pi) * mu**2),
    r"\frac{\operatorname{Var}(\theta)}{\pi} - \tau_a^2 = (1-\pi)\mu^2 \ge 0",
)

# ------------------------------- what identifies (m, tau^2) from intervals alone

# Two studies sharing standardised bounds (a, b) but differing in sampling variance. Write the
# shared bound-dependent factors as A and B; the gradients are then (A/sigma, B/sigma^2).
A, B = sp.symbols("A B", nonzero=True)
sigma_1, sigma_2 = sp.symbols("sigma_1 sigma_2", positive=True)
gradient_1 = sp.Matrix([A / sigma_1, B / sigma_1**2])
gradient_2 = sp.Matrix([A / sigma_2, B / sigma_2**2])
cross = gradient_1[0] * gradient_2[1] - gradient_1[1] * gradient_2[0]
proof.claim(
    "identical-standardised-bounds-still-identify-both-parameters-when-the-precisions-differ",
    sp.simplify(cross - A * B * (sigma_1 - sigma_2) / (sigma_1**2 * sigma_2**2)),
    r"a_{11}a_{22} - a_{12}a_{21} = \frac{AB(\sigma_1-\sigma_2)}{\sigma_1^2\sigma_2^2}",
)

if __name__ == "__main__":
    import numpy as np
    from scipy.special import ndtr
    from scipy.stats import norm

    proof.report()

    # Numeric check of the claims the implementation will code, against finite differences at
    # parameter values chosen away from every boundary.
    rng = np.random.default_rng(11)
    worst_m, worst_t = 0.0, 0.0
    for _ in range(200):
        m_v, tau2_v, s2_v = rng.normal(0.3, 0.4), rng.uniform(0.01, 0.4), rng.uniform(0.01, 0.4)
        lo, hi = np.sort(rng.normal(0.3, 1.0, size=2))
        sd = np.sqrt(s2_v + tau2_v)

        def ell(mean, between):
            scale = np.sqrt(s2_v + between)
            return np.log(ndtr((hi - mean) / scale) - ndtr((lo - mean) / scale))

        a_v, b_v = (lo - m_v) / sd, (hi - m_v) / sd
        mass = ndtr(b_v) - ndtr(a_v)
        analytic_m = -(norm.pdf(b_v) - norm.pdf(a_v)) / (sd * mass)
        analytic_t = -(b_v * norm.pdf(b_v) - a_v * norm.pdf(a_v)) / (2 * (s2_v + tau2_v) * mass)
        step = 1e-6
        numeric_m = (ell(m_v + step, tau2_v) - ell(m_v - step, tau2_v)) / (2 * step)
        numeric_t = (ell(m_v, tau2_v + step) - ell(m_v, tau2_v - step)) / (2 * step)
        worst_m = max(worst_m, abs(analytic_m - numeric_m) / max(abs(numeric_m), 1e-8))
        worst_t = max(worst_t, abs(analytic_t - numeric_t) / max(abs(numeric_t), 1e-8))
    if not (worst_m < 1e-5 and worst_t < 1e-5):
        raise AssertionError(
            f"the analytic scores disagree with finite differences: {worst_m:.2e}, {worst_t:.2e}"
        )
    print(f"  [ok] numeric: both scores match finite differences to "
          f"{max(worst_m, worst_t):.1e} relative over 200 random parameter draws")

    # Numeric check of the degenerate-interval limit: the censored log-likelihood minus log(h)
    # must approach the ordinary Gaussian log-density.
    m_v, sd_v, y_v = 0.35, 0.28, 0.51
    for width in (1e-2, 1e-4, 1e-6):
        mass = ndtr((y_v + width / 2 - m_v) / sd_v) - ndtr((y_v - width / 2 - m_v) / sd_v)
        gap = abs(np.log(mass) - np.log(width) - norm.logpdf(y_v, m_v, sd_v))
        # O(h^2) in theory, but the CDF difference cancels catastrophically at small h, so
        # the check cannot demand better than double precision allows.
        if gap > max(width**2, 1e-9):
            raise AssertionError(f"the degenerate limit is off by {gap:.2e} at width {width:g}")
    print("  [ok] numeric: a degenerate interval reproduces the Gaussian log-density to O(h^2)")

    print()
    print("Calibration target for the implementation, stated here before it exists: the")
    print("document's own 4,000-replication table (its section 8.1). 1+20 .256/.129/.157 with")
    print("coverage 93.95/77.60; 1+500 .251/.023/.114 with 94.93/0.00; 8+100 .090/.047/.103")
    print("with 94.75/23.45. Monte Carlo SE near 95% coverage is 0.34 percentage points.")
    print("A disagreement is a defect in the implementation or in that table, and has to be")
    print("resolved before anything built on this is believed.")
