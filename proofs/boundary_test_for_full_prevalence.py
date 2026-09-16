r"""Can a collection be told that its prevalence is below 1? The test, and when it has no power.

The estimator's worst failure mode has no answer at present. Where every study carries the effect,
:math:`\pi < 1` is pure error and the correction makes ``g`` worse than doing nothing -- measured
on held-out HCP subjects, an images-only pool recovered 0.85 of the true magnitude where ``g``
recovered 0.63. Where studies genuinely differ, the same machinery helps. Nothing currently tells
a user which regime a collection is in, and ``prevalence`` itself cannot: it is fitted below 1
either way.

That is a hypothesis test, :math:`H_0: \pi = 1` against :math:`H_1: \pi < 1`, and the reason it
has not been written down is that the usual recipe does not apply. :math:`\pi = 1` sits on the
*boundary* of the parameter space, so the likelihood-ratio statistic is not
:math:`\chi^2_1`. Under the boundary conditions of Self and Liang (1987,
doi:10.1080/01621459.1987.10478472), following Chernoff (1954,
doi:10.1214/aoms/1177728725), it is the half-and-half mixture

.. math:: 2\log\Lambda \;\sim\; \tfrac12\chi^2_0 + \tfrac12\chi^2_1,

because half the time the unconstrained maximum falls outside the parameter space and the
constrained and unconstrained fits coincide. The practical consequence is a factor of two: the
nominal 5% critical value is the 90th percentile of :math:`\chi^2_1`, 2.706, not the 95th, 3.841.
Using 3.841 makes the test conservative by half; using the normal-theory 1.96 on a Wald statistic
makes it wrong in the other direction, since the Wald statistic is not even well defined when the
maximum is at the boundary.

The claims below assemble the statistic from quantities the estimator already computes, and then
answer the question that decides whether it is worth building: where does its power come from?
The rank-1 theorem in ``mixture_identification.py`` already says the coordinate channel cannot
separate :math:`\mu` from :math:`\pi`, so the answer cannot be "more tables".
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

mu, pi, eps = sp.symbols("mu pi epsilon", real=True)
n_img, n_tab = sp.symbols("n_i n_t", positive=True)

proof = Proof(
    "boundary_test_for_full_prevalence",
    r"Testing $\pi = 1$ against $\pi < 1$ on the boundary",
    __doc__,
)

# ------------------------------------------------------------- the score at the boundary

# Write the mixture density as a function of pi for one observation, with f1 the active
# component's density and f0 the null component's, whatever kind of observation it is: an image
# value contributes normal densities, an indicator contributes silence probabilities. The
# algebra below does not care which, which is the point -- the test is the same statistic either
# way.
f1, f0 = sp.symbols("f_1 f_0", positive=True)
density = pi * f1 + (1 - pi) * f0
log_density = sp.log(density)

# The score in pi is the standardized difference between the two components, and at pi = 1 it
# reduces to something computable without fitting anything: 1 - f0/f1.
score_pi = sp.diff(log_density, pi)
proof.claim(
    "the-score-in-pi-is-the-component-likelihood-ratio",
    sp.simplify(score_pi - (f1 - f0) / density),
    r"\partial_\pi \ell = (f_1 - f_0)/f",
)
proof.claim(
    "and-at-the-boundary-it-needs-no-fit",
    sp.simplify(score_pi.subs(pi, 1) - (1 - f0 / f1)),
    r"\partial_\pi \ell \big|_{\pi=1} = 1 - f_0/f_1",
)

# The alternative is one-sided: pi can only move *down* from 1. So the test rejects only when the
# summed score is negative, which is what makes the null distribution a mixture rather than a
# chi-square. Half the time the score points out of the parameter space and the statistic is 0.
proof.define(r"S", sp.Symbol("S", real=True))

# ------------------------------------------------------- the boundary null distribution

# The mixture's survival function, and the critical value it implies at level alpha. With the
# chi^2_0 atom carrying mass 1/2, P(2 log Lambda > t) = (1/2) P(chi^2_1 > t) for t > 0, so the
# level-alpha cut solves P(chi^2_1 > t) = 2 alpha.
t, alpha = sp.symbols("t alpha", positive=True)
chi1_survival = sp.erfc(sp.sqrt(t) / sp.sqrt(2))
# Asserted, not derived: the half-and-half mixture is a distributional result about where the
# unconstrained maximum falls relative to the boundary, which sympy cannot establish. It is
# Chernoff's, via Self and Liang. What follows from it *is* checkable, and is checked below.
mixture_survival = proof.define(r"P(2\log\Lambda > t)", chi1_survival / 2)

# Numerically that is 2.706 at alpha = 0.05, against 3.841 from the naive chi^2_1 cut. Using the
# naive cut halves the level; the ratio of the two is worth stating because it is the whole
# practical content of the boundary correction.
naive = sp.nsolve(chi1_survival - sp.Rational(5, 100), t, 3.0)
boundary = sp.nsolve(chi1_survival - sp.Rational(10, 100), t, 2.0)
proof.define(r"t_{0.05}^{\text{naive}}", sp.N(naive, 6))
proof.define(r"t_{0.05}^{\text{boundary}}", sp.N(boundary, 6))

# ------------------------------------------------------------------ where the power is

# Under a local alternative pi = 1 - eps, expand one observation's log density about the
# boundary. Writing r = f_0/f_1 - 1, the expansion is second order in eps with coefficient
# -r^2 / 2, so the statistic is quadratic in eps and its scale is set by the spread of r.
r = (f0 - f1) / f1
expansion = sp.log((1 - eps) * f1 + eps * f0) - sp.log(f1)
proof.claim(
    "the-local-expansion-about-the-boundary-is-quadratic-in-the-departure",
    sp.simplify(
        sp.series(expansion, eps, 0, 3).removeO() - (eps * r - eps**2 * r**2 / 2)
    ),
    r"\ell(1-\epsilon) - \ell(1) = \epsilon r - \tfrac12 \epsilon^2 r^2 + O(\epsilon^3)",
)

# Under the null every observation comes from f_1, so E[r] = int f_0 - int f_1 = 0: both are
# densities. The first-order term therefore has mean zero and the expected statistic is governed
# by E[r^2], which is the Fisher information in pi at the boundary. Asserted rather than checked,
# because it is an integral identity over unspecified densities.
proof.define(r"\mathbb{E}_{H_0}[r]", 0)
I_pipi = proof.define(r"I_{\pi\pi} = \mathbb{E}_{H_0}[r^2]", sp.Symbol("I_pipi", positive=True))
proof.define(r"\mathbb{E}[2\log\Lambda]", I_pipi * eps**2)

# Power is governed by I_pipi alone -- not by the Schur complement -- because mu is a nuisance
# maximised over rather than a second boundary parameter.

# From mixture_identification.py, one indicator contributes I_pipi = n_t v^2 / D with
# v = s_a - s_0, the gap between the two components' silence probabilities. That is *not* zero --
# the rank-1 result says the indicator cannot separate mu from pi, not that it says nothing about
# pi on its own. So tables do carry power for this test even though they carry none for the
# split.
v, D = sp.symbols("v D", positive=True)
indicator_I_pipi = n_tab * v**2 / D
proof.claim(
    "the-indicator-does-inform-pi-even-though-it-cannot-separate",
    sp.simplify(sp.diff(indicator_I_pipi, n_tab) - v**2 / D),
    r"I_{\pi\pi}^{\text{ind}} = n_t v^2 / D,\ \text{increasing in } n_t",
)

# And the gap v vanishes exactly when the two components are indistinguishable to a threshold:
# s_a(mu) = s_0 means a study with the effect is no more likely to report than one without.
proof.claim(
    "power-vanishes-when-the-two-components-report-alike",
    sp.simplify(indicator_I_pipi.subs(v, 0)),
    r"v \to 0 \implies I_{\pi\pi} \to 0",
)

# Which gap is that, in terms of a study's size? The null component's silence probability does not
# depend on the sample size at all: the cut in effect-size units is z / sqrt(n) and the null
# standard deviation is 1 / sqrt(n), so their ratio is z whatever n is.
z, n_subj = sp.symbols("z n", positive=True)
cut_in_g = z / sp.sqrt(n_subj)
sigma_null = 1 / sp.sqrt(n_subj)
proof.claim(
    "the-null-components-silence-is-free-of-the-sample-size",
    sp.simplify(cut_in_g / sigma_null - z),
    r"c/\sigma_0 = z \quad\text{for every } n",
)

# So the whole of v is carried by the active component, whose standard deviation sqrt(1/n + tau^2)
# falls towards tau while the cut falls to zero. A large study with the effect reports it; a small
# one is silent either way, and its silence says nothing about which component it came from.
tau = sp.Symbol("tau", positive=True)
sigma_active = sp.sqrt(1 / n_subj + tau**2)
proof.claim(
    "the-active-components-standard-deviation-tends-to-tau-not-to-zero",
    sp.simplify(sp.limit(sigma_active, n_subj, sp.oo) - tau),
    r"\sigma_a \to \tau,\quad c \to 0 \quad\text{as } n \to \infty",
)

# Numerically (checked in the calibration below) I_pipi runs 0.124 to 21.2 per table from n = 12
# to n = 200, and is convex in n over that range. Two consequences, and they differ from the
# mu/pi split, where only the *spread* of study designs matters and the level is irrelevant:
#
#   * the level dominates -- 400 tables at n = 20 carry the same information about pi as 15 at
#     n = 120;
#   * spread helps, but only by Jensen on a convex function: about a third, against the level's
#     two orders of magnitude.

def calibration():
    """I_pipi against sample size, and the Jensen comparison the text claims."""
    import numpy as np
    from scipy.stats import norm

    true_mu, tau_val, z_val = 0.5, 0.15, 3.09

    def information(n):
        sd_a, sd_0 = np.sqrt(1.0 / n + tau_val**2), np.sqrt(1.0 / n)
        cut = z_val / np.sqrt(n)
        s_active = norm.cdf((cut - true_mu) / sd_a) - norm.cdf((-cut - true_mu) / sd_a)
        s_zero = norm.cdf(cut / sd_0) - norm.cdf(-cut / sd_0)
        gap = s_active - s_zero
        return gap**2 / max(s_active * (1 - s_active), 1e-12)

    print()
    print("I_pipi per table against sample size:")
    for n in (12, 20, 40, 80, 120, 200):
        print(f"  n = {n:3d}: {information(n):8.4f}")

    grid = np.exp(np.linspace(np.log(12), np.log(120), 400))
    spread_value = np.mean([information(n) for n in grid])
    level_value = information(grid.mean())
    print(f"\nspread 12-120 at mean n {grid.mean():.1f}: E[I] {spread_value:.4f} against "
          f"I(E[n]) {level_value:.4f}, so spread helps by "
          f"{100 * (spread_value - level_value) / level_value:+.0f}%")
    print(f"level 12 to 200 multiplies I_pipi by "
          f"{information(200) / information(12):.0f}x, which is the larger lever")
    assert spread_value > level_value, "convexity claim failed"
    assert information(200) > 100 * information(12), "level claim failed"


if __name__ == "__main__":
    print(f"{len(proof.claims)} claims verified in {proof.name}")
    for label, shown in proof.claims:
        print(f"  {label}")
    print()
    print(f"critical values at alpha = 0.05: boundary {float(boundary):.4f}, "
          f"naive chi^2_1 {float(naive):.4f}")
    print("Using the naive cut runs the test at half its nominal level.")
    print()
    print("Unlike the mu/pi split, this test does draw power from coordinate tables: the")
    print("indicator's I_pipi is n_t v^2 / D, increasing in the number of tables. Rank 1 means")
    print("the indicator cannot separate mu from pi, not that it is silent about pi alone.")
    calibration()
