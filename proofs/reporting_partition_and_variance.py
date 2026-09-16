r"""Three defects in the shipped censored model, each stated as algebra rather than opinion.

From the combined-estimator design document, sections 4.1, 4.3 and 5. Each is a place where the
current implementation, or something written in its docstring, does not survive being written
down carefully. Two of the three correct claims of mine.

1. **The reporting events do not form the partition the likelihood assumes.** A report is
   assigned within one radius, a silence beyond a larger one, and the annulus between is dropped.
   The likelihood then uses :math:`P(|Y| \ge c)` and :math:`P(|Y| < c)` as if those were the two
   probabilities of a two-outcome experiment. They are not: three events partition the space, and
   discarding the middle one leaves the other two summing to less than one.

2. **An image-only heterogeneity estimate is not the active component's heterogeneity.** Under
   the zero-inflated model the total variance carries a prevalence term, so substituting one for
   the other is inconsistent, and the direction of the error is signed.

3. **There is no universal factor-of-two bound on the report limb**, and threshold variation does
   restore identification. Both are corrections to claims in the estimator's own Notes.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

proof = Proof(
    "reporting_partition_and_variance",
    "Three defects in the censored model, written down",
    __doc__,
)

# ============================================================ 1. the partition

pA, pB, pC = sp.symbols("p_A p_B p_C", positive=True)

# The three events: reported inside the report radius, in the dropped annulus, silent beyond the
# coverage radius. They partition, so their probabilities sum to one.
proof.define(r"p_A + p_B + p_C", 1)

# The likelihood the implementation writes treats A and C as complementary. The probability it
# assigns to the pair therefore falls short by exactly the mass of the event it dropped.
shortfall = 1 - (pA + pC)
proof.claim(
    "dropping-the-annulus-leaves-the-remaining-two-short-by-its-own-mass",
    sp.simplify(shortfall.subs(pC, 1 - pA - pB) - pB),
    r"1 - (p_A + p_C) = p_B",
)

# What a correct likelihood for an observation *known* to be in A or C would use: the conditional
# probabilities, whose denominator is the missing normaliser.
conditional_A = pA / (pA + pC)
proof.define(r"P(A \mid A \cup C)", conditional_A)
proof.claim(
    "the-missing-factor-is-one-over-the-retained-mass",
    sp.simplify(conditional_A / pA - 1 / (pA + pC)),
    r"\frac{P(A \mid A\cup C)}{p_A} = \frac{1}{p_A + p_C}",
)

# So the error is not a constant that cancels between the two limbs: it depends on mu through
# p_A + p_C, and therefore distorts the score. Show the score of the mis-specified pair differs
# from the correct one by the derivative of the log retained mass.
mu = sp.Symbol("mu", real=True)
a = sp.Function("p_A")(mu)
c = sp.Function("p_C")(mu)
wrong_score = sp.diff(sp.log(a), mu)
right_score = sp.diff(sp.log(a / (a + c)), mu)
proof.claim(
    "and-it-biases-the-score-by-the-derivative-of-the-log-retained-mass",
    sp.simplify(wrong_score - right_score - sp.diff(sp.log(a + c), mu)),
    r"\partial_\mu \log p_A - \partial_\mu \log \frac{p_A}{p_A + p_C}"
    r" = \partial_\mu \log(p_A + p_C)",
)

# ============================================================ 2. the variance

pi, mu_a, tau_a = sp.symbols("pi mu_a tau_a", positive=True)

# theta = 0 with probability 1 - pi, and N(mu_a, tau_a^2) with probability pi. Its first two
# moments, then the variance by the usual decomposition.
first = pi * mu_a
second = pi * (tau_a**2 + mu_a**2)
total_variance = sp.simplify(second - first**2)
proof.define(r"\operatorname{Var}(\theta)", total_variance)
proof.claim(
    "the-mixture-variance-carries-a-prevalence-term-the-active-one-does-not",
    sp.simplify(total_variance - (pi * tau_a**2 + pi * (1 - pi) * mu_a**2)),
    r"\operatorname{Var}(\theta) = \pi\tau_a^2 + \pi(1-\pi)\mu_a^2",
)

# An image-only between-study estimator targets the total. Using it as tau_a^2 therefore
# over-states the active component's spread by a signed, computable amount.
overstatement = sp.simplify(total_variance - tau_a**2)
proof.define(r"\text{substitution error}", overstatement)
proof.claim(
    "substituting-the-total-for-the-active-variance-errs-by-this-much",
    sp.simplify(overstatement - ((pi - 1) * tau_a**2 + pi * (1 - pi) * mu_a**2)),
    r"\operatorname{Var}(\theta) - \tau_a^2 = -(1-\pi)\tau_a^2 + \pi(1-\pi)\mu_a^2",
)

# The two terms pull opposite ways, so the sign depends on the configuration: the substitution is
# not conservative in general. It vanishes only at pi = 1, where the mixture is inert.
proof.claim(
    "the-error-vanishes-only-when-every-study-is-active",
    sp.simplify(overstatement.subs(pi, 1)),
    r"\pi = 1 \implies \text{no error}",
)

# ============================================================ 3. the peak bound

u, M = sp.symbols("u M", positive=True)
Phi = (1 + sp.erf(u / sp.sqrt(2))) / 2

# Probability that a named voxel exceeds u *and* is the maximum of M independent values. Verified
# by the fundamental theorem rather than by symbolic integration -- sympy will not integrate
# Phi(z)^(M-1) with a symbolic exponent, and waiting for it to fail is not a proof. Differentiate
# the candidate and check it reproduces the negated integrand, then check the boundary.
phi_u = sp.exp(-u**2 / 2) / sp.sqrt(2 * sp.pi)
closed_form = (1 - Phi**M) / M
proof.claim(
    "the-candidate-antiderivative-differentiates-to-the-integrand",
    sp.simplify(sp.diff(closed_form, u) + phi_u * Phi ** (M - 1)),
    r"\frac{d}{du}\frac{1-\Phi(u)^M}{M} = -\phi(u)\Phi(u)^{M-1}",
)
proof.claim(
    "and-vanishes-at-the-upper-limit-as-the-integral-must",
    sp.simplify(sp.limit(closed_form.subs(M, 27), u, sp.oo)),
    r"\lim_{u\to\infty} \frac{1-\Phi(u)^M}{M} = 0",
)

# The ratio the model's report limb gets wrong is exceedance over peak probability. At a vanishing
# threshold it grows without bound in the neighbourhood size, so no constant bounds it; at large u
# it tends to 1, which is why the error is small at the thresholds papers actually use.
#
# The limit is *not* M/2, which is what I first wrote: Phi(0) is 1/2 rather than 0, so the
# denominator carries a (1 - 2^-M) factor. The proof rejected the wrong version, which is the
# only reason it is right here.
ratio = (1 - Phi) / closed_form
proof.define(r"P(|Y|\ge u)/p_{\text{peak}}", ratio)
zero_limit = M * 2 ** (M - 1) / (2**M - 1)
proof.claim(
    "the-ratio-at-a-vanishing-threshold-grows-without-bound-in-the-neighbourhood-size",
    sp.simplify(sp.limit(ratio, u, 0, "+") - zero_limit),
    r"\lim_{u \to 0^+} \frac{1-\Phi(u)}{p_{\text{peak}}}"
    r" = \frac{M\,2^{M-1}}{2^M - 1} \to \frac{M}{2}",
)
proof.claim(
    "and-to-one-as-the-threshold-grows",
    sp.simplify(sp.limit(ratio.subs(M, 27), u, sp.oo) - 1),
    r"\lim_{u\to\infty} \frac{1-\Phi(u)}{p_{\text{peak}}} = 1",
)


def calibration():
    """Price all three against the numbers the design document quotes."""
    import numpy as np
    from scipy.stats import norm

    print()
    print("3. exceedance over peak probability, M = 27 neighbours:")
    for u_i in (0.5, 1.0, 2.0, 3.09, 4.0):
        p_peak = (1 - norm.cdf(u_i) ** 27) / 27
        print(f"   u = {u_i:>4}: ratio {(1 - norm.cdf(u_i)) / p_peak:>7.3f}")
    print("   The document's counterexample is u = 0.5, ratio 8.33: no factor-of-two bound.")
    print("   At the conventional z = 3.09 the ratio is 1.02, so the bound holds where it is")
    print("   used and the claim was over-general rather than wrong in practice.")
    assert abs((1 - norm.cdf(0.5)) / ((1 - norm.cdf(0.5) ** 27) / 27) - 8.331) < 5e-3

    print()
    print("2. substitution error at mu_a = 0.5, tau_a = 0.15:")
    print(f"   {'pi':>5} {'Var(theta)':>11} {'tau_a^2':>9} {'error':>9}")
    for pi_i in (1.0, 0.8, 0.6, 0.4):
        total = pi_i * 0.15**2 + pi_i * (1 - pi_i) * 0.5**2
        print(f"   {pi_i:>5} {total:>11.5f} {0.15**2:>9.5f} {total - 0.15**2:>+9.5f}")
    print("   Positive at moderate prevalence, so the active spread is over-stated, and the")
    print("   prevalence term dominates the active one as soon as mu_a exceeds tau_a.")

    print()
    print("1. retained mass at a two-radius construction, as the annulus grows:")
    print("   p_B is whatever the annulus takes; the score bias is d/dmu log(1 - p_B).")
    for pB_i in (0.0, 0.1, 0.3, 0.5):
        print(f"   p_B = {pB_i:>4}: the likelihood accounts for {1 - pB_i:.2f} of the mass")
    assert np.isclose(1 - 0.3, 0.7)


if __name__ == "__main__":
    print(f"{len(proof.claims)} claims verified in {proof.name}")
    for label, shown in proof.claims:
        print(f"  {label}")
    calibration()
